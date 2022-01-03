/**
 * Oct 30, 2017
 */
package siCloneFiT.algorithm;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;

import SiFit.model.ComplexEvolutionModel;
import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import SiCloneFiT.likelihood.CompleteLikelihood;
import SiCloneFiT.objects.Clone;
import SiCloneFiT.objects.DoubletPosteriorSample;
import SiCloneFiT.objects.SingleCell;
import SiCloneFiT.posterior.CompleteDoubletPosterior;
import SiCloneFiT.proposal.TreeProposal;
import SiCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Oct 30, 2017
 */
public class GibbsSamplingDoublet extends GibbsSamplingSinglet {
	
	public static DoubletPosteriorSample samplePosterior(int nCell, int nMut, ArrayList<SingleCell> listSingleCells,
														 ArrayList<String> singleCellNames, HashMap<Integer, Integer[]> cellIDGTVectorMap,
														 double fn, BetaDistribution fnPriorDist, double fp, BetaDistribution fpPriorDist,
														 double _doublet, BetaDistribution doubletPriorDist,
														 double _delProb, BetaDistribution delProbPriorDist, 
														 double _LOHProb, BetaDistribution LOHProbPriorDist,
														 double _recurProb, BetaDistribution recurProbPriorDist,
														 double alpha_0, double GammaA, double GammaB, GammaDistribution alpha_0PriorDist,
														 int burnInIter, int iterT, int iterP, int sampleTreeIter, int dataFlag, int restartIndex,
														 String outDir, SCFUtilityFunctions SCF) throws FileNotFoundException, UnsupportedEncodingException{
		// Initiate Clone objects
		int startClone = 20; // Number of clones to start with
		ArrayList<Clone> listClone = new ArrayList<>();
		for (int i = 0; i < startClone; i++){
			Clone c = new Clone(i);
			c.setNameID(i);
			listClone.add(c);
		}	
		
		// Assign cells to Clones, random assignments
		for (SingleCell S : listSingleCells){
			int cellCloneID = _rng.nextInt((19 - 0) + 1) + 0;
			S.cloneID = cellCloneID;
			listClone.get(cellCloneID).memberCellList.add(S.cellID);
		}
		
//		ArrayList<Integer> cellClonesbfr = new ArrayList<>();
//        for (SingleCell s: listSingleCells){
//        	cellClonesbfr.add(s.cloneID);
//        }
//        System.out.println("begin");
//        System.out.println(cellClonesbfr);
        
		// Remove empty clones
		listClone = SCF.removeEmptyClones(listClone);
		
		ArrayList<String> cloneNames = new ArrayList<>();
		
		// Assign Genotypes to Clones, Randomly now
		for (Clone C: listClone){
			C.assignRandomGTVector(nMut, listSingleCells, SCF);
			cloneNames.add(C.cloneName);
		}
		
		// Vector of flags for cells indicating cell is doublet or not
		int[] cellDoubletFlagArr = new int[nCell]; // The vector Y, Y_j denotes if cell j is a doublet or not
		for (int j = 0; j < nCell; j++){
			double rr = _rng.nextDouble();
			if (rr <= _doublet)
				cellDoubletFlagArr[j] = 1;
			else
				cellDoubletFlagArr[j] = 0;
		}
		
		// Assign 2nd clone to cells randomly
		ArrayList<Integer> c2CellCloneList = new ArrayList<>();
		for (int j = 0; j < nCell; j++){
			SingleCell s_j = listSingleCells.get(j);
			int s_j_clone2Id = SCF.assignRandomClone2(s_j, listClone); 
			c2CellCloneList.add(s_j_clone2Id);			
		}
		
		// Obtain the genotypeObservation for clones
		ArrayList<GenotypeObservation> cloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(listClone, cloneNames, nMut);

		double currFn = fn;
		double currFp = fp;
		double currDoublet = _doublet;
		
        // Generate random tree
        STITree<Clone> randTreeMaker = SCF.generateRandomCloneTree(cloneNames);
        STITree<Clone> randTree = SCF.getCloneTree(randTreeMaker.toNewick());       
        STITree<Clone> currTree = randTree;
        
        // Evolution Model
        ComplexEvolutionModel m = new ComplexEvolutionModel(_delProb, _LOHProb, _recurProb);
        
        // Object for proposing new tree 
        TreeProposal proposeTreeObj = new TreeProposal(); 

        // Complete Likelihood object (error likelihood + genotype likelihood of tree)
        CompleteLikelihood CLObj = new CompleteLikelihood(currTree, m, currFp, currFn, dataFlag);
        
        
        // Run heuristic search 
        DoubletPosteriorSample bestMLsample = getBestMLSample(nCell, nMut, listSingleCells, singleCellNames, cellIDGTVectorMap, listClone, currTree, cloneNames, cloneGenotypeObsMatrix, c2CellCloneList, cellDoubletFlagArr, CLObj, proposeTreeObj, m, currFn, fnPriorDist, currFp, fpPriorDist, currDoublet, doubletPriorDist, _delProb, delProbPriorDist, _LOHProb, LOHProbPriorDist, _recurProb, recurProbPriorDist, alpha_0, GammaA, GammaB, alpha_0PriorDist, 100, 20, 200, dataFlag, SCF);
        listClone = bestMLsample.listClone;
        c2CellCloneList = bestMLsample.cell2ndCloneIDVector;
        cellDoubletFlagArr = bestMLsample.cellDoubletFlagArr;
        currTree = bestMLsample.cloneTree;
        currFn = bestMLsample.sampleFn;
        currFp = bestMLsample.sampleFp;
        currDoublet = bestMLsample.sampleDoublet;
        
        ArrayList<String> copyCloneNames = SCF.getCloneNames(listClone);
        cloneNames = copyCloneNames;
    	ArrayList<GenotypeObservation> copyCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(listClone, copyCloneNames, nMut);
    	cloneGenotypeObsMatrix = copyCloneGenotypeObsMatrix;
    	
    	for (Clone C: listClone){
    		for (int cellID : C.memberCellList){
    			SingleCell s = listSingleCells.get(cellID);
    			s.cloneID = C.cloneID;
    			s.cloneName = C.cloneName;
    			s.cloneGTVector = C.cloneGTVector;
    		}
    	}
        
//        ArrayList<Integer> cellClones = new ArrayList<>();
//        for (SingleCell s: listSingleCells){
//        	cellClones.add(s.cloneID);
//        }
//        System.out.println("after HS");
//        System.out.println(cellClones);
        
        // Posterior probability of current sample 
        double posterior = CompleteDoubletPosterior.computeFullModelLogPosterior(listSingleCells, cellDoubletFlagArr, listClone, c2CellCloneList, cloneGenotypeObsMatrix, CLObj, proposeTreeObj, currTree, m, nMut, nCell, dataFlag, currFn, fnPriorDist, currFp, fpPriorDist, currDoublet, doubletPriorDist, _delProb, delProbPriorDist, _LOHProb, LOHProbPriorDist, alpha_0, alpha_0PriorDist, SCF);
        System.out.println("start sample posterior = "+ posterior);
        double likelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listSingleCells, listClone, c2CellCloneList, cellDoubletFlagArr, currFp, currFn, nMut, SCF, dataFlag);
//        double likelihood = CLObj.computeFullLikelihoodDoublet(listSingleCells, cellDoubletFlagArr, listClone, c2CellCloneList, cloneGenotypeObsMatrix, currTree, m, SCF, currFp, currFn, nMut, dataFlag);
        System.out.println("start sample likelihood = "+ likelihood);
        // Determine the frequency of sampling
        int sampleRate;
        if (iterT < 1000)
        	sampleRate = 10;
        else
        	sampleRate = iterT/200;
        
        // Best configurations and likelihoods
        ArrayList<Clone> bestClone = listClone;
        ArrayList<Integer> best2ndClone = c2CellCloneList;
        int[] bestCellDoubletFlag = cellDoubletFlagArr;
        double bestL = Double.NEGATIVE_INFINITY;
        double bestP = Double.NEGATIVE_INFINITY;
        STITree<Clone> bestTree = currTree;
        double bestFn = currFn;
        double bestFp = currFp;
        double bestDoublet = currDoublet;
		
        // Files for writing the samples
        String clonalClusterFile = outDir + "samples_" + restartIndex + "_clusters_1.txt";
        String clonalCluster2File = outDir + "samples_" + restartIndex + "_clusters_2.txt";
        String clonalGenotypesFile = outDir + "samples_" + restartIndex + "_genotypes.txt";
        String clonalTreesFile = outDir + "samples_" + restartIndex + "_trees.txt";
        String cellDoubleFlagFile = outDir + "samples_" + restartIndex + "_doubletFlags.txt";
        String posteriorErrorRatesFile = outDir + "samples_" + restartIndex + "_posteriors.txt";
        String predictedFile = outDir + "sample_" + restartIndex + "_MAP_predicted_genotype.txt";
        
        PrintWriter clusterWriter = new PrintWriter(clonalClusterFile, "UTF-8");
        PrintWriter cluster2Writer = new PrintWriter(clonalCluster2File, "UTF-8");
        PrintWriter genotypeWriter = new PrintWriter(clonalGenotypesFile, "UTF-8");
        PrintWriter treeWriter = new PrintWriter(clonalTreesFile, "UTF-8");
        PrintWriter doubletFlagWriter = new PrintWriter(cellDoubleFlagFile, "UTF-8");
        PrintWriter posteriorWriter = new PrintWriter(posteriorErrorRatesFile, "UTF-8");

        int totalIter = burnInIter + iterT;
        // For loop over the number of iterations
        for (int i = 1; i < totalIter; i++){
        	/***********************************
        	 * Update main indicators of cells *
        	 ***********************************/
        	// First, we are going to sample the main indicators
        	// for all the cells. Loop over nCell
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = listSingleCells.get(j);
        		Clone C_j = SCF.getCellClone(listClone, s_j);
        		
        		// The clone is not singleton, new clone is to be proposed
        		if (C_j.isSingleton() == false){
        			double acRatio = SampleIndicators.proposeNewCloneNonSingletonCellDoublet(listClone, cloneNames, s_j, C_j, c2CellCloneList, cellDoubletFlagArr, j, currTree, proposeTreeObj, cloneGenotypeObsMatrix, m, nMut, nCell, dataFlag, CLObj, currFp, currFn, SCF, alpha_0);
        			double rr1 = _rng.nextDouble();
        			if (rr1 <= acRatio){
        				listClone = SampleIndicators.nuCloneAddedCloneList;
        				currTree = SampleIndicators.nuCloneAddedTree;
        				cloneNames = SampleIndicators.nuCloneAddedCloneNameList;
        				cloneGenotypeObsMatrix = SampleIndicators.nuCloneAddedCloneGenotypeObsMatrix;
        				SCF.updateSingleCellList(listSingleCells, SampleIndicators.nuClone, s_j);
        			}
        		}
        		// The clone is a singleton, assign the cell to another clone
        		else{
        			double acRatio = SampleIndicators.removeCloneSingletonCellDoublet(listClone, cloneNames, s_j, C_j, c2CellCloneList, cellDoubletFlagArr, j, currTree, proposeTreeObj, cloneGenotypeObsMatrix, m, nMut, nCell, dataFlag, CLObj, currFp, currFn, SCF, alpha_0);
        			double rr1 = _rng.nextDouble();
        			if (rr1 <= acRatio){
        				listClone = SampleIndicators.cloneRemovedCloneList;
        				SCF.updateCloneList(listClone, SampleIndicators.nuCloneCellj, j);
        				currTree = SampleIndicators.cloneRemovedTree;
        				cloneNames = SampleIndicators.cloneRemovedCloneNameList;
        				cloneGenotypeObsMatrix = SampleIndicators.cloneRemovedCloneGenotypeObsMatrix;
        				c2CellCloneList = SampleIndicators.cloneRemovedC2CellCloneList;
        				SCF.updateSingleCellList(listSingleCells, listClone);
        			}
        		}
        	}
        	
        	// Perform the final steps of sampling indicator c1
        	// Here No change is made for c2, Maybe TODO: if c1_j == c2_j, should we modify that?
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = listSingleCells.get(j);        		
        		Clone C_j = SCF.getCellClone(listClone, s_j);
        		if (C_j.isSingleton() == false){
        			ArrayList<Double> unnormalizedDist = new ArrayList<>();
        			for (Clone C: listClone){
        				Clone doubletMateClone = SCF.getCloneFrmList(listClone, c2CellCloneList.get(j));
        				double val = PartialMHSteps.computeCRPErrorLikelihoodDoublet(C, s_j, doubletMateClone, cellDoubletFlagArr[j], nMut, nCell, currFp, currFn, dataFlag, CLObj.ampErrLikelihoodObj);
        				unnormalizedDist.add(val);
        			}
        			Clone C_j_new = SampleCloneGenotype.sampleNewCloneCRPErrLikelihoodDist(unnormalizedDist, listClone);
        			C_j_new.memberCellList.add(s_j.cellID);
        			C_j.removeCell(s_j.cellID);
        			SCF.updateSingleCellList(listSingleCells, C_j_new, s_j);
        		}
        	}
        	
        	/*************************************
        	 * Update second indicators of cells *
        	 *************************************/
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = listSingleCells.get(j);   
        		Clone C_j = SCF.getCellClone(listClone, s_j);
        		Clone s_j_c2 = SampleIndicator2Doublet.sampleDoubletMateClone(s_j, C_j, listClone, currFp, currFn, nMut, dataFlag, CLObj.ampErrLikelihoodObj, SCF);        		
        		c2CellCloneList.set(j, s_j_c2.cloneID);
        	}
        	
        	/**************************************
        	 * Update Tree and Model of Evolution *
        	 **************************************/
        	// Update the tree by sampling using metropolis-hastings
        	if (currTree.getLeafCount() > 3){
//        		int sampleTreeIter = 100;
        		STITree<Clone> newTree = SampleTreeModel.sampleCloneTreeMH(currTree, CLObj, _delProb, _LOHProb, sampleTreeIter, delProbPriorDist, LOHProbPriorDist, cloneGenotypeObsMatrix, nMut, dataFlag, proposeTreeObj, listClone.size());
        		currTree = newTree;        		     
        	}
        	
        	/***************************
        	 * Update Clonal Genotypes *
        	 ***************************/
        	// Sample Clonal Genotypes
        	for (int k = 0; k < listClone.size(); k++){
    			Clone C_k = listClone.get(k);
    			ArrayList<Clone> otherCloneList = SCF.getOtherCloneList(listClone, k, nMut);
    			ArrayList<String> otherCloneNames = new ArrayList<>();
    			for (Clone C: otherCloneList)
    				otherCloneNames.add(C.cloneName);
    			ArrayList<GenotypeObservation> otherClonesGenotypeObs = SCF.getCloneGenotypeObs(otherCloneList, otherCloneNames, nMut);
    			Integer[] C_k_newGTVector = SampleCloneGenotype.sampleCloneGenotypeVectorFromTreeAndErrorLikelihoodDistDoublet(otherClonesGenotypeObs, C_k, listClone, c2CellCloneList, cellDoubletFlagArr, currTree, m, cellIDGTVectorMap, CLObj.ampErrLikelihoodObj, dataFlag, nMut, SCF);
    			C_k.cloneGTVector = C_k_newGTVector;
        	}
        	
        	// Update clone gebotype matrix
        	cloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(listClone, cloneNames, nMut);
        	
        	/**********************
        	 * Update Error Rates *
        	 **********************/
        	
        	/*
        	 * Update beta, i.e., value of Fn
        	 */
        	double newFn = SampleErrorRates.sampleNewFnDoublet(listSingleCells, listClone, c2CellCloneList, cellDoubletFlagArr, CLObj.ampErrLikelihoodObj, currFp, fnPriorDist, nMut, SCF, dataFlag);
        	CLObj.ampErrLikelihoodObj.setFn(newFn, dataFlag);
        	currFn = newFn;
        	
        	/*
        	 * Update alpha, i.e., value of Fp
        	 */
        	double newFp = SampleErrorRates.sampleNewFpDoublet(listSingleCells, listClone, c2CellCloneList, cellDoubletFlagArr, CLObj.ampErrLikelihoodObj, fpPriorDist, currFn, nMut, SCF, dataFlag);
        	CLObj.ampErrLikelihoodObj.setFp(newFp, dataFlag);
        	currFp = newFp;
        	
        	/*
        	 * Update delta, i.e, value of _doublet
        	 */
        	double newDelta = SampleErrorRates.sampleNewDeltaDoublet(cellDoubletFlagArr, doubletPriorDist, nCell);
        	
        	/*
        	 * Update the Doublet Indicator Vector 
        	 */
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = listSingleCells.get(j);
        		Clone C_j = SCF.getCellClone(listClone, s_j);
        		Clone doubletMateClone_j = SCF.getCloneFrmList(listClone, c2CellCloneList.get(j));
        		int s_doubletFlag = SampleCellDoubletIndicator.sampleCellDoubletIndicator(s_j, C_j, doubletMateClone_j, currFp, currFn, newDelta, CLObj.ampErrLikelihoodObj, dataFlag, nMut);
        		cellDoubletFlagArr[j] = s_doubletFlag;
        	}
        	
        	posterior = CompleteDoubletPosterior.computeFullModelLogPosterior(listSingleCells, cellDoubletFlagArr, listClone, c2CellCloneList, cloneGenotypeObsMatrix, CLObj, proposeTreeObj, currTree, m, nMut, nCell, dataFlag, currFn, fnPriorDist, currFp, fpPriorDist, newDelta, doubletPriorDist, _delProb, delProbPriorDist, _LOHProb, LOHProbPriorDist, alpha_0, alpha_0PriorDist, SCF);
        	likelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listSingleCells, listClone, c2CellCloneList, cellDoubletFlagArr, newFp, newFn, nMut, SCF, dataFlag);
//        	likelihood = CLObj.computeFullLikelihoodDoublet(listSingleCells, cellDoubletFlagArr, listClone, c2CellCloneList, cloneGenotypeObsMatrix, currTree, m, SCF, newFp, newFn, nMut, dataFlag);
        	
        	// Update the MAP solution
        	if (posterior > bestP){
        		bestL = likelihood;
        		bestP = posterior;
        		bestClone = listClone;
        		best2ndClone = c2CellCloneList;
        		bestCellDoubletFlag = cellDoubletFlagArr;
        		bestTree = currTree;
        		bestFn = currFn;
        		bestFp = currFp;
        		bestDoublet = newDelta;
        		
        	}
        	
        	alpha_0 = SCF.sampleConcentrationParam(alpha_0, listClone.size(), nCell, GammaA, GammaB);
        	
        	// Print occasionally
        	if (i % iterP == 0){
        		System.out.println("posterior = " + posterior + " " + "likelihood = "+ likelihood);
        	}
        	
        	// Write the samples to files
        	if (i > burnInIter & i % sampleRate == 0){
        		// Write the clone index of cells
        		String cellCloneIDStr = getCellCloneString(listClone, nCell);
        		clusterWriter.printf("%s\n", cellCloneIDStr);
        		// Write the 2nd clone index of cells
        		String cell2ndCloneIDStr = getSecondCellCloneString(c2CellCloneList);
        		cluster2Writer.printf("%s\n", cell2ndCloneIDStr);
        		// Write cell doublet flags
        		String cellDoubletFlagStr = getCellDoubletFlagString(cellDoubletFlagArr);
        		doubletFlagWriter.printf("%s\n", cellDoubletFlagStr);
        		// Write the clonal tree
        		String cellClonalTree = getCellClonalTreeString(listClone, nCell, singleCellNames, currTree);
        		treeWriter.printf("%s\n", cellClonalTree);
        		// Write posterior probability of sample and error rates
        		ArrayList<String> posteriorsNerrorRates = new ArrayList<>();
        		posteriorsNerrorRates.add(Double.toString(posterior));
        		posteriorsNerrorRates.add(Double.toString(likelihood));
        		posteriorsNerrorRates.add(Double.toString(currFp));
        		posteriorsNerrorRates.add(Double.toString(currFn));
        		posteriorsNerrorRates.add(Double.toString(newDelta));
        		String posteriorsNerrorRatesStr = String.join(" ", posteriorsNerrorRates);
        		posteriorWriter.printf("%s\n", posteriorsNerrorRatesStr);
        		// Write the clonal genotypes
        		String cloneGenotypeListStr = getCloneGenotypeListString(listClone);
        		genotypeWriter.printf("%s\n", cloneGenotypeListStr);
        	}
        }
        
        // Close the file writers
        clusterWriter.close();
        cluster2Writer.close();
        treeWriter.close();
        genotypeWriter.close();
        posteriorWriter.close();
        doubletFlagWriter.close();
		
        System.out.println("MAP value = " + bestP);
        HashMap<Integer, Integer[]> cellPredictedGt = new HashMap<>();
        int[] cellCloneID = new int[nCell];
        for (Clone C : bestClone){
			System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
			System.out.println(C.memberCellList);
			System.out.println(Arrays.toString(C.cloneGTVector));
			ArrayList<String> thisCloneCells = new ArrayList<>();
			
			for (int cell: C.memberCellList){
				thisCloneCells.add(singleCellNames.get(cell));
				cellPredictedGt.put(cell, C.cloneGTVector);
				cellCloneID[cell] = C.cloneID;
			}
			System.out.println(thisCloneCells);
			STINode<Clone> C_node = bestTree.getNode(C.cloneName);
			for (String cell: thisCloneCells){
				C_node.createChild(cell);
				bestTree.getNode(cell).setParentDistance(0.000000000000001);
			}
		}
        System.out.println(bestTree.toNewick());
        
        System.out.println(Arrays.toString(cellCloneID));


        Integer[][] predictedGenotypeMat = new Integer[nMut][nCell+1];
        SCF.writeGenotypeMatrix2File(predictedFile, nCell, cellPredictedGt, predictedGenotypeMat, nMut);
        
        DoubletPosteriorSample MAPSample = new DoubletPosteriorSample(bestClone, cellCloneID, best2ndClone, bestTree, bestCellDoubletFlag, bestP, bestL, bestFp, bestFn, bestDoublet);
		
		return MAPSample;
	}
	
	
	/**
	 * Run heuristic search to get best ML sample
	 * @param nCell
	 * @param nMut
	 * @param listSingleCells
	 * @param singleCellNames
	 * @param cellIDGTVectorMap
	 * @param listClone
	 * @param currTree
	 * @param cloneNames
	 * @param cloneGenotypeObsMatrix
	 * @param c2CellCloneList
	 * @param cellDoubletFlagArr
	 * @param CLObj
	 * @param proposeTreeObj
	 * @param m
	 * @param currFn
	 * @param fnPriorDist
	 * @param currFp
	 * @param fpPriorDist
	 * @param currDoublet
	 * @param doubletPriorDist
	 * @param _delProb
	 * @param delProbPriorDist
	 * @param _LOHProb
	 * @param LOHProbPriorDist
	 * @param _recurProb
	 * @param recurProbPriorDist
	 * @param alpha_0
	 * @param GammaA
	 * @param GammaB
	 * @param alpha_0PriorDist
	 * @param iterT
	 * @param iterP
	 * @param sampleTreeIter
	 * @param dataFlag
	 * @param SCF
	 * @return
	 * Created On: Dec 21, 2017
	 */
	public static DoubletPosteriorSample getBestMLSample(int nCell, int nMut, ArrayList<SingleCell> listSingleCells,
														 ArrayList<String> singleCellNames, HashMap<Integer, Integer[]> cellIDGTVectorMap,
														 ArrayList<Clone> listClone, STITree<Clone> currTree,
														 ArrayList<String> cloneNames, ArrayList<GenotypeObservation> cloneGenotypeObsMatrix,
														 ArrayList<Integer> c2CellCloneList, int[] cellDoubletFlagArr, CompleteLikelihood CLObj,
														 TreeProposal proposeTreeObj, ComplexEvolutionModel m,
														 double currFn, BetaDistribution fnPriorDist, double currFp, BetaDistribution fpPriorDist,
														 double currDoublet, BetaDistribution doubletPriorDist,
														 double _delProb, BetaDistribution delProbPriorDist, 
														 double _LOHProb, BetaDistribution LOHProbPriorDist,
														 double _recurProb, BetaDistribution recurProbPriorDist,
														 double alpha_0, double GammaA, double GammaB, GammaDistribution alpha_0PriorDist,
														 int iterT, int iterP, int sampleTreeIter, int dataFlag, SCFUtilityFunctions SCF){
		
		// Best configurations and likelihoods
        ArrayList<Clone> bestClone = listClone;
        ArrayList<Integer> best2ndClone = c2CellCloneList;
        double bestL = Double.NEGATIVE_INFINITY;
        STITree<Clone> bestTree = currTree;
        double bestFn = currFn;
        double bestFp = currFp;
        double bestDoublet = currDoublet;
        int[] bestDoubletFlagArr = new int[nCell];
        
        double currLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listSingleCells, listClone, c2CellCloneList, cellDoubletFlagArr, currFp, currFn, nMut, SCF, dataFlag);
        double nuLikelihood;
        
        // For loop over the number of iterations
        for (int i = 1; i < iterT; i++){
        	// Create copies of all important variables
        	ArrayList<Clone> copyListClone = SCF.copyCloneList(listClone, nMut);
        	ArrayList<SingleCell> copyListSingleCells = SCF.copyCellList(listSingleCells, nMut);
        	STITree<Clone> copyCurrTree = new STITree<>(currTree);
        	ArrayList<String> copyCloneNames = SCF.getCloneNames(copyListClone);
        	ArrayList<GenotypeObservation> copyCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(copyListClone, copyCloneNames, nMut);
        	ArrayList<Integer> copyC2CellCloneList = SCF.getCopyC2CellCloneList(c2CellCloneList);
        	
        	/***********************************
        	 * Update main indicators of cells *
        	 ***********************************/       	
        	// First, we are going to sample the main indicators
        	// for all the cells. Loop over nCell
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = copyListSingleCells.get(j);
        		Clone C_j = SCF.getCellClone(copyListClone, s_j);
        		
        		// The clone is not singleton, new clone is to be proposed
        		if (C_j.isSingleton() == false){
        			double acRatio = SampleIndicators.proposeNewCloneNonSingletonCellDoublet(copyListClone, copyCloneNames, s_j, C_j, copyC2CellCloneList, cellDoubletFlagArr, j, copyCurrTree, proposeTreeObj, copyCloneGenotypeObsMatrix, m, nMut, nCell, dataFlag, CLObj, currFp, currFn, SCF, alpha_0);
        			double rr1 = _rng.nextDouble();
        			if (rr1 <= acRatio){
        				copyListClone = SampleIndicators.nuCloneAddedCloneList;
        				copyCurrTree = SampleIndicators.nuCloneAddedTree;
        				copyCloneNames = SampleIndicators.nuCloneAddedCloneNameList;
        				copyCloneGenotypeObsMatrix = SampleIndicators.nuCloneAddedCloneGenotypeObsMatrix;
        				SCF.updateSingleCellList(copyListSingleCells, SampleIndicators.nuClone, s_j);
        			}
        		}
        		// The clone is a singleton, assign the cell to another clone
        		else{
        			double acRatio = SampleIndicators.removeCloneSingletonCellDoublet(copyListClone, copyCloneNames, s_j, C_j, copyC2CellCloneList, cellDoubletFlagArr, j, copyCurrTree, proposeTreeObj, copyCloneGenotypeObsMatrix, m, nMut, nCell, dataFlag, CLObj, currFp, currFn, SCF, alpha_0);
        			double rr1 = _rng.nextDouble();
        			if (rr1 <= acRatio){
        				copyListClone = SampleIndicators.cloneRemovedCloneList;
        				SCF.updateCloneList(copyListClone, SampleIndicators.nuCloneCellj, j);
        				copyCurrTree = SampleIndicators.cloneRemovedTree;
        				copyCloneNames = SampleIndicators.cloneRemovedCloneNameList;
        				copyCloneGenotypeObsMatrix = SampleIndicators.cloneRemovedCloneGenotypeObsMatrix;
        				copyC2CellCloneList = SampleIndicators.cloneRemovedC2CellCloneList;
        				SCF.updateSingleCellList(copyListSingleCells, copyListClone);
        			}
        		}
        		
        	}
        	
        	// Perform the final steps of sampling indicator c1
        	// Here No change is made for c2, Maybe TODO: if c1_j == c2_j, should we modify that?
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = copyListSingleCells.get(j);        		
        		Clone C_j = SCF.getCellClone(copyListClone, s_j);
        		if (C_j.isSingleton() == false){
        			ArrayList<Double> unnormalizedDist = new ArrayList<>();
        			for (Clone C: copyListClone){
        				Clone doubletMateClone = SCF.getCloneFrmList(copyListClone, copyC2CellCloneList.get(j));
        				double val = PartialMHSteps.computeCRPErrorLikelihoodDoublet(C, s_j, doubletMateClone, cellDoubletFlagArr[j], nMut, nCell, currFp, currFn, dataFlag, CLObj.ampErrLikelihoodObj);
        				unnormalizedDist.add(val);
        			}
        			Clone C_j_new = SampleCloneGenotype.sampleNewCloneCRPErrLikelihoodDist(unnormalizedDist, copyListClone);
        			C_j_new.memberCellList.add(s_j.cellID);
        			C_j.removeCell(s_j.cellID);
        			SCF.updateSingleCellList(copyListSingleCells, C_j_new, s_j);
        		}
        	}
        	
        	nuLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(copyListSingleCells, copyListClone, copyC2CellCloneList, cellDoubletFlagArr, currFp, currFn, nMut, SCF, dataFlag);
        	
        	// better likelihood
        	// Update all variables
        	if (nuLikelihood > currLikelihood){
        		currLikelihood = nuLikelihood;
        		listClone = copyListClone;
        		listSingleCells = copyListSingleCells;
        		cloneNames = copyCloneNames;
        		cloneGenotypeObsMatrix = copyCloneGenotypeObsMatrix;
        		c2CellCloneList = copyC2CellCloneList;
        		currTree = copyCurrTree;
        	}
        	
        	/*************************************
        	 * Update second indicators of cells *
        	 *************************************/
        	// Sample indicator 2, c2 for each cell
        	copyC2CellCloneList = SCF.getCopyC2CellCloneList(c2CellCloneList);
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = listSingleCells.get(j);   
        		Clone C_j = SCF.getCellClone(listClone, s_j);
        		Clone s_j_c2 = SampleIndicator2Doublet.sampleDoubletMateClone(s_j, C_j, listClone, currFp, currFn, nMut, dataFlag, CLObj.ampErrLikelihoodObj, SCF);        		
        		copyC2CellCloneList.set(j, s_j_c2.cloneID);
        	}
        	
        	nuLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listSingleCells, listClone, copyC2CellCloneList, cellDoubletFlagArr, currFp, currFn, nMut, SCF, dataFlag);
        	// Better likelihood, update corresponding variable
        	if (nuLikelihood > currLikelihood){
        		c2CellCloneList = copyC2CellCloneList;
        	}
        	
        	/**************************************
        	 * Update Tree and Model of Evolution *
        	 **************************************/
        	// Update the tree by sampling using metropolis-hastings
        	if (currTree.getLeafCount() > 2){
        		STITree<Clone> newTree = SampleTreeModel.sampleCloneTreeMH(currTree, CLObj, _delProb, _LOHProb, sampleTreeIter, delProbPriorDist, LOHProbPriorDist, cloneGenotypeObsMatrix, nMut, dataFlag, proposeTreeObj, listClone.size());
        		copyCurrTree = newTree;        		     
        	}
        	else{
        		copyCurrTree = currTree;
        	}
        	
        	/***************************
        	 * Update Clonal Genotypes *
        	 ***************************/
        	// Create copies of all important variables
        	copyListClone = SCF.copyCloneList(listClone, nMut);
        	
        	// Sample Clonal Genotypes
        	for (int k = 0; k < copyListClone.size(); k++){
    			Clone C_k = copyListClone.get(k);
    			ArrayList<Clone> otherCloneList = SCF.getOtherCloneList(copyListClone, k, nMut);
    			ArrayList<String> otherCloneNames = new ArrayList<>();
    			for (Clone C: otherCloneList)
    				otherCloneNames.add(C.cloneName);
    			ArrayList<GenotypeObservation> otherClonesGenotypeObs = SCF.getCloneGenotypeObs(otherCloneList, otherCloneNames, nMut);
    			Integer[] C_k_newGTVector = SampleCloneGenotype.sampleCloneGenotypeVectorFromTreeAndErrorLikelihoodDistDoublet(otherClonesGenotypeObs, C_k, copyListClone, c2CellCloneList, cellDoubletFlagArr, copyCurrTree, m, cellIDGTVectorMap, CLObj.ampErrLikelihoodObj, dataFlag, nMut, SCF);
    			C_k.cloneGTVector = C_k_newGTVector;
        	}
        	
        	copyCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(copyListClone, cloneNames, nMut);
        	nuLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listSingleCells, copyListClone, c2CellCloneList, cellDoubletFlagArr, currFp, currFn, nMut, SCF, dataFlag);

        	// better likelihood
        	// Update all variables
        	if (nuLikelihood > currLikelihood){
        		currLikelihood = nuLikelihood;
        		listClone = copyListClone;
        		cloneGenotypeObsMatrix = copyCloneGenotypeObsMatrix;
        		currTree = copyCurrTree;
        	}

        	/**********************
        	 * Update Error Rates *
        	 **********************/
        	
        	/*
        	 * Update beta, i.e., value of Fn
        	 */
        	double newFn = SampleErrorRates.sampleNewFnDoublet(listSingleCells, listClone, c2CellCloneList, cellDoubletFlagArr, CLObj.ampErrLikelihoodObj, currFp, fnPriorDist, nMut, SCF, dataFlag);
        	nuLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listSingleCells, listClone, c2CellCloneList, cellDoubletFlagArr, currFp, newFn, nMut, SCF, dataFlag);
        	// better likelihood
        	// Update all variables
        	if (nuLikelihood > currLikelihood){
        		currLikelihood = nuLikelihood;
        		CLObj.ampErrLikelihoodObj.setFn(newFn, dataFlag);
            	currFn = newFn;
        	}
        	else{
        		CLObj.ampErrLikelihoodObj.setFn(currFn, dataFlag);
        	}
        	
        	/*
        	 * Update alpha, i.e., value of Fp
        	 */
        	double newFp = SampleErrorRates.sampleNewFpDoublet(listSingleCells, listClone, c2CellCloneList, cellDoubletFlagArr, CLObj.ampErrLikelihoodObj, fpPriorDist, currFn, nMut, SCF, dataFlag);
        	nuLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listSingleCells, listClone, c2CellCloneList, cellDoubletFlagArr, newFp, currFn, nMut, SCF, dataFlag);
        	// better likelihood
        	// Update all variables
        	if (nuLikelihood > currLikelihood){
        		currLikelihood = nuLikelihood;
        		CLObj.ampErrLikelihoodObj.setFp(newFp, dataFlag);
            	currFp = newFp;
        	}
        	else{
        		CLObj.ampErrLikelihoodObj.setFp(currFp, dataFlag);
        	}
        	
        	/*
        	 * Update delta, i.e, value of _doublet
        	 */
        	double newDelta = SampleErrorRates.sampleNewDeltaDoublet(cellDoubletFlagArr, doubletPriorDist, nCell);

        	/*
        	 * Update the Doublet Indicator Vector 
        	 */
        	// Create a Copy of the vector first
        	int[] copyCellDoubletFlagArr = SCF.getCopyCellDoubletFlagArr(cellDoubletFlagArr);
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = listSingleCells.get(j);
        		Clone C_j = SCF.getCellClone(listClone, s_j);
        		Clone doubletMateClone_j = SCF.getCloneFrmList(listClone, c2CellCloneList.get(j));
        		int s_doubletFlag = SampleCellDoubletIndicator.sampleCellDoubletIndicator(s_j, C_j, doubletMateClone_j, currFp, currFn, newDelta, CLObj.ampErrLikelihoodObj, dataFlag, nMut);
        		copyCellDoubletFlagArr[j] = s_doubletFlag;
        	}
        	
        	nuLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listSingleCells, listClone, c2CellCloneList, copyCellDoubletFlagArr, currFp, currFn, nMut, SCF, dataFlag);

        	if (nuLikelihood > currLikelihood){
        		currLikelihood = nuLikelihood;
        		cellDoubletFlagArr = copyCellDoubletFlagArr;
        		currDoublet = newDelta;
        	}
        	
        	if (currLikelihood > bestL){
        		bestL = currLikelihood;
        		bestClone = listClone;
        		best2ndClone = c2CellCloneList;
        		bestTree = currTree;
        		bestFn = currFn;
        		bestFp = currFp;
        		bestDoubletFlagArr = cellDoubletFlagArr;
        		bestDoublet = currDoublet;
        	}
        	
        	// Print occasionally
        	if (i % iterP == 0){
        		System.out.println("likelihood = "+ bestL);
        	}
        	
        	alpha_0 = SCF.sampleConcentrationParam(alpha_0, listClone.size(), nCell, GammaA, GammaB);
        	
        }
        
        DoubletPosteriorSample bestMLSample = new DoubletPosteriorSample(bestClone, best2ndClone, bestTree, bestDoubletFlagArr, bestL, bestFp, bestFn, bestDoublet);
        System.out.println("best likelihood = "+ bestL);
		return bestMLSample;
	}
	
	/**
	 * 
	 * @param listClone
	 * @param nCell
	 * @return
	 * Created On: Dec 19, 2017
	 */
	private static String getCellCloneString(ArrayList<Clone> listClone, int nCell){
		int[] cellCloneID = new int[nCell];
		ArrayList<String> cellCloneIDStrList = new ArrayList<>();
		for (Clone C : listClone){
			for (int cell: C.memberCellList){
				cellCloneID[cell] = C.cloneID;
			}
		}
		for (int i = 0; i < nCell; i++)
			cellCloneIDStrList.add(Integer.toString(cellCloneID[i]));
		String cellCloneIDStr = String.join(" ", cellCloneIDStrList);
		return cellCloneIDStr;
	}
	
	/**
	 * Return string representation of cell 2nd clones
	 * @param listClone
	 * @return
	 * Created On: Dec 19, 2017
	 */
	private static String getSecondCellCloneString(ArrayList<Integer> listClone){
		ArrayList<String> cellCloneIDStrList = new ArrayList<>();
		for (Integer cellClone: listClone){
			cellCloneIDStrList.add(Integer.toString(cellClone));
		}
		String cellCloneIDStr = String.join(" ", cellCloneIDStrList);
		return cellCloneIDStr;
	}
	
	/**
	 * Return string representation of cell doublet flags
	 * @param cellDoubletFlagArr
	 * @return
	 * Created On: Dec 19, 2017
	 */
	private static String getCellDoubletFlagString(int[] cellDoubletFlagArr){
		ArrayList<String> cellDoubletFlagStringList = new ArrayList<>();
		for (int i : cellDoubletFlagArr){
			cellDoubletFlagStringList.add(Integer.toString(i));
		}
		String cellDoubletFlagString = String.join(" ", cellDoubletFlagStringList);
		return cellDoubletFlagString;
	}
	
	/**
	 * 
	 * @param listClone
	 * @return
	 * Created On: Dec 19, 2017
	 */
	private static String getCloneGenotypeListString(ArrayList<Clone> listClone){
		ArrayList<String> cloneGenotypeStringList = new ArrayList<>();
		for (Clone C : listClone){
			cloneGenotypeStringList.add(getCloneGenotypeString(C));
		}
		String cloneGenotypeString = String.join(",", cloneGenotypeStringList);
		return cloneGenotypeString;
	}
	
	/**
	 * Get the genotype string for a single clone
	 * @param C
	 * @return
	 * Created On: Oct 30, 2017
	 */
	private static String getCloneGenotypeString(Clone C) {
		ArrayList<String> cloneGtStrList = new ArrayList<>();
		for (int g : C.cloneGTVector){
			cloneGtStrList.add(Integer.toString(g));
		}
		String cloneGtStr = String.join(" ", cloneGtStrList);
		return cloneGtStr;
	}
	
	/**
	 * 
	 * @param listClone
	 * @param nCell
	 * @param singleCellNames
	 * @param currTree
	 * @return
	 * Created On: Dec 19, 2017
	 */
	private static String getCellClonalTreeString(ArrayList<Clone> listClone, int nCell, ArrayList<String> singleCellNames, STITree<Clone> currTree){
		STITree<Clone> copyCurrTree = new STITree<>(currTree);
		for (Clone C : listClone){
			ArrayList<String> thisCloneCells = new ArrayList<>();
			for (int cell: C.memberCellList){
				thisCloneCells.add(singleCellNames.get(cell));
			}
			STINode<Clone> C_node = copyCurrTree.getNode(C.cloneName);
			for (String cell: thisCloneCells){
				C_node.createChild(cell);
				copyCurrTree.getNode(cell).setParentDistance(0.000000000000001);
			}
		}
		return copyCurrTree.toNewick();
	}
	

	/**
	 * @param args
	 * Created On: Oct 30, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String s = GibbsSamplingDoublet.test();
		System.out.println(s);
	}

}
