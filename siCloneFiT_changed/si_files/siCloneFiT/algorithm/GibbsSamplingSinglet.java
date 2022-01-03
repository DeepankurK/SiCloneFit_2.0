/**
 * Oct 26, 2017
 */
package algorithm;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;

import SiFit.model.ComplexEvolutionModel;
import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import siCloneFiT.likelihood.CompleteLikelihood;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.objects.SingletPosteriorSample;
import siCloneFiT.posterior.CompleteSingletPosterior;
import siCloneFiT.proposal.TreeProposal;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Oct 26, 2017
 */
public class GibbsSamplingSinglet {
	
	public static Random _rng = new Random();
	
	/**
	 * Sample from the posterior distribution
	 * @param nCell
	 * @param nMut
	 * @param listSingleCells
	 * @param singleCellNames
	 * @param cellIDGTVectorMap
	 * @param fn
	 * @param fnPriorDist
	 * @param fp
	 * @param fpPriorDist
	 * @param _delProb
	 * @param delProbPriorDist
	 * @param _LOHProb
	 * @param LOHProbPriorDist
	 * @param _recurProb
	 * @param recurProbPriorDist
	 * @param burnInIter
	 * @param iterT
	 * @param dataFlag
	 * @param restartIndex
	 * @param SCF
	 * @return
	 * Created On: Oct 26, 2017
	 * @throws UnsupportedEncodingException 
	 * @throws FileNotFoundException 
	 */
	public static SingletPosteriorSample samplePosterior(int nCell, int nMut, ArrayList<SingleCell> listSingleCells,ArrayList<SingleCell> listSingleCells_tot,
								  						 ArrayList<String> singleCellNames, HashMap<Integer, Integer[]> cellIDGTVectorMap,HashMap<Integer, Integer[]> cellIDTotVectorMap,
								  						 double fn, BetaDistribution fnPriorDist, double fp, BetaDistribution fpPriorDist,
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
		
		listClone = SCF.removeEmptyClones(listClone); // Remove empty clone
		
		ArrayList<String> cloneNames = new ArrayList<>();
		
		// Assign Genotypes to Clones, Randomly now
		for (Clone C: listClone){
			C.assignRandomGTVector(nMut, listSingleCells, SCF);
//			System.out.println(C.cloneName);
//			System.out.println(Arrays.toString(C.cloneGTVector));
			cloneNames.add(C.cloneName);
		}
		// Obtain the genotypeObservation for clones
		ArrayList<GenotypeObservation> cloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(listClone, cloneNames, nMut);
		ArrayList<Double> currFnList = new ArrayList<Double>();       
        for(int i=1;i<=nMut;i++){
            currFnList.add(fn);
        }
		double currFp = fp;
		
        // Generate random tree
        STITree<Clone> randTreeMaker = SCF.generateRandomCloneTree(cloneNames);
        STITree<Clone> randTree = SCF.getCloneTree(randTreeMaker.toNewick());       
        STITree<Clone> currTree = randTree;
        
        // Evolution Model
        ComplexEvolutionModel m = new ComplexEvolutionModel(_delProb, _LOHProb, _recurProb);
        
        // Object for proposing new tree 
        TreeProposal proposeTreeObj = new TreeProposal(); 

        // Complete Likelihood object (error likelihood + genotype likelihood of tree)
        CompleteLikelihood CLObj = new CompleteLikelihood(currTree, m, currFp, currFnList, dataFlag);
		
        // Posterior probability of current sample 
        double posterior = CompleteSingletPosterior.computeFullModelLogPosterior(listSingleCells,listSingleCells_tot, listClone, cloneGenotypeObsMatrix, CLObj, proposeTreeObj, currTree, m, nMut, nCell, dataFlag, currFnList, fnPriorDist, currFp, fpPriorDist, _delProb, delProbPriorDist, _LOHProb, LOHProbPriorDist, alpha_0, alpha_0PriorDist);
        System.out.println("start sample posterior = "+ posterior);
        
        double likelihood = CLObj.computeFullLikelihood(listSingleCells,listSingleCells_tot, listClone, currFp, currFnList, cloneGenotypeObsMatrix, nMut);
       
         // Determine the frequency of sampling
        int sampleRate;
        if (iterT < 1000)
        	sampleRate = 10;
        else
        	sampleRate = iterT/1000;
        
        // Best configurations and likelihoods
        ArrayList<Clone> bestClone = listClone;
        double bestL = Double.NEGATIVE_INFINITY;
        double bestP = Double.NEGATIVE_INFINITY;
        STITree<Clone> bestTree = currTree;
        ArrayList<Double> bestFn = currFnList;
        double bestFp = currFp;
        double bestDelProb = _delProb;
        double bestLOHProb = _LOHProb;
        
        // Files for writing the samples
        String clonalClusterFile = outDir + "samples_" + restartIndex + "_clusters.txt";
        String clonalGenotypesFile = outDir + "samples_" + restartIndex + "_genotypes.txt";
        String clonalTreesFile = outDir + "samples_" + restartIndex + "_trees.txt";
        String posteriorErrorRatesFile = outDir + "samples_" + restartIndex + "_posteriors.txt";
        String predictedFile = outDir + "samples_" + restartIndex + "_MAP_predicted_genotype.txt";
        PrintWriter clusterWriter = new PrintWriter(clonalClusterFile, "UTF-8");
        PrintWriter genotypeWriter = new PrintWriter(clonalGenotypesFile, "UTF-8");
        PrintWriter treeWriter = new PrintWriter(clonalTreesFile, "UTF-8");
        PrintWriter posteriorWriter = new PrintWriter(posteriorErrorRatesFile, "UTF-8");
        
        int totalIter = burnInIter + iterT;
        // For loop over the number of iterations
        for (int i = 1; i < totalIter; i++){
        	
        	/******************************
        	 * Update indicators of cells *
        	 ******************************/
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = listSingleCells.get(j);
                SingleCell s_j_tot = listSingleCells_tot.get(j);
        		Clone C_j = SCF.getCellClone(listClone, s_j);
        		
        		// The clone is not singleton, new clone is to be proposed
        		if (C_j.isSingleton() == false){
        			double acRatio = SampleIndicators.proposeNewCloneNonSingletonCell(listClone, cloneNames, s_j,s_j_tot, C_j, currTree, proposeTreeObj, cloneGenotypeObsMatrix, m, nMut, nCell, dataFlag, CLObj, currFp, currFnList, SCF, alpha_0);        			
        			double rr1 = _rng.nextDouble();
        			if (rr1 <= acRatio){
//        				System.out.println("not singleton "+C_j.cloneName + " is changed, " + "cell " + j + " is in "+ SampleIndicators.nuClone.cloneName);
        				listClone = SampleIndicators.nuCloneAddedCloneList;
        				currTree = SampleIndicators.nuCloneAddedTree;
        				cloneNames = SampleIndicators.nuCloneAddedCloneNameList;
        				cloneGenotypeObsMatrix = SampleIndicators.nuCloneAddedCloneGenotypeObsMatrix;
        				SCF.updateSingleCellList(listSingleCells, SampleIndicators.nuClone, s_j);
        			}
        		}
        		// The clone is a singleton, assign the cell to another clone
        		else{
        			double acRatio = SampleIndicators.removeCloneSingletonCell(listClone, cloneNames, s_j,s_j_tot, C_j, currTree, proposeTreeObj, cloneGenotypeObsMatrix, m, nMut, nCell, dataFlag, CLObj, currFp, currFnList, SCF, alpha_0);

        			double rr1 = _rng.nextDouble();
        			if (rr1 <= acRatio){
        				listClone = SampleIndicators.cloneRemovedCloneList;
        				SCF.updateCloneList(listClone, SampleIndicators.nuCloneCellj, j);
        				currTree = SampleIndicators.cloneRemovedTree;
        				cloneNames = SampleIndicators.cloneRemovedCloneNameList;
        				cloneGenotypeObsMatrix = SampleIndicators.cloneRemovedCloneGenotypeObsMatrix;
        				SCF.updateSingleCellList(listSingleCells, listClone);
        			}
        		}
        	}
        	// Now update the indicators of non-singletons
        	// Based on Chinese restaurant process
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = listSingleCells.get(j);  
                SingleCell s_j_tot = listSingleCells_tot.get(j);      		
        		Clone C_j = SCF.getCellClone(listClone, s_j);
        		if (C_j.isSingleton() == false){
        			ArrayList<Double> unnormalizedDist = new ArrayList<>();
        			for (Clone C: listClone){
        				double val = PartialMHSteps.computeCRPErrorLikelihood(C, s_j,s_j_tot, nMut, nCell, currFp, currFnList, CLObj.ampErrLikelihoodObj);
        				unnormalizedDist.add(val);
        			}
        			Clone C_j_new = SampleCloneGenotype.sampleNewCloneCRPErrLikelihoodDist(unnormalizedDist, listClone);
        			C_j_new.memberCellList.add(s_j.cellID);
        			C_j.removeCell(s_j.cellID);
        			SCF.updateSingleCellList(listSingleCells, C_j_new, s_j);
        		}
        	}
        	
        	/**************************************
        	 * Update Tree and Model of Evolution *
        	 **************************************/
        	// Update the tree by sampling using metropolis-hastings
        	if (currTree.getLeafCount() > 3){
//        		int sampleTreeIter = 300;
        		STITree<Clone> newTree = SampleTreeModel.sampleCloneTreeMH(currTree, CLObj, _delProb, _LOHProb, sampleTreeIter, delProbPriorDist, LOHProbPriorDist, cloneGenotypeObsMatrix, nMut, dataFlag, proposeTreeObj, listClone.size());
        		currTree = newTree;        		     
        	}
        	
        	/***************************
        	 * Update Clonal Genotypes *
        	 ***************************/
        	for (int k = 0; k < listClone.size(); k++){
    			Clone C_k = listClone.get(k);
    			ArrayList<Clone> otherCloneList = SCF.getOtherCloneList(listClone, k, nMut);
    			ArrayList<String> otherCloneNames = new ArrayList<>();
    			for (Clone C: otherCloneList)
    				otherCloneNames.add(C.cloneName);
    			ArrayList<GenotypeObservation> otherClonesGenotypeObs = SCF.getCloneGenotypeObs(otherCloneList, otherCloneNames, nMut);
    			Integer[] C_k_newGTVector = SampleCloneGenotype.sampleCloneGenotypeVectorFromTreeAndErrorLikelihoodDist(otherClonesGenotypeObs, C_k, currTree, m, cellIDGTVectorMap, cellIDTotVectorMap,CLObj.ampErrLikelihoodObj, dataFlag, nMut);
    			C_k.cloneGTVector = C_k_newGTVector;
    		}
        	
        	/**********************
        	 * Update Error Rates *
        	 **********************/
            ArrayList<Double> newFn = new ArrayList<Double>(); 
        	newFn = SampleErrorRates.sampleNewFnHS(listSingleCells, listSingleCells_tot,listClone, CLObj.ampErrLikelihoodObj, currFp, fnPriorDist, nMut);
        	currFnList = newFn;

        	double newFp = SampleErrorRates.sampleNewFpHS(listSingleCells, listSingleCells_tot, listClone, CLObj.ampErrLikelihoodObj, fpPriorDist, currFnList, nMut);
        	currFp = newFp;
        	
        	cloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(listClone, cloneNames, nMut);
        	
        	posterior = CompleteSingletPosterior.computeFullModelLogPosterior(listSingleCells,listSingleCells_tot,listClone, cloneGenotypeObsMatrix, CLObj, proposeTreeObj, currTree, m, nMut, nCell, dataFlag, currFnList, fnPriorDist, newFp, fpPriorDist, CLObj.treeGenotypeLikelihoodObj.model.delProb, delProbPriorDist, CLObj.treeGenotypeLikelihoodObj.model.LOHProb, LOHProbPriorDist, alpha_0, alpha_0PriorDist);
        	likelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listSingleCells,listSingleCells_tot, listClone, currFp, currFnList, nMut);
        	
        	// Update the MAP solution
//        	if (likelihood > bestL){
        	if (posterior > bestP){
        		bestL = likelihood;
        		bestP = posterior;
        		bestClone = listClone;
        		bestTree = currTree;
        		bestFn= currFnList;
        		bestFp = currFp;
        		bestDelProb = CLObj.treeGenotypeLikelihoodObj.model.delProb;
        		bestLOHProb = CLObj.treeGenotypeLikelihoodObj.model.LOHProb;        		
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
        		// Write the clonal tree
        		String cellClonalTree = getCellClonalTreeString(listClone, nCell, singleCellNames, currTree);
        		treeWriter.printf("%s\n", cellClonalTree);
        		// Write posterior probability of sample and error rates
        		ArrayList<String> posteriorsNerrorRates = new ArrayList<>();
        		posteriorsNerrorRates.add(Double.toString(posterior));
        		posteriorsNerrorRates.add(Double.toString(likelihood));
        		posteriorsNerrorRates.add(Double.toString(currFp));
                for(int q=0;q<nMut;q++){
        		posteriorsNerrorRates.add(Double.toString(currFnList.get(q)));
                }
        		String posteriorsNerrorRatesStr = String.join(" ", posteriorsNerrorRates);
        		posteriorWriter.printf("%s\n", posteriorsNerrorRatesStr);
        		// Write the clonal genotypes
        		String cloneGenotypeListStr = getCloneGenotypeListString(listClone);
        		genotypeWriter.printf("%s\n", cloneGenotypeListStr);
        	}
        }
        
        // Close the file writers
        clusterWriter.close();
        treeWriter.close();
        posteriorWriter.close();
        genotypeWriter.close();
        
        System.out.println("MAP value = " + bestP);
        System.out.println("ML value = " + bestL);
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
		
		
//		SingletPosteriorSample MAPSample = new SingletPosteriorSample(bestClone, cellCloneID, bestTree, bestP, bestL, bestFp, bestFn); 
		SingletPosteriorSample MAPSample = new SingletPosteriorSample(bestClone, cellCloneID, bestTree, bestP, bestL, bestFp, bestFn, bestDelProb, bestLOHProb); 
		return MAPSample;
	}
	
	/**
	 * get the genotype string list for all clones
	 * @param listClone
	 * @return
	 * Created On: Oct 30, 2017
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
	 * return a string that has the info of the clone of cells
	 * meant for printing it in file
	 * @param listClone
	 * @param nCell
	 * @return
	 * Created On: Oct 26, 2017
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
	 * Return the newick string of cell clonal tree
	 * to be printed in the file
	 * @param listClone
	 * @param nCell
	 * @param singleCellNames
	 * @param currTree
	 * @return
	 * Created On: Oct 26, 2017
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
	
	public static String test(){
		String S = "hello";
		return S;
	}

	/**
	 * @param args
	 * Created On: Oct 26, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
