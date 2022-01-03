/**
 * Oct 8, 2017
 */
package siCloneFiT.main;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;

import com.sun.xml.internal.bind.v2.runtime.unmarshaller.XsiNilLoader.Single;

import SiFit.model.ComplexEvolutionModel;
import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import siCloneFiT.algorithm.PartialMHSteps;
import siCloneFiT.algorithm.SampleCellDoubletIndicator;
import siCloneFiT.algorithm.SampleCloneGenotype;
import siCloneFiT.algorithm.SampleErrorRates;
import siCloneFiT.algorithm.SampleIndicator2Doublet;
import siCloneFiT.algorithm.SampleIndicators;
import siCloneFiT.algorithm.SampleTreeModel;
import siCloneFiT.algorithm.SamplingAlgos;
import siCloneFiT.io.GenotypeMatrixReader;
import siCloneFiT.likelihood.CompleteLikelihood;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.DoubletGenotype;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.proposal.TreeProposal;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Oct 8, 2017
 */
public class RunSiCloneFiTdoubletHS {
	
	public static int nCell;							// Number of cells
	public static int nMut;								// Number of mutations
	public static double _fp = 0.0174;					// FP rate, alpha
	public static double _fnStart = 0.1375;				// FN rate, beta
	public static double _doublet = 0.1;						// doublet rate, delta
	public static double _delProb = 0.05;				// Deletion Rate
	public static double _LOHProb = 0.01;				// LOH Rate
	public static double _recurProb = 0.05;				// Recurrent mutation Rate
	public static double GammaA = 1;					// Hyper-parameter for alpha_0
	public static double GammaB = 1;					// Hyper-parameter for alpha_0
	public static int iterT = 50000;					// Total number of iterations
	public static int iterP = 1000;						// Number of iterations after which likelihood will be printed
	public static int dataFlag = 0;						// Flag to indicate the data type, binary or ternary
	public static int modelFlag = 1;					// Flag to indicate which model to use
	public static String varMatFilename = null;			// Filename of the input matrix
	public static String trueTreeFilename = null;		// File containing the true tree in newick form
	public static String cellNamesFilename = null;		// File containing the cell names


	/**
	 * 
	 */
	public RunSiCloneFiTdoubletHS() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param args
	 * Created On: Oct 8, 2017
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		// Basic functionalities
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
		Random _rng = new Random();
		
		// Input files
//		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/examples/25cells.txt";
		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/testData/100_cells/dataset9/noisy_genotype_dataset9.txt";
		
//		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/CO8_rslts/CO8.master_matrix.sifit";
//		String cellNames = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/CO8_rslts/CO8_cellnames.txt";

		// Read the genotype matrix
		GenotypeMatrixReader GMR = new GenotypeMatrixReader(varFile);
		nMut = 50;
		nCell = 1000;
		HashMap<String, Integer[]> cellGTVectorMap = GMR.getCellGTVectorMap(varFile, null, nMut, nCell);
		ArrayList<String> singleCellNames = GMR.scNames;

		// Construct a list of single cell objects
		ArrayList<SingleCell> listSingleCells = SCF.constructSingleCellList(cellGTVectorMap, singleCellNames, nCell);
		HashMap<Integer, Integer[]> cellIDGTVectorMap = new HashMap<>(); // HashMap<cellID, observedVector>
		for (SingleCell S : listSingleCells){
			cellIDGTVectorMap.put(S.cellID, S.observedGTVector);
		}
		
		int[] cellDoubletFlagArr = new int[nCell]; // The vector Y, Y_j denotes if cell j is a doublet or not
		for (int j = 0; j < nCell; j++){
			double rr = _rng.nextDouble();
			if (rr <= _doublet)
				cellDoubletFlagArr[j] = 1;
			else
				cellDoubletFlagArr[j] = 0;
		}
//		System.out.println(Arrays.toString(cellDoubletFlagArr));
		
		
		
		// Initiate Clone objects
		int startClone = 5; // Number of clones to start with
		ArrayList<Clone> listClone = new ArrayList<>();
		for (int i = 0; i < startClone; i++){
			Clone c = new Clone(i);
			c.setNameID(i);
			listClone.add(c);
		}
		
		// Assign cells to Clones, random assignments
		for (SingleCell S : listSingleCells){
			int cellCloneID = _rng.nextInt((4 - 0) + 1) + 0;
			S.cloneID = cellCloneID;
			listClone.get(cellCloneID).memberCellList.add(S.cellID);
		}
		
		// Remove empty clones
		listClone = SCF.removeEmptyClones(listClone);
		
		for (Clone C : listClone){
			System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
			System.out.println(C.memberCellList);
		}
		
		ArrayList<String> cloneNames = new ArrayList<>();
		
		// Assign Genotypes to Clones, Randomly now
		for (Clone C: listClone){
			C.assignRandomGTVector(nMut, listSingleCells, SCF);
			cloneNames.add(C.cloneName);
		}
		
		// Assign 2nd clone to cells randomly
		ArrayList<Integer> c2CellCloneList = new ArrayList<>();
		for (int j = 0; j < nCell; j++){
			SingleCell s_j = listSingleCells.get(j);
			int s_j_clone2Id = SCF.assignRandomClone2(s_j, listClone); 
			c2CellCloneList.add(s_j_clone2Id);			
		}
//		System.out.println(c2CellCloneList);
		// Obtain the genotypeObservation for clones
		ArrayList<GenotypeObservation> cloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(listClone, cloneNames, nMut);

		// Concentration parameter
		double alpha_0 = 1; 
		GammaDistribution alpha_0PriorDist = new GammaDistribution(GammaA, GammaB);
		
		// Work on the parameters
		// Distribution for Fn
        double fnPriorMean = _fnStart;
        double fnPriorSD = fnPriorMean*0.5;
        double betaPriora = ((1 - fnPriorMean)*fnPriorMean*fnPriorMean/(fnPriorSD*fnPriorSD)) - fnPriorMean;
        double betaPriorb = betaPriora * ((1/fnPriorMean) - 1);
        BetaDistribution fnPriorDist = new BetaDistribution(betaPriora, betaPriorb);
        double currFn = _fnStart;
        
        // Distribution for Fp
        double fpPriorMean = _fp;
        double fpPriorSD = fpPriorMean * 0.5;
        double fpBetaPriora = ((1 - fpPriorMean)*fpPriorMean*fpPriorMean/(fpPriorSD*fpPriorSD)) - fpPriorMean;
        double fpBetaPriorb = fpBetaPriora * ((1/fpPriorMean) - 1);
        BetaDistribution fpPriorDist = new BetaDistribution(fpBetaPriora, fpBetaPriorb);
        double currFp = _fp;
        
        // Distribution for delta (doublet rate)
        double _doubletPriora = SCF.getBetaPriora(_doublet, _doublet*0.5);
        double _doubletPriorb = _doubletPriora * ((1/_doublet) - 1);
        BetaDistribution doubletPriorDist = new BetaDistribution(_doubletPriora, _doubletPriorb);
        double currDoublet = _doublet;
         
        // Parameterize prior distribution for _delProb
        double _delProbPriora = SCF.getBetaPriora(_delProb, _delProb*0.1);
        double _delProbPriorb = _delProbPriora * ((1/_delProb) - 1);
        BetaDistribution delProbPriorDist = new BetaDistribution(_delProbPriora, _delProbPriorb); 
        
        // Parameterize prior distribution for _LOHProb
        double _LOHProbPriora = SCF.getBetaPriora(_LOHProb, _LOHProb*0.1);
        double _LOHProbPriorb = _LOHProbPriora * ((1/_LOHProb) - 1);
        BetaDistribution LOHProbPriorDist = new BetaDistribution(_LOHProbPriora, _LOHProbPriorb);
        
        // Generate random tree
        STITree<Clone> randTreeMaker = SCF.generateRandomCloneTree(cloneNames);
        STITree<Clone> randTree = SCF.getCloneTree(randTreeMaker.toNewick());
        
        STITree<Clone> currTree = randTree;
        
        System.out.println(currTree.toNewick());
        
        // Evolution Model
        ComplexEvolutionModel m = new ComplexEvolutionModel(_delProb, _LOHProb, _recurProb);
        
        // Object for proposing new tree 
        TreeProposal proposeTreeObj = new TreeProposal(); 

        // Complete Likelihood object (error likelihood + genotype likelihood of tree)
        CompleteLikelihood CLObj = new CompleteLikelihood(randTree, m, currFp, currFn, dataFlag);
        
        double currLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listSingleCells, listClone, c2CellCloneList, cellDoubletFlagArr, currFp, currFn, nMut, SCF, dataFlag);
        System.out.println("likelihood = "+ currLikelihood);
        
        int cloneIDGenerator = listClone.get(listClone.size()-1).cloneID; // This is where the clone ID starts from
        
        // Best configurations and likelihoods
        ArrayList<Clone> bestClone = listClone;
        double bestL = Double.NEGATIVE_INFINITY;
        STITree<Clone> bestTree = currTree;
        double bestFn = currFn;
        double bestFp = currFp;
        double bestDoublet = currDoublet;
        int[] bestDoubletFlagArr = new int[nCell];
        
        // Number of iterations for testing
        iterT = 100;
        
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
        		
        	double nuLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(copyListSingleCells, copyListClone, copyC2CellCloneList, cellDoubletFlagArr, currFp, currFn, nMut, SCF, dataFlag);
        	
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
        	
//        	for (Clone C : listClone){
//				System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
//				System.out.println(C.memberCellList);
//			}
//        	
//        	System.out.println(c2CellCloneList);
        	
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
//        	System.out.println(currLikelihood);
//        	System.out.println(nuLikelihood);
//        	System.out.println(copyC2CellCloneList);
        	// Better likelihood, update corresponding variable
        	if (nuLikelihood > currLikelihood){
        		c2CellCloneList = copyC2CellCloneList;
        	}
//        	System.out.println(c2CellCloneList);
        	
        	
        	/**************************************
        	 * Update Tree and Model of Evolution *
        	 **************************************/
        	// Update the tree by sampling using metropolis-hastings
        	if (currTree.getLeafCount() > 3){
        		int sampleTreeIter = 100;
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
//        	System.out.println(Arrays.toString(copyCellDoubletFlagArr));
//        	System.out.println("iteration = " + i + " " + Arrays.toString(cellDoubletFlagArr));
        	// better likelihood
        	// Update all variables
        	if (nuLikelihood > currLikelihood){
//        		System.out.println("hello");
        		currLikelihood = nuLikelihood;
        		cellDoubletFlagArr = copyCellDoubletFlagArr;
        		currDoublet = newDelta;
        	}
        	
        	if (currLikelihood > bestL){
        		bestL = currLikelihood;
        		bestClone = listClone;
        		bestTree = currTree;
        		bestFn = currFn;
        		bestFp = currFp;
        		bestDoubletFlagArr = cellDoubletFlagArr;
        		bestDoublet = currDoublet;
        	}
//        	System.out.println("iteration = " + i + " " + Arrays.toString(cellDoubletFlagArr));
        	System.out.println("iteration = " + i + " likelihood = " + currLikelihood);
        	
        	alpha_0 = SCF.sampleConcentrationParam(alpha_0, listClone.size(), nCell, GammaA, GammaB);
        	
        	
        	
        }
        
        System.out.println("best likelihood = " + bestL);
        System.out.println(bestTree.toNewick());
        for (Clone C : bestClone){
			System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
			System.out.println(C.memberCellList);
			System.out.println(Arrays.toString(C.cloneGTVector));
			ArrayList<String> thisCloneCells = new ArrayList<>();
			
			for (int cell: C.memberCellList){
				thisCloneCells.add(singleCellNames.get(cell));
			}
			System.out.println(thisCloneCells);
			STINode<Clone> C_node = bestTree.getNode(C.cloneName);
			for (String cell: thisCloneCells){
				C_node.createChild(cell);
				bestTree.getNode(cell).setParentDistance(0.000000000000001);
			}
		}
        System.out.println(bestTree.toNewick());
        System.out.println("Fn = " + bestFn);
        System.out.println("Fp = " + bestFp);
        System.out.println("delta = " + bestDoublet);
        System.out.println(Arrays.toString(bestDoubletFlagArr));


	}

}
