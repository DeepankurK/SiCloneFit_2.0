/**
 * Oct 4, 2017
 */
package siCloneFiT.main;

import java.io.IOException;
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
import siCloneFiT.algorithm.PartialMHSteps;
import siCloneFiT.algorithm.SampleCloneGenotype;
import siCloneFiT.algorithm.SampleErrorRates;
import siCloneFiT.algorithm.SampleIndicators;
import siCloneFiT.algorithm.SampleTreeModel;
import siCloneFiT.io.GenotypeMatrixReader;
import siCloneFiT.likelihood.CompleteLikelihood;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.proposal.TreeProposal;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Oct 4, 2017
 */
public class RunSiCloneFiTHS {
	public static int nCell;							// Number of cells
	public static int nMut;								// Number of mutations
	public static double _fp = 0.3;					// FP rate, alpha
	public static double _fnStart = 0.3;				// FN rate, beta
	public static double _delProb = 0.05;				// Deletion Rate
	public static double _LOHProb = 0.01;				// LOH Rate
	public static double _recurProb = 0.05;				// Recurrent mutation Rate
	public static double GammaA = 5;					// Hyper-parameter for alpha_0
	public static double GammaB = 5;					// Hyper-parameter for alpha_0
	public static int iterT = 50000;					// Total number of iterations
	public static int iterP = 1000;						// Number of iterations after which likelihood will be printed
	public static int dataFlag = 1;						// Flag to indicate the data type, binary or ternary
	public static int modelFlag = 1;					// Flag to indicate which model to use
	public static String varMatFilename = null;			// Filename of the input matrix
	public static String trueTreeFilename = null;		// File containing the true tree in newick form
	public static String cellNamesFilename = null;		// File containing the cell names

	/**
	 * 
	 */
	public RunSiCloneFiTHS() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param args
	 * Created On: Oct 4, 2017
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		// Basic functionalities
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
		Random _rng = new Random();
		
		// Input files
//		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/examples/25cells.txt";
//		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/testData/100_cells/dataset4/noisy_genotype_dataset4.txt";
//		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/test/dataset4/noisy_genotype_dataset4.txt";
		
//		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/CRC0827_varmat.txt";
//		String cellNames = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/cellnames_CRC0827.txt";
		
		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/CO8_rslts/CO8.master_matrix.sifit";
		String cellNames = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/CO8_rslts/CO8_cellnames.txt";
		
		// Read the genotype matrix
		GenotypeMatrixReader GMR = new GenotypeMatrixReader(varFile);
		nMut = 36;
		nCell = 182;
		HashMap<String, Integer[]> cellGTVectorMap = GMR.getCellGTVectorMap(varFile, cellNames, nMut, nCell);
		ArrayList<String> singleCellNames = GMR.scNames;

		// Construct a list of single cell objects
		ArrayList<SingleCell> listSingleCells = SCF.constructSingleCellList(cellGTVectorMap, singleCellNames, nCell);
		HashMap<Integer, Integer[]> cellIDGTVectorMap = new HashMap<>(); // HashMap<cellID, observedVector>
		for (SingleCell S : listSingleCells){
			cellIDGTVectorMap.put(S.cellID, S.observedGTVector);
		}
		
		// Initiate Clone objects
		int startClone = 30; // Number of clones to start with
		ArrayList<Clone> listClone = new ArrayList<>();
		for (int i = 0; i < startClone; i++){
			Clone c = new Clone(i);
			c.setNameID(i);
			listClone.add(c);
		}
		
		// Assign cells to Clones, random assignments
		for (SingleCell S : listSingleCells){
			int cellCloneID = _rng.nextInt((29 - 0) + 1) + 0;
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
		// Obtain the genotypeObservation for clones
		ArrayList<GenotypeObservation> cloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(listClone, cloneNames, nMut);
		
		// Concentration parameter
		double alpha_0 = 10; 
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
        double l = CLObj.computeFullLikelihood(listSingleCells, listClone, currFp, currFn, cloneGenotypeObsMatrix, nMut);
        
        
        int cloneIDGenerator = listClone.get(listClone.size()-1).cloneID; // This is where the clone ID starts from
        
        // Best configurations and likelihoods
        ArrayList<Clone> bestClone = listClone;
        double bestL = Double.NEGATIVE_INFINITY;
        STITree<Clone> bestTree = currTree;
        double bestFn = currFn;
        double bestFp = currFp;
        
        // Only the error likelihood
        double currLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listSingleCells, listClone, currFp, currFn, nMut);
        // Error + Tree likelihood
//        double currLikelihood = CLObj.computeFullLikelihood(listSingleCells, listClone, currTree, m, currFp, currFn, cloneGenotypeObsMatrix, nMut, dataFlag);
        
        System.out.println("current likelihood = "+ currLikelihood);
        
        // Number of iterations for testing
        iterT = 300;
        
        // For loop over the number of iterations
        for (int i = 1; i < iterT; i++){
        	
        	// Create copies of all important variables
        	ArrayList<Clone> copyListClone = SCF.copyCloneList(listClone, nMut);
        	ArrayList<SingleCell> copyListSingleCells = SCF.copyCellList(listSingleCells, nMut);
        	STITree<Clone> copyCurrTree = new STITree<>(currTree);
        	ArrayList<String> copyCloneNames = SCF.getCloneNames(copyListClone);
        	ArrayList<GenotypeObservation> copyCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(copyListClone, copyCloneNames, nMut);
        	
        	// First, we are going to sample the indicators
        	// for all the cells. Loop over nCell
        	for (int j = 0; j < nCell; j++){
        		
        		SingleCell s_j = copyListSingleCells.get(j);
        		Clone C_j = SCF.getCellClone(copyListClone, s_j);
//        		System.out.println("C_j " + C_j.cloneName + " s_j "+ s_j.name);
        		// The clone is not singleton, new clone is to be proposed
        		if (C_j.isSingleton() == false){
        			double acRatio = SampleIndicators.proposeNewCloneNonSingletonCell(copyListClone, copyCloneNames, s_j, C_j, copyCurrTree, proposeTreeObj, copyCloneGenotypeObsMatrix, m, nMut, nCell, dataFlag, CLObj, currFp, currFn, SCF, alpha_0);
        			
        			double rr1 = _rng.nextDouble();
        			if (rr1 <= acRatio){
//        				System.out.println("not singleton "+C_j.cloneName + " is changed, " + "cell " + j + " is in "+ SampleIndicators.nuClone.cloneName);
        				copyListClone = SampleIndicators.nuCloneAddedCloneList;
        				copyCurrTree = SampleIndicators.nuCloneAddedTree;
        				copyCloneNames = SampleIndicators.nuCloneAddedCloneNameList;
        				copyCloneGenotypeObsMatrix = SampleIndicators.nuCloneAddedCloneGenotypeObsMatrix;
        				SCF.updateSingleCellList(copyListSingleCells, SampleIndicators.nuClone, s_j);
        			}
        		}
        		// The clone is a singleton, assign the cell to another clone
        		else{
//        			System.out.println("singleton "+C_j.cloneName + " is to be changed");        			
        			double acRatio = SampleIndicators.removeCloneSingletonCell(copyListClone, copyCloneNames, s_j, C_j, copyCurrTree, proposeTreeObj, copyCloneGenotypeObsMatrix, m, nMut, nCell, dataFlag, CLObj, currFp, currFn, SCF, alpha_0);

        			double rr1 = _rng.nextDouble();
        			if (rr1 <= acRatio){
//        				System.out.println("singleton "+C_j.cloneName + " is changed, " + "cell " + j + " is in "+ SampleIndicators.nuCloneCellj.cloneName);
//        				System.out.println("cell "+j + " belongs to "+C_j_new.cloneName);
        				copyListClone = SampleIndicators.cloneRemovedCloneList;
        				SCF.updateCloneList(copyListClone, SampleIndicators.nuCloneCellj, j);
        				copyCurrTree = SampleIndicators.cloneRemovedTree;
        				copyCloneNames = SampleIndicators.cloneRemovedCloneNameList;
        				copyCloneGenotypeObsMatrix = SampleIndicators.cloneRemovedCloneGenotypeObsMatrix;
        				SCF.updateSingleCellList(copyListSingleCells, copyListClone);
        			}
        		}
        	}
        	
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = copyListSingleCells.get(j);
        		
        		Clone C_j = SCF.getCellClone(copyListClone, s_j);
//        		System.out.println(s_j.cellID + " in clone " + C_j.cloneName);
        		if (C_j.isSingleton() == false){
//        			System.out.println("not singleton "+C_j.cloneName + " is to be changed");
        			ArrayList<Double> unnormalizedDist = new ArrayList<>();
        			for (Clone C: copyListClone){
        				double val = PartialMHSteps.computeCRPErrorLikelihood(C, s_j, nMut, nCell, currFp, currFn, CLObj.ampErrLikelihoodObj);
        				unnormalizedDist.add(val);
        			}
//        			System.out.println(unnormalizedDist);
        			Clone C_j_new = SampleCloneGenotype.sampleNewCloneCRPErrLikelihoodDist(unnormalizedDist, copyListClone);
        			C_j_new.memberCellList.add(s_j.cellID);
        			C_j.removeCell(s_j.cellID);
        			SCF.updateSingleCellList(copyListSingleCells, C_j_new, s_j);
        		}
        	}
        	
        	// Only the error likelihood
        	double nuLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(copyListSingleCells, copyListClone, currFp, currFn, nMut);
        	// Tree + error likelihood
//        	double nuLikelihood = CLObj.computeFullLikelihood(copyListSingleCells, copyListClone, copyCurrTree, m, currFp, currFn, copyCloneGenotypeObsMatrix, nMut, dataFlag);
        	
        	// better likelihood
        	// Update all variables
        	if (nuLikelihood > currLikelihood){
        		currLikelihood = nuLikelihood;
        		listClone = copyListClone;
        		listSingleCells = copyListSingleCells;
        		cloneNames = copyCloneNames;
        		cloneGenotypeObsMatrix = copyCloneGenotypeObsMatrix;
        		currTree = copyCurrTree;
        	}
        	
        	// Update the tree by sampling using metropolis-hastings
        	if (currTree.getLeafCount() > 3){
        		int sampleTreeIter = 100;
        		STITree<Clone> newTree = SampleTreeModel.sampleCloneTreeMH(currTree, CLObj, _delProb, _LOHProb, sampleTreeIter, delProbPriorDist, LOHProbPriorDist, cloneGenotypeObsMatrix, nMut, dataFlag, proposeTreeObj, listClone.size());
        		copyCurrTree = newTree;        		     
        	}
        	else{
        		copyCurrTree = currTree;
        	}
        	
//        	// Tree + error likelihood
//        	nuLikelihood = CLObj.computeFullLikelihood(listSingleCells, listClone, copyCurrTree, m, currFp, currFn, cloneGenotypeObsMatrix, nMut, dataFlag);
//        	// better likelihood
//        	// Update all variables
//        	if (nuLikelihood > currLikelihood){
//        		currLikelihood = nuLikelihood;
//        		currTree = copyCurrTree;
//        	}
        	
        	
        	// Create copies of all important variables
        	copyListClone = SCF.copyCloneList(listClone, nMut);
        	
    		for (int k = 0; k < copyListClone.size(); k++){
    			Clone C_k = copyListClone.get(k);
    			ArrayList<Clone> otherCloneList = SCF.getOtherCloneList(copyListClone, k, nMut);
    			ArrayList<String> otherCloneNames = new ArrayList<>();
    			for (Clone C: otherCloneList)
    				otherCloneNames.add(C.cloneName);
    			ArrayList<GenotypeObservation> otherClonesGenotypeObs = SCF.getCloneGenotypeObs(otherCloneList, otherCloneNames, nMut);
    			Integer[] C_k_newGTVector = SampleCloneGenotype.sampleCloneGenotypeVectorFromTreeAndErrorLikelihoodDist(otherClonesGenotypeObs, C_k, copyCurrTree, m, cellIDGTVectorMap, CLObj.ampErrLikelihoodObj, dataFlag, nMut);
    			C_k.cloneGTVector = C_k_newGTVector;
    			//        		break;
    		}
    		
    		copyCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(copyListClone, cloneNames, nMut);
    		
    		// Only error likelihood
    		nuLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listSingleCells, copyListClone, currFp, currFn, nMut);
    		// Tree + error likelihood
//    		nuLikelihood = CLObj.computeFullLikelihood(listSingleCells, copyListClone, currTree, m, currFp, currFn, copyCloneGenotypeObsMatrix, nMut, dataFlag);
    		
    		// better likelihood
        	// Update all variables
        	if (nuLikelihood > currLikelihood){
        		currLikelihood = nuLikelihood;
        		listClone = copyListClone;
        		cloneGenotypeObsMatrix = copyCloneGenotypeObsMatrix;
        		currTree = copyCurrTree;
        	}
        	
        	// Update the error values
        	double newFn = SampleErrorRates.sampleNewFnHS(listSingleCells, listClone, CLObj.ampErrLikelihoodObj, currFp, fnPriorDist, nMut);
        	// Only error likelihood
        	nuLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixFpFn(listSingleCells, listClone, currFp, newFn, nMut);
        	// Tree + error likelihood
//        	nuLikelihood = CLObj.computeFullLikelihood(listSingleCells, listClone, currTree, m, currFp, newFn, cloneGenotypeObsMatrix, nMut, dataFlag);
        	
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
        	
        	double newFp = SampleErrorRates.sampleNewFpHS(listSingleCells, listClone, CLObj.ampErrLikelihoodObj, fpPriorDist, currFn, nMut);
        	// Only error likelihood
        	nuLikelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixFpFn(listSingleCells, listClone, newFp, currFn, nMut);
        	// Tree + error likelihood
//        	nuLikelihood = CLObj.computeFullLikelihood(listSingleCells, listClone, currTree, m, newFp, currFn, cloneGenotypeObsMatrix, nMut, dataFlag);
        	
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
        	
        	
        	if (currLikelihood > bestL){
        		bestL = currLikelihood;
        		bestClone = listClone;
        		bestTree = currTree;
        		bestFn = currFn;
        		bestFp = currFp;
        	}
        	
//        	System.out.println("iteration = " + i + " likelihood = " + currLikelihood);
        	System.out.println(currLikelihood);
        	alpha_0 = SCF.sampleConcentrationParam(alpha_0, listClone.size(), nCell, GammaA, GammaB);
        	
        	
//        	for (Clone C : listClone){
//				System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
//				System.out.println(C.memberCellList);
//			}
        	
        	
        	
        }
        
        
        System.out.println("best likelihood = " + bestL);
        System.out.println(bestTree.toNewick());
        int[] cellCloneID = new int[nCell];
        for (Clone C : bestClone){
			System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
			System.out.println(C.memberCellList);
			System.out.println(Arrays.toString(C.cloneGTVector));
			ArrayList<String> thisCloneCells = new ArrayList<>();
			
			for (int cell: C.memberCellList){
				thisCloneCells.add(singleCellNames.get(cell));
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
        System.out.println("Fn = " + bestFn);
        System.out.println("Fp = " + bestFp);
        System.out.println(Arrays.toString(cellCloneID));
        
	}

}











