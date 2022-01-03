/**
 * 
 */
package siCloneFiT.main;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import com.sun.scenario.effect.impl.state.LinearConvolveKernel;

import SiFit.PerturbFn;
import SiFit.TopologyBranchPerturbations;
import SiFit.model.ComplexEvolutionModel;
import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import jeigen.DenseMatrix;
import siCloneFiT.algorithm.PartialMHSteps;
import siCloneFiT.algorithm.SampleCloneGenotype;
import siCloneFiT.algorithm.SampleCloneGenotypesMH;
import siCloneFiT.algorithm.SampleErrorRates;
import siCloneFiT.algorithm.SampleIndicators;
import siCloneFiT.algorithm.SampleTreeModel;
import siCloneFiT.algorithm.SamplingAlgos;
import siCloneFiT.io.GenotypeMatrixReader;
import siCloneFiT.likelihood.CompleteLikelihood;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.posterior.CompleteSingletPosterior;
import siCloneFiT.proposal.ModelParamProposal;
import siCloneFiT.proposal.TreeProposal;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 *
 */
public class RunSiCloneFiT {
	
	public static int nCell;							// Number of cells
	public static int nMut;								// Number of mutations
	public static double _fp = 0.0174;					// FP rate, alpha
	public static double _fnStart = 0.1375;				// FN rate, beta
	public static double _delProb = 0.05;				// Deletion Rate
	public static double _LOHProb = 0.01;				// LOH Rate
	public static double _recurProb = 0.05;				// Recurrent mutation Rate
	public static double GammaA = 1;					// Hyper-parameter for alpha_0
	public static double GammaB = 1;					// Hyper-parameter for alpha_0
	public static int iterT = 50000;					// Total number of iterations
	public static int iterP = 1000;						// Number of iterations after which likelihood will be printed
	public static int dataFlag = 1;						// Flag to indicate the data type, binary or ternary
	public static int modelFlag = 1;					// Flag to indicate which model to use
	public static String varMatFilename = null;			// Filename of the input matrix
	public static String trueTreeFilename = null;		// File containing the true tree in newick form
	public static String cellNamesFilename = null;		// File containing the cell names


	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Doing everything in main now, will modularize it later, the whole gibbs sampling is here now
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
		Random _rng = new Random();
//		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/examples/25cells.txt";
//		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/testData/100cells_noisy_genotype.txt";
//		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/testData/100_cells/dataset5/noisy_genotype_dataset5.txt";
		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/test/dataset3/noisy_genotype_dataset3.txt";
//		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/CO8_rslts/CO8.master_matrix.sifit";
//		String varFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/CRC0827_varmat.txt";
		
//		String cellNames = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/CO8_rslts/CO8_cellnames.txt";
//		String cellNames = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/cellnames_CRC0827.txt";
		// Read the genotype matrix
		GenotypeMatrixReader GMR = new GenotypeMatrixReader(varFile);
		nMut = 30;
		nCell = 100;
		HashMap<String, Integer[]> cellGTVectorMap = GMR.getCellGTVectorMap(varFile, null, nMut, nCell);
		ArrayList<String> singleCellNames = GMR.scNames;
		
		// Construct a list of single cell objects
		ArrayList<SingleCell> listSingleCells = SCF.constructSingleCellList(cellGTVectorMap, singleCellNames, nCell);
		HashMap<Integer, Integer[]> cellIDGTVectorMap = new HashMap<>(); // HashMap<cellID, observedVector>
		ArrayList<Integer> cellList = new ArrayList<>();
		for (SingleCell S : listSingleCells){
			cellIDGTVectorMap.put(S.cellID, S.observedGTVector);
			cellList.add(S.cellID);
		}
//		for (Integer i : cellIDGTVectorMap.keySet()){
//			System.out.println(i);
//			System.out.println(Arrays.toString(cellIDGTVectorMap.get(i)));
//		}
		
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
		
//		for (Clone C : listClone){
//			System.out.println(C.cloneID);
//			System.out.println(C.memberCellList);
//		}
		
		listClone = SCF.removeEmptyClones(listClone);
		
		for (Clone C : listClone){
			System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
			System.out.println(C.memberCellList);
		}
		
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
        
        // Parameterize prior distribution for _delProb
        double _delProbPriora = SCF.getBetaPriora(_delProb, _delProb*0.1);
        double _delProbPriorb = _delProbPriora * ((1/_delProb) - 1);
        BetaDistribution delProbPriorDist = new BetaDistribution(_delProbPriora, _delProbPriorb); 
        
        // Parameterize prior distribution for _LOHProb
        double _LOHProbPriora = SCF.getBetaPriora(_LOHProb, _LOHProb*0.1);
        double _LOHProbPriorb = _LOHProbPriora * ((1/_LOHProb) - 1);
        BetaDistribution LOHProbPriorDist = new BetaDistribution(_LOHProbPriora, _LOHProbPriorb);

        // Parameterize prior distribution for recurProb
        double recurProbPriora = SCF.getBetaPriora(_recurProb, _recurProb*0.3); // More variance for recurProb
        double recurProbPriorb = recurProbPriora * ((1/_recurProb) - 1);
	
        // Generate random tree
        STITree<Clone> randTreeMaker = SCF.generateRandomCloneTree(cloneNames);
        STITree<Clone> randTree = SCF.getCloneTree(randTreeMaker.toNewick());
        
        STITree<Clone> currTree = randTree;
        
        System.out.println(currTree.toNewick());
        
//        PerturbFn proposeErrorRateObj = new PerturbFn(_fnStart);
        
        // Evolution Model
        ComplexEvolutionModel m = new ComplexEvolutionModel(_delProb, _LOHProb, _recurProb);
        
        // Object for proposing new tree 
        TreeProposal proposeTreeObj = new TreeProposal(); 

        // Complete Likelihood object (error likelihood + genotype likelihood of tree)
        CompleteLikelihood CLObj = new CompleteLikelihood(randTree, m, currFp, currFn, dataFlag);
        double l = CLObj.computeFullLikelihood(listSingleCells, listClone, currFp, currFn, cloneGenotypeObsMatrix, nMut);
        System.out.println("likelihood = "+ l);
        
        // Sample posterior 
        double p = CompleteSingletPosterior.computeFullModelLogPosterior(listSingleCells, listClone, cloneGenotypeObsMatrix, CLObj, proposeTreeObj, currTree, m, nMut, nCell, dataFlag, currFn, fnPriorDist, currFp, fpPriorDist, _delProb, delProbPriorDist, _LOHProb, LOHProbPriorDist, alpha_0, alpha_0PriorDist);
        System.out.println("posterior = "+ p);
        
        int cloneIDGenerator = listClone.get(listClone.size()-1).cloneID; // This is where the clone ID starts from
//        System.out.println(cloneIDGenerator);
        
        
        ArrayList<Clone> bestClone = listClone;
        double bestL = Double.NEGATIVE_INFINITY;
        STITree<Clone> bestTree = currTree;
        
        
        int[] shuffleCells = new int[nCell];
        for (int j = 0; j < nCell; j++){
        	shuffleCells[j] = j;
        }
        
//        /**
        
        // Number of iterations for testing
        iterT = 150;
         
        
        // For loop over the number of iterations
        for (int i = 1; i < iterT; i++){
        	
//        	shuffleArray(shuffleCells);
//        	
//        	// Sample the indicators based on TreeCRP
//        	for (int j: shuffleCells){
////        	for (int j = 0; j < nCell; j++){
//        		
//        		SingleCell s_j = listSingleCells.get(j);
//        		Clone C_j = SCF.getCellClone(listClone, s_j);
//        		// Current clone of cell j is not singleton, cell can be assigned to another clone 
//        		if (C_j.isSingleton() == false){
////        			System.out.println("cell " + j + " is changing clone");
////        			System.out.println(C_j.cloneName);
//        			ArrayList<Double> indicatorDist = new ArrayList<>();
//        			
//        			// Compute the probabilities for the current clones
//        			for (Clone C: listClone){
//        				double val = PartialMHSteps.computeCRPErrorLikelihood(C, s_j, nMut, nCell, currFp, currFn, CLObj.ampErrLikelihoodObj, alpha_0);
//        				indicatorDist.add(val);
//        			}
//        			
//        			ArrayList<STITree<Clone>> nuCloneTreeList = new ArrayList<>(); // List for storing all possible trees due to addition of new clone to the current tree
//        			ArrayList<Integer[]> nuCloneGTvectorList = new ArrayList<>();  // List for stroing all possible GT vectors for the new clone depending on its position in the tree
//        			
//        			// Create a new Clone
//        			cloneIDGenerator = listClone.get(listClone.size()-1).cloneID + 1;
//        			Clone C_j_new = new Clone(cloneIDGenerator);
//        			C_j_new.setName(cloneIDGenerator);        			
//        			C_j_new.addCell(s_j.cellID);
//        			
//        			// Compute the probabilities, trees and GTvectors for instantiation fo the new clones
//        			for (Clone C: listClone){
//        				STITree<Clone> newTree = proposeTreeObj.addNode(currTree, C_j_new, C);
//        				Integer[] C_j_new_GT_vectorC = SampleCloneGenotype.sampleCloneGenotypeVectorFromTreeLikelihoodDist(cloneGenotypeObsMatrix, C_j_new, newTree, m, dataFlag, nMut);
//        				double val = PartialMHSteps.computeNewCloneTreeCRPErrorLikelihood(C_j_new_GT_vectorC, s_j, nMut, nCell, currFp, currFn, CLObj.ampErrLikelihoodObj, alpha_0, listClone.size());
//        				indicatorDist.add(val);
//        				nuCloneTreeList.add(newTree);
//        				nuCloneGTvectorList.add(C_j_new_GT_vectorC);
//        			}
//        			
////        			System.out.println(indicatorDist);
//        			ArrayList<Double> normalizedIndicatorDist = SamplingAlgos.normalizeDist(indicatorDist);
////        			System.out.println(normalizedIndicatorDist);
//        			int nuCloneIndex = SamplingAlgos.sampleInverseCDFDiscrete(normalizedIndicatorDist);
////        			System.out.println(nuCloneIndex);
//        			
//        			// cell is assigned to an existing clone
//        			if (nuCloneIndex <= listClone.size() - 1){
//        				listClone.get(nuCloneIndex).addCell(s_j.cellID);
//        				C_j.removeCell(s_j.cellID);
//            			SCF.updateSingleCellList(listSingleCells, listClone.get(nuCloneIndex), s_j);
//        			}
//        			// cell is assigned to a new clone
//        			// Update the list of clones, tree 
//        			else{
//        				C_j_new.cloneGTVector = nuCloneGTvectorList.get(nuCloneIndex - listClone.size());
//        				currTree = nuCloneTreeList.get(nuCloneIndex - listClone.size());
//        				ArrayList<Clone> nuListClone = SCF.getNuCloneList(listClone, C_j_new, C_j, s_j, nMut);
//        				ArrayList<String> nuCloneNames = SCF.getNuCloneNameList(cloneNames, C_j_new.cloneName);
//        				ArrayList<GenotypeObservation> nuCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(nuListClone, nuCloneNames, nMut);
//        				listClone = nuListClone;
//        				cloneNames = nuCloneNames;
//        				cloneGenotypeObsMatrix = nuCloneGenotypeObsMatrix;
//        				SCF.updateSingleCellList(listSingleCells, C_j_new, s_j);
//        			}
//        		}
//        		else{
//        			
//        			ArrayList<Clone> nuCloneList = new ArrayList<>();
//        			for (Clone C : listClone){
//        				if (C.cloneID != C_j.cloneID){
//        					Clone C_copy = SCF.copyClone(C, nMut);
//        					nuCloneList.add(C_copy);
//        				}
//        			}
//        			// Propose a new tree by removing the old clone
//        			currTree = proposeTreeObj.removeNode(currTree, C_j);
//        			
////        			System.out.println("new number of clones = " + nuCloneList.size());
////        			System.out.println(newTree.toNewick());
////        			if (newTree.getLeafCount() == 1)
////        				System.out.println(currTree.toNewick());
//        			
//        			// Update the id and name of clones in the new Clone List
//        			// Since, one clone is being removed, nClones = nClones - 1 (K = K-1)
//        			for (Clone C: nuCloneList){        				
//        				int C_id = nuCloneList.indexOf(C);
//        				if (C_id != C.cloneID){
//        					String cloneName = "C" + Integer.toString(C_id);      					
//        					currTree.getNode(C.cloneName).setName(cloneName);
//        					C.setId(C_id);
//        					C.setName(C_id);				
//        				}
//        			}
//        			
//        			listClone = nuCloneList;
//        			SCF.updateSingleCellList(listSingleCells, listClone);
//        			cloneNames = SCF.getNuCloneNameList(listClone);
//        			cloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(listClone, cloneNames, nMut);
//        			ArrayList<Double> indicatorDist = new ArrayList<>();
//        			
//        			// Compute the probabilities for the current clones
//        			for (Clone C: listClone){
//        				double val = PartialMHSteps.computeCRPErrorLikelihood(C, s_j, nMut, nCell, currFp, currFn, CLObj.ampErrLikelihoodObj, alpha_0);
//        				indicatorDist.add(val);
//        			}
//        			
//        			ArrayList<STITree<Clone>> nuCloneTreeList = new ArrayList<>(); // List for storing all possible trees due to addition of new clone to the current tree
//        			ArrayList<Integer[]> nuCloneGTvectorList = new ArrayList<>();  // List for storing all possible GT vectors for the new clone depending on its position in the tree
//        			
//        			// Create a new Clone
//        			cloneIDGenerator = listClone.get(listClone.size()-1).cloneID + 1;
//        			Clone C_j_new = new Clone(cloneIDGenerator);
//        			C_j_new.setName(cloneIDGenerator);        			
//        			C_j_new.addCell(s_j.cellID);
//        			
//        			// Compute the probabilities, trees and GTvectors for instantiation of the new clones
//        			for (Clone C: listClone){
//        				STITree<Clone> newTree = proposeTreeObj.addNode(currTree, C_j_new, C);
//        				Integer[] C_j_new_GT_vectorC = SampleCloneGenotype.sampleCloneGenotypeVectorFromTreeLikelihoodDist(cloneGenotypeObsMatrix, C_j_new, newTree, m, dataFlag, nMut);
//        				double val = PartialMHSteps.computeNewCloneTreeCRPErrorLikelihood(C_j_new_GT_vectorC, s_j, nMut, nCell, currFp, currFn, CLObj.ampErrLikelihoodObj, alpha_0, listClone.size());
//        				indicatorDist.add(val);
//        				nuCloneTreeList.add(newTree);
//        				nuCloneGTvectorList.add(C_j_new_GT_vectorC);
//        			}
//        			
////        			System.out.println(indicatorDist);
//        			ArrayList<Double> normalizedIndicatorDist = SamplingAlgos.normalizeDist(indicatorDist);
////        			System.out.println(normalizedIndicatorDist);
//        			int nuCloneIndex = SamplingAlgos.sampleInverseCDFDiscrete(normalizedIndicatorDist);
////        			System.out.println(nuCloneIndex);
//        			
//        			// cell is assigned to an existing clone
//        			if (nuCloneIndex <= listClone.size() - 1){
//        				listClone.get(nuCloneIndex).addCell(s_j.cellID);
////        				C_j.removeCell(s_j.cellID);
//            			SCF.updateSingleCellList(listSingleCells, listClone.get(nuCloneIndex), s_j);
//        			}
//        			// cell is assigned to a new clone
//        			// Update the list of clones, tree 
//        			else{
//        				C_j_new.cloneGTVector = nuCloneGTvectorList.get(nuCloneIndex - listClone.size());
//        				currTree = nuCloneTreeList.get(nuCloneIndex - listClone.size());
//        				ArrayList<Clone> nuListClone = SCF.getNuCloneList(listClone, C_j_new, nMut);
//        				ArrayList<String> nuCloneNames = SCF.getNuCloneNameList(cloneNames, C_j_new.cloneName);
//        				ArrayList<GenotypeObservation> nuCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(nuListClone, nuCloneNames, nMut);
//        				listClone = nuListClone;
//        				cloneNames = nuCloneNames;
//        				cloneGenotypeObsMatrix = nuCloneGenotypeObsMatrix;
//        				SCF.updateSingleCellList(listSingleCells, C_j_new, s_j);
//        			}
//        		}
////        		break;
//        	}
        	
//        	// Print to check
//			for (Clone C : listClone){
//				System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
//				System.out.println(C.memberCellList);
//			}
        	
        	
        	// Random shuffling
//        	shuffleArray(shuffleCells);
//        	System.out.println(Arrays.toString(shuffleCells));
        	
//        	/**
        	
        	// First, we are going to sample the indicators
        	// for all the cells. Loop over nCell
//        	for (int j : shuffleCells){
        	for (int j = 0; j < nCell; j++){
//        		System.out.println("cell " + j + " is changing clone, nu clone is being proposed" );
        		SingleCell s_j = listSingleCells.get(j);
        		Clone C_j = SCF.getCellClone(listClone, s_j);
        		
        		// The clone is not singleton, new clone is to be proposed
        		if (C_j.isSingleton() == false){
        			double acRatio = SampleIndicators.proposeNewCloneNonSingletonCell(listClone, cloneNames, s_j, C_j, currTree, proposeTreeObj, cloneGenotypeObsMatrix, m, nMut, nCell, dataFlag, CLObj, currFp, currFn, SCF, alpha_0);
        			
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
//        			System.out.println(listClone.size());
        			
//        			System.out.println("singleton "+C_j.cloneName + " is to be changed");
        			
        			double acRatio = SampleIndicators.removeCloneSingletonCell(listClone, cloneNames, s_j, C_j, currTree, proposeTreeObj, cloneGenotypeObsMatrix, m, nMut, nCell, dataFlag, CLObj, currFp, currFn, SCF, alpha_0);

        			double rr1 = _rng.nextDouble();
        			if (rr1 <= acRatio){
//        				System.out.println("singleton "+C_j.cloneName + " is changed, " + "cell " + j + " is in "+ SampleIndicators.nuCloneCellj.cloneName);
//        				System.out.println("cell "+j + " belongs to "+C_j_new.cloneName);
//        				C_j_new.addCell(j);
        				listClone = SampleIndicators.cloneRemovedCloneList;
        				SCF.updateCloneList(listClone, SampleIndicators.nuCloneCellj, j);
        				currTree = SampleIndicators.cloneRemovedTree;
        				cloneNames = SampleIndicators.cloneRemovedCloneNameList;
        				cloneGenotypeObsMatrix = SampleIndicators.cloneRemovedCloneGenotypeObsMatrix;
        				SCF.updateSingleCellList(listSingleCells, listClone);
//        				SCF.updateSingleCellList(listSingleCells, C_j_new, s_j);
        			}
        			
        		}
        		
        	}
        	// Print to check
//			for (Clone C : listClone){
//				System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
//				System.out.println(C.memberCellList);
//			}
//			System.out.println(currTree.toNewick());
//        	for (int j : shuffleCells){
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = listSingleCells.get(j);
        		
        		Clone C_j = SCF.getCellClone(listClone, s_j);
//        		System.out.println(s_j.cellID + " in clone " + C_j.cloneName);
        		if (C_j.isSingleton() == false){
//        			System.out.println("not singleton "+C_j.cloneName + " is to be changed");
        			ArrayList<Double> unnormalizedDist = new ArrayList<>();
        			for (Clone C: listClone){
        				double val = PartialMHSteps.computeCRPErrorLikelihood(C, s_j, nMut, nCell, currFp, currFn, CLObj.ampErrLikelihoodObj);
        				unnormalizedDist.add(val);
        			}
//        			System.out.println(unnormalizedDist);
        			Clone C_j_new = SampleCloneGenotype.sampleNewCloneCRPErrLikelihoodDist(unnormalizedDist, listClone);
        			C_j_new.memberCellList.add(s_j.cellID);
        			C_j.removeCell(s_j.cellID);
        			SCF.updateSingleCellList(listSingleCells, C_j_new, s_j);
//        			break;
        		}
        	}
        	
        	
//        	**/
        	
        	
        	
        	
        	
        	
//        	for (Clone C : listClone){
//				System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
//				System.out.println(C.memberCellList);
//			}
        	
        	// Update the tree by sampling using metropolis-hastings
        	if (currTree.getLeafCount() > 3){
        		int sampleTreeIter = 100;
        		STITree<Clone> newTree = SampleTreeModel.sampleCloneTreeMH(currTree, CLObj, _delProb, _LOHProb, sampleTreeIter, delProbPriorDist, LOHProbPriorDist, cloneGenotypeObsMatrix, nMut, dataFlag, proposeTreeObj, listClone.size());
        		currTree = newTree;        		     
        	}
        	
        	// Update the genotype of the clones
        	
//        	ArrayList<Integer[]> clonesNewGenotypeMatrix = SampleCloneGenotypesMH.sampleAllCloneGenotypesPosteriorMetropolis(listClone, cloneNames, listClone.size(), currTree, m, cellIDGTVectorMap, CLObj.ampErrLikelihoodObj, dataFlag, nMut, 500, SCF);
//        	for (int k = 0; k < listClone.size(); k++){
//        		listClone.get(k).cloneGTVector = clonesNewGenotypeMatrix.get(k);
//        	}
        	
//        	if (i%20 == 0)
//        		System.out.println("iter " + i + CLObj.ampErrLikelihoodObj.getErrorMat());
        	
//        	if (listClone.size() > 1){
//        	System.out.println("nClones = " + listClone.size());
        		for (int k = 0; k < listClone.size(); k++){
        			Clone C_k = listClone.get(k);
//        			System.out.println("now "+C_k.cloneName);
        			ArrayList<Clone> otherCloneList = SCF.getOtherCloneList(listClone, k, nMut);
        			ArrayList<String> otherCloneNames = new ArrayList<>();
        			for (Clone C: otherCloneList)
        				otherCloneNames.add(C.cloneName);
        			ArrayList<GenotypeObservation> otherClonesGenotypeObs = SCF.getCloneGenotypeObs(otherCloneList, otherCloneNames, nMut);
        			Integer[] C_k_newGTVector = SampleCloneGenotype.sampleCloneGenotypeVectorFromTreeAndErrorLikelihoodDist(otherClonesGenotypeObs, C_k, currTree, m, cellIDGTVectorMap, CLObj.ampErrLikelihoodObj, dataFlag, nMut);
        			C_k.cloneGTVector = C_k_newGTVector;
        			//        		break;
        		}
//        	}
        	
//        	System.out.println("iteration "+ i);
        	// Update the error values
//        	double newFn = SampleErrorRates.sampleNewFnPosteriorMH(listSingleCells, listClone, CLObj.ampErrLikelihoodObj, currFp, currFn, fnPriorDist, nMut, 1000);
        	double newFn = SampleErrorRates.sampleNewFnHS(listSingleCells, listClone, CLObj.ampErrLikelihoodObj, currFp, fnPriorDist, nMut);
//        	System.out.println("newFn = "+newFn);
        	CLObj.ampErrLikelihoodObj.setFn(newFn, dataFlag);
        	currFn = newFn;
        	
        	double newFp = SampleErrorRates.sampleNewFpHS(listSingleCells, listClone, CLObj.ampErrLikelihoodObj, fpPriorDist, currFn, nMut);
//        	double newFp = SampleErrorRates.sampleNewFpPosteriorMH(listSingleCells, listClone, CLObj.ampErrLikelihoodObj, currFp, currFn, fpPriorDist, nMut, 1000); 
//        	System.out.println("newFp = "+newFp);
        	CLObj.ampErrLikelihoodObj.setFp(newFp, dataFlag);
        	currFp = newFp;
        	
        	cloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(listClone, cloneNames, nMut);
        	
        	double posterior = CompleteSingletPosterior.computeFullModelLogPosterior(listSingleCells, listClone, cloneGenotypeObsMatrix, CLObj, proposeTreeObj, currTree, m, nMut, nCell, dataFlag, newFn, fnPriorDist, newFp, fpPriorDist, CLObj.treeGenotypeLikelihoodObj.model.delProb, delProbPriorDist, CLObj.treeGenotypeLikelihoodObj.model.LOHProb, LOHProbPriorDist, alpha_0, alpha_0PriorDist);
        	
        	double likelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listSingleCells, listClone, currFp, currFn, nMut);
//        	double posterior = CLObj.computeFullLikelihood(listSingleCells, listClone, currTree, m, _fp, _fnStart, cloneGenotypeObsMatrix, nMut);
//        	double partition = CLObj.computePartitionLogProb(listClone, alpha_0, nCell);
//        	double treePrior = proposeTreeObj.getTreeLogPriorProb(currTree);
//        	System.out.println("posterior =" +posterior + " " + "partition = " + partition+ " " + "treeP=" + treePrior);
//        	System.out.println("likelihood =" +likelihood);
        	System.out.println(posterior + " " + likelihood);
        	if (likelihood > bestL){
        		bestL = likelihood;
        		bestClone = listClone;
        		bestTree = currTree;
        	}
        	
//        	ArrayList<Integer> cellCloneList = new ArrayList<>();
//        	for (SingleCell s : listSingleCells)
//        		cellCloneList.add(s.cloneID);
//        	System.out.println(cellCloneList);
        	
        	alpha_0 = SCF.sampleConcentrationParam(alpha_0, listClone.size(), nCell, GammaA, GammaB);

//        	System.out.println(cloneNames);
        }
        
        System.out.println("best likelihood = " + bestL);
        System.out.println(bestTree.toNewick());
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
        
        String outFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/test/dataset3/predicted_genotype_dataset3.txt";
        
        Integer[][] predictedGenotypeMat = new Integer[nMut][nCell+1];
        SCF.writeGenotypeMatrix2File(outFile, cellList, cellPredictedGt, predictedGenotypeMat, nMut);
        
//        System.out.println(currTree.toNewick());
//        for (Clone C : listClone){
//			System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
//			System.out.println(C.memberCellList);
//			System.out.println(Arrays.toString(C.cloneGTVector));
//			ArrayList<String> thisCloneCells = new ArrayList<>();
//			
//			for (int cell: C.memberCellList){
//				thisCloneCells.add(singleCellNames.get(cell));
//			}
//			System.out.println(thisCloneCells);
//			STINode<Clone> C_node = currTree.getNode(C.cloneName);
//			for (String cell: thisCloneCells){
//				C_node.createChild(cell);
//				currTree.getNode(cell).setParentDistance(0.000000000000001);
//			}
//		}
//        System.out.println(currTree.toNewick());
        
//      **/  
        
        
        
        
        
        
        /**
        
        
//        System.out.println(cloneIDGenerator);
//        TopologyBranchPerturbations TBP = new TopologyBranchPerturbations();
        
        // For loop over the number of iterations
        for (int i = 1; i < 50; i++){
        	
        	for (int j = 0; j < nCell; j++){
//        		System.out.println("cell " + j + " is changing clone, nu clone is being proposed" );
//        		System.out.println(currTree.toNewick());
        		// Get the clone where cell j belongs to
        		SingleCell s_j = listSingleCells.get(j);
//        		System.out.println(Arrays.toString(s_j.observedGTVector));
        		Clone C_j = SCF.getCellClone(listClone, s_j);
//        		System.out.println(Arrays.toString(C_j.cloneGTVector));
        		
        		// The clone is not singleton, new clone is to be proposed
        		if (C_j.isSingleton() == false){
//        			System.out.println("not singleton "+C_j.cloneName + " is to be changed");
        			// Create a new Clone
        			cloneIDGenerator += 1;
        			Clone C_j_new = new Clone(cloneIDGenerator);
        			C_j_new.setName(cloneIDGenerator);
        			
        			C_j_new.addCell(j);
//        			s_j.cloneID = C_j_new.cloneID;
        			// Place the new clone in the tree and modify the tree
        			STITree<Clone> newTree = proposeTreeObj.addNode(currTree, C_j_new);
//        			System.out.println(newTree.toNewick());
        			
        			// We need to sample the genotype of the clone
        			Integer[] C_j_new_GT_vector = SampleCloneGenotype.sampleCloneGenotypeVectorFromTreeLikelihoodDist(cloneGenotypeObsMatrix, C_j_new, newTree, m, 1, nMut);
        			C_j_new.cloneGTVector = C_j_new_GT_vector;
//        			System.out.println(Arrays.toString(C_j_new_GT_vector));
        			
        			// Calculate each term of the acceptance ratio
        			double errorLikelihoodRatio = CLObj.ampErrLikelihoodObj.computeNewCloneErrorLikelihoodRatio(C_j_new.cloneGTVector, C_j.cloneGTVector, s_j.observedGTVector, _fp, _fnStart, nMut);
//        			System.out.println("errorLikelihoodRatio ="+ errorLikelihoodRatio);
        			ArrayList<Clone> nuListClone = SCF.getNuCloneList(listClone, C_j_new, C_j, s_j, nMut);
        			ArrayList<String> nuCloneNames = SCF.getNuCloneNameList(cloneNames, C_j_new.cloneName);
        			ArrayList<GenotypeObservation> nuCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(nuListClone, nuCloneNames, nMut);
        			double treeLikelihoodRatio = CLObj.treeGenotypeLikelihoodObj.computeTreeCloneLikelihoodRatio(nuCloneGenotypeObsMatrix, cloneGenotypeObsMatrix, newTree, currTree, nMut, dataFlag);
//        			System.out.println("treeLikelihoodRatio ="+ treeLikelihoodRatio);
        			double treePriorRatio = proposeTreeObj.getTreePriorRatio(nuCloneNames.size(), 1) * Math.exp(proposeTreeObj.getClusteringLogPriorRatio(nuListClone, listClone, alpha_0));
//        			System.out.println("treePriorRatio =" + treePriorRatio);
        			double hastingsRatio = proposeTreeObj.getHastingsRatioAddClone(alpha_0, nuCloneNames.size(), nuCloneNames.size() -1);
//        			System.out.println("hastingsRatio =" + hastingsRatio);
//            		System.out.println("jacobian = "+ proposeTreeObj.getJacobian());
//        			double acRatio = PartialMHSteps.getAcceptanceRatio(alpha_0, nCell, errorLikelihoodRatio, treeLikelihoodRatio, treePriorRatio, hastingsRatio, proposeTreeObj.getJacobian(), 1);
        			double acRatio = Math.min(1, (alpha_0/(nCell -1))*errorLikelihoodRatio);
//        			System.out.println(acRatio);
        			
        			double rr1 = _rng.nextDouble();
        			if (rr1 <= acRatio){
//        				System.out.println("not singleton "+C_j.cloneName + " is changed");
        				listClone = nuListClone;
        				currTree = newTree;
        				cloneNames = nuCloneNames;
        				cloneGenotypeObsMatrix = nuCloneGenotypeObsMatrix;
        				SCF.updateSingleCellList(listSingleCells, C_j_new, s_j);
        			}
        			// Print to check
//        			for (Clone C : listClone){
//        				System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
//        				System.out.println(C.memberCellList);
//        			}
//        			System.out.println(cloneNames);
//        			System.out.println(currTree.toNewick());
//        			for (SingleCell s: listSingleCells){
//        				System.out.printf("cell %s, ID = %d has clone = %s, %d\n", s.name, s.cellID, s.cloneName, s.cloneID);
//        			}
//        			
//        			break;
        			
        		}
        		// The clone is a singleton, assign the cell to another clone
        		else{
        			
//        		if (C_j.isSingleton() == true){
        			// Print to check
//        			System.out.println("singleton "+C_j.cloneName + " is to be changed");
//        			for (Clone C : listClone){
//        				System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
//        				System.out.println(C.memberCellList);
//        			}
        			// Compute the CRP distribution
        			ArrayList<Double> CRPdist = new ArrayList<>();
        			ArrayList<Clone> nuCloneList = new ArrayList<>();
        			for (Clone C : listClone){
        				if (C.cloneID != C_j.cloneID){
        					double C_CRPprob = (double)C.memberCellList.size()/(nCell - 1);
        					CRPdist.add(C_CRPprob);
        					Clone C_copy = SCF.copyClone(C, nMut);
        					nuCloneList.add(C_copy);
        				}
        			}
//        			System.out.println(CRPdist);
        			// Sample from the CRP distribution
        			Clone C_j_new = SamplingAlgos.sampleInverseCDFDiscrete(CRPdist, nuCloneList);
        			
        			// Update the list of single cells of this clone
//        			C_j_new.addCell(j);
//        			s_j.cloneID = C_j_new.cloneID;
        			
        			// Propose a new tree by removing the old clone
        			STITree<Clone> newTree = proposeTreeObj.removeNode(currTree, C_j);
//        			System.out.println(newTree.toNewick());
//        			System.out.println(C_j_new.cloneID + " " + C_j_new.cloneName);
        			// Calculate each term of the acceptance ratio
//        			System.out.println(Arrays.toString(C_j.cloneGTVector));
//        			System.out.println(Arrays.toString(C_j_new.cloneGTVector));
//        			System.out.println(Arrays.toString(s_j.observedGTVector));
        			
        			double errorLikelihoodRatio = CLObj.ampErrLikelihoodObj.computeNewCloneErrorLikelihoodRatio(C_j_new.cloneGTVector, C_j.cloneGTVector, s_j.observedGTVector, _fp, _fnStart, nMut);
//        			System.out.println(errorLikelihoodRatio);
        			ArrayList<String> nuCloneNames = SCF.getNuCloneNameList(nuCloneList);
        			ArrayList<GenotypeObservation> nuCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(nuCloneList, nuCloneNames, nMut);
        			double treeLikelihoodRatio = CLObj.treeGenotypeLikelihoodObj.computeTreeCloneLikelihoodRatio(nuCloneGenotypeObsMatrix, cloneGenotypeObsMatrix, newTree, currTree, nMut, dataFlag);
//        			System.out.println(treeLikelihoodRatio);
        			double treePriorRatio = proposeTreeObj.getTreePriorRatio(nuCloneNames.size(), 0) * Math.exp(proposeTreeObj.getClusteringLogPriorRatio(nuCloneList, listClone, alpha_0));
//        			System.out.println(treePriorRatio);
        			double hastingsRatio = proposeTreeObj.getHastingsRatioRemoveClone(nuCloneNames.size(), nuCloneNames.size()+1);
//        			System.out.println("hastings ratio = "+ hastingsRatio);
        			double acRatio = PartialMHSteps.getAcceptanceRatio(alpha_0, nCell, errorLikelihoodRatio, treeLikelihoodRatio, treePriorRatio, hastingsRatio, proposeTreeObj.getJacobian(), 0);
//        			double acRatio = Math.min(1, ((nCell -1)/alpha_0)*errorLikelihoodRatio);
//        			System.out.println(acRatio);
        			double rr1 = _rng.nextDouble();
        			if (rr1 <= acRatio){
//        				System.out.println("singleton "+C_j.cloneName + " is changed");
//        				System.out.println("cell "+j + " belongs to "+C_j_new.cloneName);
        				C_j_new.addCell(j);
        				listClone = nuCloneList;
        				currTree = newTree;
        				cloneNames = nuCloneNames;
        				cloneGenotypeObsMatrix = nuCloneGenotypeObsMatrix;
        				SCF.updateSingleCellList(listSingleCells, C_j_new, s_j);
        			}
//        			for (Clone C : listClone){
//        				System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
//        				System.out.println(C.memberCellList);
//        			}

        			
//        			break;
        		}
        	}
        	
//        	System.out.println("final clones");
//        	for (Clone C : listClone){
//				System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
//				System.out.println(C.memberCellList);
//			}
        	
        	for (int j = 0; j < nCell; j++){
        		SingleCell s_j = listSingleCells.get(j);
        		
        		Clone C_j = SCF.getCellClone(listClone, s_j);
//        		System.out.println(s_j.cellID + " in clone " + C_j.cloneName);
        		if (C_j.isSingleton() == false){
//        			System.out.println("not singleton "+C_j.cloneName + " is to be changed");
        			ArrayList<Double> unnormalizedDist = new ArrayList<>();
        			for (Clone C: listClone){
        				double val = PartialMHSteps.computeCRPErrorLikelihood(C, s_j, nMut, nCell, _fp, _fnStart, CLObj.ampErrLikelihoodObj);
        				unnormalizedDist.add(val);
        			}
//        			System.out.println(unnormalizedDist);
        			Clone C_j_new = SampleCloneGenotype.sampleNewCloneCRPErrLikelihoodDist(unnormalizedDist, listClone);
        			C_j_new.memberCellList.add(s_j.cellID);
        			C_j.removeCell(s_j.cellID);
        			SCF.updateSingleCellList(listSingleCells, C_j_new, s_j);
//        			break;
        		}
        	}
//        	System.out.println("final clones after 1 iteration");
//        	for (Clone C : listClone){
//				System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
//				System.out.println(C.memberCellList);
//				System.out.println(Arrays.toString(C.cloneGTVector));
//			}
        	
        	// Update the tree by sampling using metropolis-hastings
        	if (currTree.getLeafCount() > 3){
        		int sampleTreeIter = 500;
        		STITree<Clone> newTree = SampleTreeModel.sampleCloneTreeMH(currTree, CLObj, _delProb, _LOHProb, sampleTreeIter, delProbPriorDist, LOHProbPriorDist, cloneGenotypeObsMatrix, nMut, dataFlag, proposeTreeObj, listClone.size());
        		currTree = newTree;        		     
        	}
        	
        	// Update the genotype of the clones
        	for (int k = 0; k < listClone.size(); k++){
        		Clone C_k = listClone.get(k);
//        		System.out.println("now "+C_k.cloneName);
        		ArrayList<Clone> otherCloneList = SCF.getOtherCloneList(listClone, k, nMut);
        		ArrayList<String> otherCloneNames = new ArrayList<>();
        		for (Clone C: otherCloneList)
        			otherCloneNames.add(C.cloneName);
        		ArrayList<GenotypeObservation> otherClonesGenotypeObs = SCF.getCloneGenotypeObs(otherCloneList, otherCloneNames, nMut);
        		Integer[] C_k_newGTVector = SampleCloneGenotype.sampleCloneGenotypeVectorFromTreeAndErrorLikelihoodDist(otherClonesGenotypeObs, C_k, currTree, m, cellIDGTVectorMap, CLObj.ampErrLikelihoodObj, dataFlag, nMut);
        		C_k.cloneGTVector = C_k_newGTVector;
//        		break;
        	}
        	
        	// Update the error values
        	double newFn = SampleErrorRates.sampleNewFn(listSingleCells, listClone, CLObj.ampErrLikelihoodObj, _fp, fnPriorDist, nMut);
//        	System.out.println("newFn = "+newFn);
        	CLObj.ampErrLikelihoodObj.setFn(newFn, dataFlag);
        	_fnStart = newFn;
        	
        	double newFp = SampleErrorRates.sampleNewFp(listSingleCells, listClone, CLObj.ampErrLikelihoodObj, fpPriorDist, _fnStart, nMut);
//        	System.out.println("newFp = "+newFp);
        	CLObj.ampErrLikelihoodObj.setFp(newFp, dataFlag);
        	_fp = newFp;
        	
        	
//        	if (i == 9){
//        		double[] fnGridPosterior = new double[99];
//        		double maxFnPosterior = Double.MIN_VALUE;
//        		double beta = 0.01;
//        		for (int f = 0; f<99; f++){
//        			double likelihood = Math.exp(CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listSingleCells, listClone, _fp, beta, nMut));
//        			double prior = fnPriorDist.density(beta);
//        			fnGridPosterior[f] = likelihood * prior;      
//        			if (fnGridPosterior[f] > maxFnPosterior)
//        				maxFnPosterior = fnGridPosterior[f];
//        			beta += 0.01;
//        		}
//        		UniformRealDistribution fn_fX = new UniformRealDistribution(0.01, 0.99);
//        		UniformRealDistribution fn_gY = new UniformRealDistribution(0, maxFnPosterior);
//        		double x = fn_fX.sample();
//        		double x_posterior = Math.exp(CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listSingleCells, listClone, _fp, x, nMut)) * fnPriorDist.density(x);
//        		double y = fn_gY.sample();
//        		double newFn;
//        		while(true){
//        			if (y <= x_posterior){
//        				newFn = x;
//        				break;
//        			}
//        			else{
//        				x = fn_fX.sample();
//        				x_posterior = Math.exp(CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listSingleCells, listClone, _fp, x, nMut)) * fnPriorDist.density(x);
//        				y = fn_gY.sample();
//        			}
//        			
//        		}
        		
//        		System.out.println(Arrays.toString(fnGridPosterior));
        		
//        	}
//        	double newFn = proposeErrorRateObj.proposeNewFn(_fnStart, fnPriorSD);
//        	double newFnLogLikelihoodRatio = CLObj.ampErrLikelihoodObj.computeNewFnErrorLikelihoodRatio(listSingleCells, listClone, _fp, _fnStart, newFn, nMut);
//        	if (newFnLogLikelihoodRatio > 1){
//        		_fnStart = newFn;
//        		CLObj.ampErrLikelihoodObj.setFn(newFn, dataFlag);
//        		System.out.println("Fn value updated");
//        	}
        	
//        	for (Clone C : listClone){
//				System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
//				System.out.println(C.memberCellList);
//				System.out.println(Arrays.toString(C.cloneGTVector));
//			}
        	
        	cloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(listClone, cloneNames, nMut);
        	
        	double posterior = CompleteSingletPosterior.computeFullModelLogPosterior(listSingleCells, listClone, cloneGenotypeObsMatrix, CLObj, proposeTreeObj, currTree, m, nMut, nCell, newFn, fnPriorDist, newFp, fpPriorDist, CLObj.treeGenotypeLikelihoodObj.model.delProb, delProbPriorDist, CLObj.treeGenotypeLikelihoodObj.model.LOHProb, LOHProbPriorDist, alpha_0, alpha_0PriorDist);
        	
//        	double posterior = CLObj.computeFullLikelihood(listSingleCells, listClone, currTree, m, _fp, _fnStart, cloneGenotypeObsMatrix, nMut);
//        	double partition = CLObj.computePartitionLogProb(listClone, alpha_0, nCell);
//        	double treePrior = proposeTreeObj.getTreeLogPriorProb(currTree);
//        	System.out.println("posterior =" +posterior + " " + "partition = " + partition+ " " + "treeP=" + treePrior);
        	double likelihood = CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listSingleCells, listClone, _fp, _fnStart, nMut);
        	System.out.println("likelihood =" +likelihood);
        	
        	alpha_0 = SCF.sampleConcentrationParam(alpha_0, listClone.size(), nCell, GammaA, GammaB);
        	
        	
        	
        	
        	
        }
        
        
        for (Clone C : listClone){
			System.out.printf("Clone %s, ID = %d\n", C.cloneName, C.cloneID);
			System.out.println(C.memberCellList);
			System.out.println(Arrays.toString(C.cloneGTVector));
			ArrayList<String> thisCloneCells = new ArrayList<>();
			
			for (int cell: C.memberCellList){
				thisCloneCells.add(singleCellNames.get(cell));
			}
			System.out.println(thisCloneCells);
			STINode<Clone> C_node = currTree.getNode(C.cloneName);
			for (String cell: thisCloneCells){
				C_node.createChild(cell);
				currTree.getNode(cell).setParentDistance(0.000000000000001);
			}
		}
        System.out.println(currTree.toNewick());
        
        **/
        
        

	}
	
	private static void shuffleArray(int[] array)
	{
	    int index, temp;
	    Random random = new Random();
	    for (int i = array.length - 1; i > 0; i--)
	    {
	        index = random.nextInt(i + 1);
	        temp = array[index];
	        array[index] = array[i];
	        array[i] = temp;
	    }
	}

}
