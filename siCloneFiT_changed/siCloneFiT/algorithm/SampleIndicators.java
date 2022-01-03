/**
 * Sep 21, 2017
 */
package siCloneFiT.algorithm;

import java.util.ArrayList;

import SiFit.model.ComplexEvolutionModel;
import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import siCloneFiT.likelihood.CompleteLikelihood;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.DoubletGenotype;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.proposal.TreeProposal;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Sep 21, 2017
 */
public class SampleIndicators {
	
	public static ArrayList<Clone> nuCloneAddedCloneList;
	public static ArrayList<String> nuCloneAddedCloneNameList;
	public static STITree<Clone> nuCloneAddedTree;
	public static ArrayList<GenotypeObservation> nuCloneAddedCloneGenotypeObsMatrix;
	public static Clone nuClone;
	
	public static ArrayList<Clone> cloneRemovedCloneList;
	public static ArrayList<String> cloneRemovedCloneNameList;
	public static STITree<Clone> cloneRemovedTree;
	public static ArrayList<GenotypeObservation> cloneRemovedCloneGenotypeObsMatrix;
	public static Clone nuCloneCellj;
	public static ArrayList<Integer> cloneRemovedC2CellCloneList;
	
	
	/**
	 * Create a new Clone for the current cell, 
	 * Creates updated versions of listClone, cloneGenotypeObsMatrix and cloneTree
	 * @param currListClone, current list of Clone
	 * @param cloneNames, list of clone names
	 * @param s_j, Single cell j, for which we are sampling the indicator
	 * @param s_j_tot, Single cell j in the total, from which we are sampling the indicator 
	 * @param C_j, current clone of s_j
	 * @param currTree, current tree
	 * @param proposeTreeObj, Tree proposal object
	 * @param cloneGenotypeObsMatrix, current clone genotypes
	 * @param m, evolution model
	 * @param nMut
	 * @param nCell
	 * @param df, daataFlag
	 * @param CLObj, Likelihood calculation object
	 * @param _fp
	 * @param _fnStart
	 * @param SCF
	 * @param alpha_0
	 * @return
	 * Created On: Sep 21, 2017
	 */
	public static double proposeNewCloneNonSingletonCell(ArrayList<Clone> currListClone, ArrayList<String> cloneNames,
														  SingleCell s_j, SingleCell s_j_tot, Clone C_j,
														  STITree<Clone> currTree, TreeProposal proposeTreeObj,
														  ArrayList<GenotypeObservation> cloneGenotypeObsMatrix,
														  ComplexEvolutionModel m, int nMut, int nCell, int df,
														  CompleteLikelihood CLObj, double _fp, ArrayList<Double> _fnStart,
														  SCFUtilityFunctions SCF, double alpha_0){
		// Create a new Clone
		int cloneIDGenerator = currListClone.get(currListClone.size()-1).cloneID + 1;
		Clone C_j_new = new Clone(cloneIDGenerator);
		C_j_new.setName(cloneIDGenerator);        			
		C_j_new.addCell(s_j.cellID);
		
		// Place the new clone in the tree and modify the tree
		STITree<Clone> newTree = proposeTreeObj.addNode(currTree, C_j_new);
		
		// Sample the genotype of the clone, new Metropolis try 
		Integer[] C_j_new_GT_vector = SampleCloneGenotype.sampleCloneGenotypeVectorPriorTreeMetropolis(cloneGenotypeObsMatrix, currListClone.size(), C_j_new, s_j, newTree, m, df, nMut, 200);
		
		// We need to sample the genotype of the clone
//		Integer[] C_j_new_GT_vector = SampleCloneGenotype.sampleCloneGenotypeVectorFromTreeLikelihoodDist(cloneGenotypeObsMatrix, C_j_new, newTree, m, df, nMut);
		C_j_new.cloneGTVector = C_j_new_GT_vector;
		
		// Calculate errorLikelihoodRatio
		double errorLikelihoodRatio = CLObj.ampErrLikelihoodObj.computeNewCloneErrorLikelihoodRatio(C_j_new.cloneGTVector, C_j.cloneGTVector, s_j.observedGTVector, s_j_tot.observedGTVector,_fp, _fnStart, nMut);
//		System.out.println("errorLikelihoodRatio ="+ errorLikelihoodRatio);
		
		// Construct new ListClone, new ListCloneNames and new CloneGenotypeObsMatrix
		ArrayList<Clone> nuListClone = SCF.getNuCloneList(currListClone, C_j_new, C_j, s_j, nMut);
		ArrayList<String> nuCloneNames = SCF.getNuCloneNameList(cloneNames, C_j_new.cloneName);
		ArrayList<GenotypeObservation> nuCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(nuListClone, nuCloneNames, nMut);
		
		// The updated objects are stored as static entries
		nuCloneAddedCloneList = nuListClone;
		nuCloneAddedCloneNameList = nuCloneNames;
		nuCloneAddedCloneGenotypeObsMatrix = nuCloneGenotypeObsMatrix;
		nuCloneAddedTree = newTree;
		nuClone = C_j_new;

		// These are calculation of other terms 
		double treeLikelihoodRatio = CLObj.treeGenotypeLikelihoodObj.computeTreeCloneLikelihoodRatio(nuCloneGenotypeObsMatrix, cloneGenotypeObsMatrix, newTree, currTree, nMut, df);
//		System.out.println("treeLikelihoodRatio ="+ treeLikelihoodRatio);
		double treePriorRatio = proposeTreeObj.getTreePriorRatio(nuCloneNames.size(), 1) * Math.exp(proposeTreeObj.getClusteringLogPriorRatio(nuListClone, currListClone, alpha_0));
//		double treePriorRatio = proposeTreeObj.getTreePriorRatio(nuCloneNames.size(), 1); // ignoring the clustering prior prob now
//		double treePriorRatio = proposeTreeObj.getTreePriorRatio(nuCloneNames.size(), nCell, 1);
		
		
//		double treePriorRatio = Math.exp(proposeTreeObj.getClusteringLogPriorRatio(nuListClone, currListClone, alpha_0));
		
//		double treePriorRatio = 1;
//		System.out.println("treePriorRatio =" + treePriorRatio);
		double hastingsRatio = proposeTreeObj.getHastingsRatioAddClone(alpha_0, nuCloneNames.size(), nuCloneNames.size() -1);
//		System.out.println("hastingsRatio =" + hastingsRatio);
//		System.out.println("jacobian = "+ proposeTreeObj.getJacobian());
		double acRatio = PartialMHSteps.getAcceptanceRatio(alpha_0, nCell, errorLikelihoodRatio, treeLikelihoodRatio, treePriorRatio, hastingsRatio, proposeTreeObj.getJacobian(), 1);

//		double acRatio = Math.min(1, (alpha_0/(nuListClone.size()*(nCell -1)))*errorLikelihoodRatio);
		return acRatio;
	}
	
	/**
	 * Create a new Clone for the current cell, in the presence of doublet
	 * Creates updated versions of listClone, cloneGenotypeObsMatrix and cloneTree
	 * @param currListClone
	 * @param cloneNames
	 * @param s_j
	 * @param C_j
	 * @param c2CellCloneList
	 * @param cellDoubletFlagArr
	 * @param j
	 * @param currTree
	 * @param proposeTreeObj
	 * @param cloneGenotypeObsMatrix
	 * @param m
	 * @param nMut
	 * @param nCell
	 * @param df
	 * @param CLObj
	 * @param _fp
	 * @param _fnStart
	 * @param SCF
	 * @param alpha_0
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static double proposeNewCloneNonSingletonCellDoublet(ArrayList<Clone> currListClone, ArrayList<String> cloneNames,
																SingleCell s_j, Clone C_j, ArrayList<Integer> c2CellCloneList,
																int[] cellDoubletFlagArr, int j,
																STITree<Clone> currTree, TreeProposal proposeTreeObj,
																ArrayList<GenotypeObservation> cloneGenotypeObsMatrix,
																ComplexEvolutionModel m, int nMut, int nCell, int df,
																CompleteLikelihood CLObj, double _fp, double _fnStart,
																SCFUtilityFunctions SCF, double alpha_0){
		// Create a new Clone
		int cloneIDGenerator = currListClone.get(currListClone.size()-1).cloneID + 1;
		Clone C_j_new = new Clone(cloneIDGenerator);
		C_j_new.setName(cloneIDGenerator);        			
		C_j_new.addCell(s_j.cellID);

		// Place the new clone in the tree and modify the tree
		STITree<Clone> newTree = proposeTreeObj.addNode(currTree, C_j_new);

		// We need to sample the genotype of the clone
		Integer[] C_j_new_GT_vector = SampleCloneGenotype.sampleCloneGenotypeVectorFromTreeLikelihoodDist(cloneGenotypeObsMatrix, C_j_new, newTree, m, df, nMut);
		C_j_new.cloneGTVector = C_j_new_GT_vector;
		
		// Calculate errorLikelihoodRatio
		double errorLikelihoodRatio;
		// cell is not doublet
		if (cellDoubletFlagArr[j] == 0){
			errorLikelihoodRatio = CLObj.ampErrLikelihoodObj.computeNewCloneErrorLikelihoodRatio(C_j_new.cloneGTVector, C_j.cloneGTVector, s_j.observedGTVector, _fp, _fnStart, nMut);
		}
		// cell is doublet
		else{
			Clone s_j_doubletPairClone = SCF.getCloneFrmList(currListClone, c2CellCloneList.get(j));
			Integer[] newDoubletCellGtVector = DoubletGenotype.getDoubletGtVector(C_j_new.cloneGTVector, s_j_doubletPairClone.cloneGTVector, nMut, df);
			Integer[] oldDoubletCellGtVector = DoubletGenotype.getDoubletGtVector(C_j.cloneGTVector, s_j_doubletPairClone.cloneGTVector, nMut, df);
			errorLikelihoodRatio = CLObj.ampErrLikelihoodObj.computeNewCloneErrorLikelihoodRatio(newDoubletCellGtVector, oldDoubletCellGtVector, s_j.observedGTVector, _fp, _fnStart, nMut);
		}
		
		// Construct new ListClone, new ListCloneNames and new CloneGenotypeObsMatrix
		ArrayList<Clone> nuListClone = SCF.getNuCloneList(currListClone, C_j_new, C_j, s_j, nMut);
		ArrayList<String> nuCloneNames = SCF.getNuCloneNameList(cloneNames, C_j_new.cloneName);
		ArrayList<GenotypeObservation> nuCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(nuListClone, nuCloneNames, nMut);
		
		// The updated objects are stored as static entries
		nuCloneAddedCloneList = nuListClone;
		nuCloneAddedCloneNameList = nuCloneNames;
		nuCloneAddedCloneGenotypeObsMatrix = nuCloneGenotypeObsMatrix;
		nuCloneAddedTree = newTree;
		nuClone = C_j_new;

		// These are calculation of other terms 
		double treeLikelihoodRatio = CLObj.treeGenotypeLikelihoodObj.computeTreeCloneLikelihoodRatio(nuCloneGenotypeObsMatrix, cloneGenotypeObsMatrix, newTree, currTree, nMut, df);
		double treePriorRatio = proposeTreeObj.getTreePriorRatio(nuCloneNames.size(), 1) * Math.exp(proposeTreeObj.getClusteringLogPriorRatio(nuListClone, currListClone, alpha_0));
//		System.out.println("treePriorRatio =" + treePriorRatio);
		double hastingsRatio = proposeTreeObj.getHastingsRatioAddClone(alpha_0, nuCloneNames.size(), nuCloneNames.size() -1);
//		System.out.println("hastingsRatio =" + hastingsRatio);
//		System.out.println("jacobian = "+ proposeTreeObj.getJacobian());
		double acRatio = PartialMHSteps.getAcceptanceRatio(alpha_0, nCell, errorLikelihoodRatio, treeLikelihoodRatio, treePriorRatio, hastingsRatio, proposeTreeObj.getJacobian(), 1);

//		double acRatio = Math.min(1, (alpha_0/(nuListClone.size()*(nCell -1)))*errorLikelihoodRatio);
		return acRatio;
	}
	
	/**
	 * Remove a singleton Clone from the list of clones and
	 * assign the cell s_j to a new clone
	 * @param listClone
	 * @param cloneNames
	 * @param s_j
	 * @param C_j
	 * @param currTree
	 * @param proposeTreeObj
	 * @param cloneGenotypeObsMatrix
	 * @param m
	 * @param nMut
	 * @param nCell
	 * @param df
	 * @param CLObj
	 * @param _fp
	 * @param _fnStart
	 * @param SCF
	 * @param alpha_0
	 * @return
	 * Created On: Sep 21, 2017
	 */
	public static double removeCloneSingletonCell(ArrayList<Clone> listClone, ArrayList<String> cloneNames,
												  SingleCell s_j,  SingleCell s_j_tot,Clone C_j, 
												  STITree<Clone> currTree, TreeProposal proposeTreeObj,
												  ArrayList<GenotypeObservation> cloneGenotypeObsMatrix,
												  ComplexEvolutionModel m, int nMut, int nCell, int df,
												  CompleteLikelihood CLObj, double _fp, ArrayList<Double> _fnStart,
												  SCFUtilityFunctions SCF, double alpha_0){
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
		
		// Propose a new tree by removing the old clone
		STITree<Clone> newTree = proposeTreeObj.removeNode(currTree, C_j);
		
//		System.out.println("new number of clones = " + nuCloneList.size());
//		System.out.println(newTree.toNewick());
//		if (newTree.getLeafCount() == 1)
//			System.out.println(currTree.toNewick());
		
		// Update the id and name of clones in the new Clone List
		// Since, one clone is being removed, nClones = nClones - 1 (K = K-1)
		for (Clone C: nuCloneList){        				
			int C_id = nuCloneList.indexOf(C);
			if (C_id != C.cloneID){
				String cloneName = "C" + Integer.toString(C_id);
				
				newTree.getNode(C.cloneName).setName(cloneName);
				C.setId(C_id);
				C.setName(C_id);				
			}
		}
		
		// Sample from the CRP distribution
		Clone C_j_new = SamplingAlgos.sampleInverseCDFDiscrete(CRPdist, nuCloneList);
		
		// Calculate errorLikelihoodRatio
		double errorLikelihoodRatio = CLObj.ampErrLikelihoodObj.computeNewCloneErrorLikelihoodRatio(C_j_new.cloneGTVector, C_j.cloneGTVector, s_j.observedGTVector,s_j_tot.observedGTVector, _fp, _fnStart, nMut);

		// Construct new ListCloneNames and new CloneGenotypeObsMatrix
		ArrayList<String> nuCloneNames = SCF.getNuCloneNameList(nuCloneList);
		ArrayList<GenotypeObservation> nuCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(nuCloneList, nuCloneNames, nMut);
		
		// The updated objects are stored as static entries
		cloneRemovedCloneList = nuCloneList;
		cloneRemovedCloneNameList = nuCloneNames;
		cloneRemovedCloneGenotypeObsMatrix = nuCloneGenotypeObsMatrix;
		cloneRemovedTree = newTree;
		nuCloneCellj = C_j_new;
		
		double treeLikelihoodRatio = CLObj.treeGenotypeLikelihoodObj.computeTreeCloneLikelihoodRatio(nuCloneGenotypeObsMatrix, cloneGenotypeObsMatrix, newTree, currTree, nMut, df);
//		System.out.println(treeLikelihoodRatio);
		double treePriorRatio = proposeTreeObj.getTreePriorRatio(nuCloneNames.size(), 0) * Math.exp(proposeTreeObj.getClusteringLogPriorRatio(nuCloneList, listClone, alpha_0));
//		double treePriorRatio = proposeTreeObj.getTreePriorRatio(nuCloneNames.size(), 0);
//		double treePriorRatio = proposeTreeObj.getTreePriorRatio(nuCloneNames.size(), nCell, 0);
		
//		double treePriorRatio = Math.exp(proposeTreeObj.getClusteringLogPriorRatio(nuCloneList, listClone, alpha_0));
//		double treePriorRatio = 1;
//		System.out.println(treePriorRatio);
		double hastingsRatio = proposeTreeObj.getHastingsRatioRemoveClone(nuCloneNames.size(), nuCloneNames.size()+1);
//		System.out.println("hastings ratio = "+ hastingsRatio);
		double acRatio = PartialMHSteps.getAcceptanceRatio(alpha_0, nCell, errorLikelihoodRatio, treeLikelihoodRatio, treePriorRatio, hastingsRatio, proposeTreeObj.getJacobian(), 0);
//		double acRatio = Math.min(1, ((nCell -1)/alpha_0)*errorLikelihoodRatio);
		return acRatio;
	}
	
	/**
	 * Remove a singleton Clone from the list of clones and
	 * assign the cell s_j to a new clone
	 * Make corresponding changes to c2CellCloneList and other variables
	 * @param listClone
	 * @param cloneNames
	 * @param s_j
	 * @param C_j
	 * @param Inc2CellCloneList
	 * @param cellDoubletFlagArr
	 * @param j
	 * @param currTree
	 * @param proposeTreeObj
	 * @param cloneGenotypeObsMatrix
	 * @param m
	 * @param nMut
	 * @param nCell
	 * @param df
	 * @param CLObj
	 * @param _fp
	 * @param _fnStart
	 * @param SCF
	 * @param alpha_0
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static double removeCloneSingletonCellDoublet(ArrayList<Clone> listClone, ArrayList<String> cloneNames,
														  SingleCell s_j, Clone C_j, ArrayList<Integer> Inc2CellCloneList,
														  int[] cellDoubletFlagArr, int j,
														  STITree<Clone> currTree, TreeProposal proposeTreeObj,
														  ArrayList<GenotypeObservation> cloneGenotypeObsMatrix,
														  ComplexEvolutionModel m, int nMut, int nCell, int df,
														  CompleteLikelihood CLObj, double _fp, double _fnStart,
														  SCFUtilityFunctions SCF, double alpha_0){
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
		
		// Propose a new tree by removing the old clone
		STITree<Clone> newTree = proposeTreeObj.removeNode(currTree, C_j);

		// Update the second clone List, intermediate representation
		// The cells whose 2nd clone is now C_j, are replaced by -1, will be set to C_j_new after it is sampled
		// The 2nd clone for s_j is now set to -2, will be sampled later (after C_j_new is sampled)
		// If the index of 2nd clone is higher than C_j.cloneID, then those cloneIDs are decreased by 1 as clone is removed
		
		ArrayList<Integer> c2CellCloneList = SCF.getCopyC2CellCloneList(Inc2CellCloneList);
//		System.out.println(c2CellCloneList);
//		System.out.println("s_j = " + j); 
//		System.out.println("C_j.cloneID = "+ C_j.cloneID);
		c2CellCloneList.set(j, -2);
		for (int i = 0; i < nCell; i++){
			if (c2CellCloneList.get(i) == C_j.cloneID)
				c2CellCloneList.set(i, -1);
			else if (c2CellCloneList.get(i) > C_j.cloneID)
				c2CellCloneList.set(i, c2CellCloneList.get(i)-1);
		}
//		System.out.println(c2CellCloneList);
		// Update the id and name of clones in the new Clone List
		// Since, one clone is being removed, nClones = nClones - 1 (K = K-1)
		for (Clone C: nuCloneList){        				
			int C_id = nuCloneList.indexOf(C);
			if (C_id != C.cloneID){
				String cloneName = "C" + Integer.toString(C_id);
				
				newTree.getNode(C.cloneName).setName(cloneName);
				C.setId(C_id);
				C.setName(C_id);				
			}
		}
		
		// Sample from the CRP distribution
		Clone C_j_new = SamplingAlgos.sampleInverseCDFDiscrete(CRPdist, nuCloneList);
//		System.out.println("C_j_new = "+C_j_new.cloneID);
		// Update the second clone List,
		// The cells whose 2nd clone were C_j are now replaced by C_j_new 
		for (int i = 0; i < nCell; i++){
			if (c2CellCloneList.get(i) == -1)
				c2CellCloneList.set(i, C_j_new.cloneID);
		}
//		System.out.println(c2CellCloneList);
		// Sample the second clone for cell s_j
		// Compute the distribution
		// P(c2_j = k (not equal to c1_j) | X_j, c1_j, G, Y_j = 1, \alpha, \beta) ~ P(X_j| Y_j = 1, G_c1 + G_c2, \alpha, \beta)
		ArrayList<Double> c2_s_j_unnormalizedDist = new ArrayList<>(); 
		ArrayList<Clone> s_j_c2_otherCloneList = new ArrayList<>();
		for (Clone C: nuCloneList){
			if (C.cloneID != C_j_new.cloneID){
				Clone C_copy = SCF.copyClone(C, nMut);
				s_j_c2_otherCloneList.add(C_copy);
				Integer[] C_copy_doubletGtVector = DoubletGenotype.getDoubletGtVector(C_j_new.cloneGTVector, C_copy.cloneGTVector, nMut, df);
				c2_s_j_unnormalizedDist.add(CLObj.ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellGTVector(C_copy_doubletGtVector, s_j.observedGTVector, _fp, _fnStart, nMut));
			}
		}
		Clone s_j_c2 = SamplingAlgos.sampleInverseCDFDiscrete(SamplingAlgos.normalizeDist(c2_s_j_unnormalizedDist), s_j_c2_otherCloneList);
//		System.out.println("s_j_c2 = "+ s_j_c2.cloneID);
		c2CellCloneList.set(j, s_j_c2.cloneID);
//		System.out.println(c2CellCloneList);
		// Calculate errorLikelihoodRatio
		double errorLikelihoodRatio;
		// cell is not doublet
		if (cellDoubletFlagArr[j] == 0){
			errorLikelihoodRatio = CLObj.ampErrLikelihoodObj.computeNewCloneErrorLikelihoodRatio(C_j_new.cloneGTVector, C_j.cloneGTVector, s_j.observedGTVector, _fp, _fnStart, nMut);
		}
		else{
			Clone s_j_doubletPairClone = SCF.getCloneFrmList(nuCloneList, c2CellCloneList.get(j));
			Integer[] newDoubletCellGtVector = DoubletGenotype.getDoubletGtVector(C_j_new.cloneGTVector, s_j_doubletPairClone.cloneGTVector, nMut, df);
			Integer[] oldDoubletCellGtVector = DoubletGenotype.getDoubletGtVector(C_j.cloneGTVector, s_j_doubletPairClone.cloneGTVector, nMut, df);
			errorLikelihoodRatio = CLObj.ampErrLikelihoodObj.computeNewCloneErrorLikelihoodRatio(newDoubletCellGtVector, oldDoubletCellGtVector, s_j.observedGTVector, _fp, _fnStart, nMut);
		}
		
		// Construct new ListCloneNames and new CloneGenotypeObsMatrix
		ArrayList<String> nuCloneNames = SCF.getNuCloneNameList(nuCloneList);
		ArrayList<GenotypeObservation> nuCloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(nuCloneList, nuCloneNames, nMut);
		
		// The updated objects are stored as static entries
		cloneRemovedCloneList = nuCloneList;
		cloneRemovedCloneNameList = nuCloneNames;
		cloneRemovedCloneGenotypeObsMatrix = nuCloneGenotypeObsMatrix;
		cloneRemovedTree = newTree;
		nuCloneCellj = C_j_new;
		cloneRemovedC2CellCloneList = c2CellCloneList;
		
		double treeLikelihoodRatio = CLObj.treeGenotypeLikelihoodObj.computeTreeCloneLikelihoodRatio(nuCloneGenotypeObsMatrix, cloneGenotypeObsMatrix, newTree, currTree, nMut, df);
//		System.out.println(treeLikelihoodRatio);
		double treePriorRatio = proposeTreeObj.getTreePriorRatio(nuCloneNames.size(), 0) * Math.exp(proposeTreeObj.getClusteringLogPriorRatio(nuCloneList, listClone, alpha_0));
//		System.out.println(treePriorRatio);
		double hastingsRatio = proposeTreeObj.getHastingsRatioRemoveClone(nuCloneNames.size(), nuCloneNames.size()+1);
//		System.out.println("hastings ratio = "+ hastingsRatio);
		double acRatio = PartialMHSteps.getAcceptanceRatio(alpha_0, nCell, errorLikelihoodRatio, treeLikelihoodRatio, treePriorRatio, hastingsRatio, proposeTreeObj.getJacobian(), 0);
//		double acRatio = Math.min(1, ((nCell -1)/alpha_0)*errorLikelihoodRatio);
		return acRatio;
	}

	/**
	 * @param args
	 * Created On: Sep 21, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
