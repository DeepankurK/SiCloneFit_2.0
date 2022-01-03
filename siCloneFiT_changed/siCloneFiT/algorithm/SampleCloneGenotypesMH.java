/**
 * Oct 14, 2017
 */
package siCloneFiT.algorithm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import SiFit.model.ComplexEvolutionModel;
import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import siCloneFiT.likelihood.SCSAmplificationErrorLikelihood;
import siCloneFiT.likelihood.TreeGenotypeLikelihood;
import siCloneFiT.objects.ChangeCloneGtMatrix;
import siCloneFiT.objects.Clone;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Oct 14, 2017
 */
public class SampleCloneGenotypesMH {
	
	public static Random _rng = new Random();
	
	/**
	 * Sample the genotype of all clones together using metropolis-hastings
	 * @param listClone
	 * @param cloneNames
	 * @param nClone
	 * @param t
	 * @param m
	 * @param cellIDGTVectorMap
	 * @param errorLikelihoodObj
	 * @param df
	 * @param nMut
	 * @param sampleIter
	 * @param SCF
	 * @return
	 * Created On: Oct 14, 2017
	 */
	public static ArrayList<Integer[]> sampleAllCloneGenotypesPosteriorMetropolis(ArrayList<Clone> listClone, ArrayList<String> cloneNames, int nClone,
																				  Tree t, ComplexEvolutionModel m, HashMap<Integer, Integer[]> cellIDGTVectorMap, 
																				  SCSAmplificationErrorLikelihood errorLikelihoodObj, int df, int nMut, int sampleIter,
																				  SCFUtilityFunctions SCF){
		ArrayList<Integer[]> cloneGenotypeMatrix = getCopyCloneGenotypeMatrix(listClone, nMut);
		// Local list for storing samples
		ArrayList<ArrayList<Integer[]>> sampledCloneGtMatList = new ArrayList<>();
		ArrayList<ArrayList<GenotypeObservation>> sampledCloneGenotypeObsMatrix = new ArrayList<>();
		
		
		// Tree-Likelihood Object (Acting as prior here)
		TreeGenotypeLikelihood currTreeLikelihoodObj = new TreeGenotypeLikelihood(t, m);
		
		// Obtain a local copy of cloneGenotypeObs
		ArrayList<GenotypeObservation> cloneGenotypeObsMatrix = SCF.getCloneGenotypeObs(listClone, cloneNames, nMut);
		
		// Add current matrix to list
		sampledCloneGtMatList.add(cloneGenotypeMatrix);
		sampledCloneGenotypeObsMatrix.add(cloneGenotypeObsMatrix);
		
		for (int st = 1; st < sampleIter; st++){
			ArrayList<Integer[]> currCloneGenotypeMatrix = sampledCloneGtMatList.get(st - 1);
//			double currLikelihood = errorLikelihoodObj.computeAmpErrorLoglikelihoodCellMatrix(listClone, currCloneGenotypeMatrix, cellIDGTVectorMap, nMut);
//			System.out.println("likelihood = " + currLikelihood);
			ArrayList<GenotypeObservation> currCloneGenotypeObsMatrix = sampledCloneGenotypeObsMatrix.get(st-1);
			ChangeCloneGtMatrix newGtMatObj = proposeNewCloneGtMatSingleChange(currCloneGenotypeMatrix, nMut, df, cloneNames);
			Clone cloneToChange = listClone.get(newGtMatObj.changeCloneID);
			// Calculate likelihood ratio
			double newVectorLogErrorLikelihoodRatio = errorLikelihoodObj.computeCloneNewGenotypeErrorLogLikelihoodRatio(newGtMatObj.currGt, newGtMatObj.newGt, newGtMatObj.changeGtIndex, cloneToChange, cellIDGTVectorMap);
			
			// Calculate prior ratio, comes from the tree
			GenotypeObservation currGtAllClonesGtObs = currCloneGenotypeObsMatrix.get(newGtMatObj.changeGtIndex);
			GenotypeObservation newGtAllClonesGtObs = getNewGenotypeObs(currGtAllClonesGtObs, newGtMatObj.changeCloneName, newGtMatObj.newGt, nClone);
			double newGtVectorLogTreePriorRatio = currTreeLikelihoodObj.getNewGtObsLogLikelihoodRatio(newGtAllClonesGtObs, currGtAllClonesGtObs, df);
			double acRatio = Math.min(1.0, Math.exp(newVectorLogErrorLikelihoodRatio + newGtVectorLogTreePriorRatio));
			double rr1 = _rng.nextDouble();
			if (rr1 <= acRatio){
//				System.out.println("new genotype selected");
				sampledCloneGtMatList.add(newGtMatObj.newCloneGtMatrix);
				currCloneGenotypeObsMatrix.set(newGtMatObj.changeGtIndex, newGtAllClonesGtObs);
				sampledCloneGenotypeObsMatrix.add(currCloneGenotypeObsMatrix);
			}
			else{
				sampledCloneGtMatList.add(currCloneGenotypeMatrix);
				sampledCloneGenotypeObsMatrix.add(currCloneGenotypeObsMatrix);
			}
		}
		
		
		
		
		return sampledCloneGtMatList.get(sampleIter-1);
	}
	
	/**
	 * Get new genotype observation
	 * @param currGtObs
	 * @param cloneName
	 * @param genotype
	 * @param nClone
	 * @return
	 * Created On: Oct 14, 2017
	 */
	private static GenotypeObservation getNewGenotypeObs(GenotypeObservation currGtObs, String cloneName, int genotype, int nClone){
		ArrayList<String> cloneNames = new ArrayList<>();
		Integer[] genotypes = new Integer[nClone];
		int i = 0;
		for (Map.Entry<String, Integer> entry: currGtObs.cellGenotypeMap.entrySet()){
			cloneNames.add(entry.getKey());
			if (entry.getKey() == cloneName)
				genotypes[i] = genotype;
			else
				genotypes[i] = entry.getValue();
			i++;
		}
		GenotypeObservation newGtObs = new GenotypeObservation(cloneNames, genotypes);
		return newGtObs;
	}
	
	/**
	 * Create a new Clone genotype matrix by changing a single mutation of a clone, all of them randomly chosen
	 * @param currGtMat
	 * @param nMut
	 * @param df
	 * @param cloneNames
	 * @return
	 * Created On: Oct 14, 2017
	 */
	public static ChangeCloneGtMatrix proposeNewCloneGtMatSingleChange(ArrayList<Integer[]> currGtMat, int nMut, int df, ArrayList<String> cloneNames){
		ArrayList<Integer[]> newCloneGtMat = copyCloneGenotypeMatrix(currGtMat, nMut);
		int changeCloneID = _rng.nextInt(cloneNames.size());
		String changeCloneName = cloneNames.get(changeCloneID);
		int changeGtId = _rng.nextInt(nMut);
		int currGt = currGtMat.get(changeCloneID)[changeGtId];
		int newGt = proposeNewgt(currGt, df);
		newCloneGtMat.get(changeCloneID)[changeGtId] = newGt;
		ChangeCloneGtMatrix newGtMatObj = new ChangeCloneGtMatrix(newCloneGtMat, currGt, newGt, changeGtId, changeCloneName, changeCloneID);
		return newGtMatObj;
	}
	
	/**
	 * Copy a clone genotype matrix
	 * @param currGtMat
	 * @param nMut
	 * @return
	 * Created On: Oct 14, 2017
	 */
	private static ArrayList<Integer[]> copyCloneGenotypeMatrix(ArrayList<Integer[]> currGtMat, int nMut){
		ArrayList<Integer[]> cloneGenotypeMatrix = new ArrayList<>();
		for (Integer[] C: currGtMat){
			Integer[] C_vector = new Integer[nMut];
			for (int i = 0; i < nMut; i++){
				C_vector[i] = C[i];
			}
			cloneGenotypeMatrix.add(C_vector);
		}
		return cloneGenotypeMatrix;
	}
	
	/**
	 * Get a clone Genotype matrix from list of clones 
	 * @param listClone
	 * @param cloneNames
	 * @param nMut
	 * @return
	 * Created On: Oct 14, 2017
	 */
	private static ArrayList<Integer[]> getCopyCloneGenotypeMatrix(ArrayList<Clone> listClone, int nMut){
		ArrayList<Integer[]> cloneGenotypeMatrix = new ArrayList<>();
		for (Clone C: listClone){
			Integer[] C_vector = new Integer[nMut];
			for (int i = 0; i < nMut; i++){
				C_vector[i] = C.cloneGTVector[i];
			}
			cloneGenotypeMatrix.add(C_vector);
		}
		return cloneGenotypeMatrix;
	}
	
	/**
	 * Propose a new genotype for a given current genotype
	 * @param currGt
	 * @param df
	 * @return
	 * Created On: Oct 13, 2017
	 */
	private static int proposeNewgt(int currGt, int df){
		// binary genotype
		if (df == 0){
			if (currGt == 0)
				return 1;
			else
				return 0;
		}
		// ternary genotype
		else{
			if (currGt == 0){
				return 1 + _rng.nextInt(2);
			}
			else if (currGt == 2){
				return _rng.nextInt(2);
			}
			else{
				double rr = _rng.nextDouble();
				if (rr <= 0.5)
					return 0;
				else
					return 1;
			}
				
		}
	}

	/**
	 * @param args
	 * Created On: Oct 14, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
