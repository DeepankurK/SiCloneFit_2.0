/**
 * Aug 8, 2017
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
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.CloneGtVector;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Aug 8, 2017
 */
public class SampleCloneGenotype {
	
	public static Random _rng = new Random();
	
	/**
	 * Compute the Posterior Probability of a genotype gt for a clone C for a locus mutIndex
	 * F(G_ki = gt | T) * \prod_{j \in Clone_k} E(X_ji | G_ki = gt) 
	 * @param OtherLeavesGenotypes, The genotypes of other clones for the locus
	 * @param C, clone for which, likelihood is being calculated
	 * @param genotype, the candidate genotype
	 * @param t, Tree
	 * @param mutIndex, index of the locus, 0,...,n-1
	 * @param cellGTVectorMap, The HashMap of observed genotype matrix
	 * @param TGL, TreeGenotypeLikelihood object
	 * @param errorLikelihoodObj, SCSAmplificationErrorLikelihood object
	 * @return
	 * Created On: Aug 8, 2017
	 */
	public static double getCloneGenotypePosterior(GenotypeObservation OtherLeavesGenotypes,
												   Clone C, Integer genotype, Tree t, int mutIndex,
												   HashMap<Integer, Integer[]> cellGTVectorMap,
												   TreeGenotypeLikelihood TGL, 
												   SCSAmplificationErrorLikelihood errorLikelihoodObj){
		
		// This is the likelihood from the tree, F(G_ki | T), working as prior
		double treeGenotypeLogLikelihood = TGL.getSingleLeafGenotypeProbGivenTreeOtherGenotypes(OtherLeavesGenotypes, C.cloneName, genotype);
		// This is the log-likelihood from the error model
		double errorLogLikelihood = errorLikelihoodObj.computeCloneErrorLogLikelihood(genotype, mutIndex, C, cellGTVectorMap);
		double logPosterior = Math.log(treeGenotypeLogLikelihood) + errorLogLikelihood;
		
		return Math.exp(logPosterior);
	}
	
	/**
	 * Compute the Posterior Probability of a genotype gt for a clone C for a locus mutIndex
	 * F(G_ki = gt | T) * \prod_{j \in Clone_k} E(X_ji | G_ki = gt) 
	 * @param OtherLeavesGenotypes
	 * @param C
	 * @param genotype
	 * @param t
	 * @param mutIndex
	 * @param cellGTVectorMap
	 * @param TGL
	 * @param errorLikelihoodObj
	 * @param df
	 * @return
	 * Created On: Oct 2, 2017
	 */
	public static double getCloneGenotypePosterior(GenotypeObservation OtherLeavesGenotypes,
												   Clone C, Integer genotype, Tree t, int mutIndex,
												   HashMap<Integer, Integer[]> cellGTVectorMap,
												   TreeGenotypeLikelihood TGL, 
												   SCSAmplificationErrorLikelihood errorLikelihoodObj, int df){

		// This is the likelihood from the tree, F(G_ki | T), working as prior
		double treeGenotypeLogLikelihood = TGL.getSingleLeafGenotypeProbGivenTreeOtherGenotypes(OtherLeavesGenotypes, C.cloneName, genotype, df);
//		if (Double.isNaN(treeGenotypeLogLikelihood))
//			System.out.println("tree likelihood is NaN");
		// This is the log-likelihood from the error model
		double errorLogLikelihood = errorLikelihoodObj.computeCloneErrorLogLikelihood(genotype, mutIndex, C, cellGTVectorMap);
//		if (Double.isNaN(errorLogLikelihood))
//			System.out.println("error likelihood is NaN");
//		if (mutIndex % 10 == 0)
//			System.out.println("mut = "+ mutIndex + " genotype = " + genotype +  " treeL = " + Math.log(treeGenotypeLogLikelihood) + " errorL = " + errorLogLikelihood);
		double logPosterior = errorLogLikelihood + Math.log(treeGenotypeLogLikelihood) ;

		return Math.exp(logPosterior);
//		return treeGenotypeLogLikelihood;
	}
	
	/**
	 * Compute the Posterior Probability of a genotype gt for a clone C for a locus mutIndex
	 * Considers the presence of doublets
	 * F(G_ki = gt | T,M) * \prod_{j \in Clone_k} E(X_ji | G_ki = gt, c2_j_i, \alpha, \beta) 
	 * @param OtherLeavesGenotypes
	 * @param listClone
	 * @param c2CellCloneList
	 * @param cellDoubletFlagArr
	 * @param C
	 * @param genotype
	 * @param t
	 * @param mutIndex
	 * @param cellGTVectorMap
	 * @param TGL
	 * @param errorLikelihoodObj
	 * @param SCF
	 * @param df
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static double getCloneGenotypePosteriorDoublet(GenotypeObservation OtherLeavesGenotypes, ArrayList<Clone> listClone, 
														   ArrayList<Integer> c2CellCloneList, int[] cellDoubletFlagArr,
														   Clone C, Integer genotype, Tree t, int mutIndex,
														   HashMap<Integer, Integer[]> cellGTVectorMap,
														   TreeGenotypeLikelihood TGL, SCSAmplificationErrorLikelihood errorLikelihoodObj, 
														   SCFUtilityFunctions SCF, int df){
		// This is the likelihood from the tree, F(G_ki | T), working as prior
		double treeGenotypeLogLikelihood = TGL.getSingleLeafGenotypeProbGivenTreeOtherGenotypes(OtherLeavesGenotypes, C.cloneName, genotype, df);
		
		// This is the log-likelihood from the error model
		double errorLogLikelihood = errorLikelihoodObj.computeCloneErrorLogLikelihoodDoublet(genotype, mutIndex, C, listClone, c2CellCloneList, cellDoubletFlagArr, df, cellGTVectorMap, SCF);
		double logPosterior = Math.log(treeGenotypeLogLikelihood) + errorLogLikelihood;

		return Math.exp(logPosterior);
	}
	
	/**
	 * Compute the Posterior Probability for all genotypes for a clone C for a locus (mutIndex)
	 * Considers the doublets
	 * @param OtherLeavesGenotypes
	 * @param listClone
	 * @param c2CellCloneList
	 * @param cellDoubletFlagArr
	 * @param C
	 * @param t
	 * @param mutIndex
	 * @param cellGTVectorMap
	 * @param TGL
	 * @param errorLikelihoodObj
	 * @param SCF
	 * @param df
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static ArrayList<Double> getCloneGenotypePosteriorDoubletList(GenotypeObservation OtherLeavesGenotypes, ArrayList<Clone> listClone, 
																		   ArrayList<Integer> c2CellCloneList, int[] cellDoubletFlagArr,
																		   Clone C, Tree t, int mutIndex,
																		   HashMap<Integer, Integer[]> cellGTVectorMap,
																		   TreeGenotypeLikelihood TGL, SCSAmplificationErrorLikelihood errorLikelihoodObj, 
																		   SCFUtilityFunctions SCF, int df){
		ArrayList<Double> genotypePosteriorList = new ArrayList<>();
		int possibleGenotypes;
		if (df == 0)
			possibleGenotypes = 2;
		else
			possibleGenotypes = 3;
		for (int i = 0; i < possibleGenotypes; i++){
			double genotype_i_posterior = getCloneGenotypePosteriorDoublet(OtherLeavesGenotypes, listClone, c2CellCloneList, cellDoubletFlagArr, C, i, t, mutIndex, cellGTVectorMap, TGL, errorLikelihoodObj, SCF, df);
			genotypePosteriorList.add(genotype_i_posterior);
		}
		return genotypePosteriorList;	
	}
	
	/**
	 * Compute the Posterior Probability for all genotypes for a clone C for a locus (mutIndex)
	 * @param OtherLeavesGenotypes
	 * @param C
	 * @param t
	 * @param mutIndex
	 * @param cellGTVectorMap
	 * @param TGL
	 * @param errorLikelihoodObj
	 * @param df
	 * @return
	 * Created On: Aug 8, 2017
	 */
	public static ArrayList<Double> getCloneGenotypePosteriorList(GenotypeObservation OtherLeavesGenotypes,
												   Clone C, Tree t, int mutIndex,
												   HashMap<Integer, Integer[]> cellGTVectorMap,
												   TreeGenotypeLikelihood TGL, 
												   SCSAmplificationErrorLikelihood errorLikelihoodObj, int df){
		ArrayList<Double> genotypePosteriorList = new ArrayList<>();
		int possibleGenotypes;
		if (df == 0)
			possibleGenotypes = 2;
		else
			possibleGenotypes = 3;
		for (int i = 0; i < possibleGenotypes; i++){
			double genotype_i_posterior = getCloneGenotypePosterior(OtherLeavesGenotypes, C, i, t, mutIndex, cellGTVectorMap, TGL, errorLikelihoodObj, df);
			genotypePosteriorList.add(genotype_i_posterior);
		}
		return genotypePosteriorList;			
	}
	
	/**
	 * Obtain the genotype likelihood distribution (unnormalized) for a single locus
	 * for the new clone being proposed. Probabilities calculated based on the new proposed tree
	 * and given the genotypes of the other clones. 
	 * @param OtherLeavesGenotypes, current genotype of other clones
	 * @param nuClone, new clone whose genotype likelihood distribution is to be computed
	 * @param t, the new tree with the new clone attached to one leaf
	 * @param mutIndex, index of the locus
	 * @param m, evolution model
	 * @param df, data Flag
	 * @return
	 * Created On: Aug 12, 2017
	 */
	public static ArrayList<Double> getSingleMutCloneGenotypeLikelihoodDistFromTree(GenotypeObservation OtherLeavesGenotypes,
			   Clone nuClone, TreeGenotypeLikelihood newTreeLikelihoodObj, int mutIndex, int df){
		ArrayList<Double> genotypeLikelihoodDist = new ArrayList<>();
		int possibleGenotypes;
		if (df == 0)
			possibleGenotypes = 2;
		else
			possibleGenotypes = 3;
		
		for (int i = 0; i < possibleGenotypes; i++){
			double genotype_i_likelihood = newTreeLikelihoodObj.getSingleLeafGenotypeProbGivenTreeOtherGenotypes(OtherLeavesGenotypes, nuClone.cloneName, i);
			genotypeLikelihoodDist.add(genotype_i_likelihood);
		}
		return genotypeLikelihoodDist;
	}
	
	/**
	 * Sample a genotype from the genotypeLikelihoodDist
	 * @param genotypeLikelihoodDist
	 * @return
	 * Created On: Aug 12, 2017
	 */
	public static int sampleSingleMutCloneGenotypeFromTreeLikelihoodDist(ArrayList<Double> genotypeLikelihoodDist){
		// Normalize the distribution
		ArrayList<Double> normalizedGenotypeLikelihoodDist = SamplingAlgos.normalizeDist(genotypeLikelihoodDist);
//		System.out.println(normalizedGenotypeLikelihoodDist);
		// Used Inverse CDF for sampling
		return SamplingAlgos.sampleInverseCDFDiscrete(normalizedGenotypeLikelihoodDist);
	}
	
	public static Clone sampleNewCloneCRPErrLikelihoodDist(ArrayList<Double> CRPErrLikelihood, ArrayList<Clone> listClone){
		ArrayList<Double> normalizedCRPErrLikelihood = SamplingAlgos.normalizeDist(CRPErrLikelihood);
		return SamplingAlgos.sampleInverseCDFDiscrete(normalizedCRPErrLikelihood, listClone);
	}
	
	/**
	 * Sample the Genotype Vector of a clone using the Tree and Error Likelihood
	 * @param OtherLeavesGenotypesList
	 * @param nuClone
	 * @param t
	 * @param m
	 * @param cellGTVectorMap
	 * @param errorLikelihoodObj
	 * @param df
	 * @param nMut
	 * @return
	 * Created On: Aug 16, 2017
	 */
	public static Integer[] sampleCloneGenotypeVectorFromTreeAndErrorLikelihoodDist(
			ArrayList<GenotypeObservation> OtherLeavesGenotypesList, 
			Clone nuClone, Tree t, ComplexEvolutionModel m,
			HashMap<Integer, Integer[]> cellGTVectorMap, 
			SCSAmplificationErrorLikelihood errorLikelihoodObj, int df, int nMut){
		Integer[] cloneGenotypeVector = new Integer[nMut];
		// Construct the tree likelihood object
		TreeGenotypeLikelihood newTreeLikelihoodObj = new TreeGenotypeLikelihood(t, m);
		for (int i = 0; i < nMut; i++){
			ArrayList<Double> genotypeLikelihoodDist_locus_i = getCloneGenotypePosteriorList(OtherLeavesGenotypesList.get(i), nuClone, t, i, cellGTVectorMap, newTreeLikelihoodObj, errorLikelihoodObj, df);
			
//			for (double d : genotypeLikelihoodDist_locus_i){
//				if (Double.isNaN(d))
//					System.out.println("mut " + i + " is NaN");
//			}
			
			cloneGenotypeVector[i] = sampleSingleMutCloneGenotypeFromTreeLikelihoodDist(genotypeLikelihoodDist_locus_i);
		}
		return cloneGenotypeVector;
	}
	
	/**
	 * Sample genotype vector for the new clone from prior distribution from tree
	 * @param OtherLeavesGenotypesList
	 * @param nOtherClones
	 * @param C
	 * @param s
	 * @param t
	 * @param m
	 * @param df
	 * @param nMut
	 * @param sampleIter
	 * @return
	 * Created On: Oct 14, 2017
	 */
	public static Integer[] sampleCloneGenotypeVectorPriorTreeMetropolis(ArrayList<GenotypeObservation> OtherLeavesGenotypesList, int nOtherClones, Clone C, 
																		 SingleCell s, Tree t, ComplexEvolutionModel m, int df, int nMut, int sampleIter){
		// Local lists for storing samples
		ArrayList<Integer[]> sampledGTVectorList = new ArrayList<>();
		int possibleGenotypes;
		if (df == 0)
			possibleGenotypes = 2;
		else
			possibleGenotypes = 3;
		Integer[] cloneGtVector = new Integer[nMut];
		for (int i = 0; i< nMut; i++){
			if (s.observedGTVector[i] != 3)
				cloneGtVector[i] = s.observedGTVector[i];
			else
				cloneGtVector[i] = _rng.nextInt(possibleGenotypes);
		}
		
		// Add current vector to list
		sampledGTVectorList.add(cloneGtVector);

		// Tree-Likelihood Object (Acting as prior here)
		TreeGenotypeLikelihood currTreeLikelihoodObj = new TreeGenotypeLikelihood(t, m);
		for (int st = 1; st < sampleIter; st++){
			Integer[] currGtVector = sampledGTVectorList.get(st - 1);
			CloneGtVector newGtVectorObj = proposeNewGtVector(currGtVector, nMut, df);
			// Calculate prior ratio, comes from the tree
			GenotypeObservation otherClonesGtObservation = OtherLeavesGenotypesList.get(newGtVectorObj.changeGtIndex);
			GenotypeObservation currGtAllClonesGtObs = getFullGenotypeObs(otherClonesGtObservation, nOtherClones, C.cloneName, newGtVectorObj.currGt);
			GenotypeObservation newGtAllClonesGtObs = getFullGenotypeObs(otherClonesGtObservation, nOtherClones, C.cloneName, newGtVectorObj.newGt);
			double newGtVectorLogTreePriorRatio = currTreeLikelihoodObj.getNewGtObsLogLikelihoodRatio(newGtAllClonesGtObs, currGtAllClonesGtObs, df);
			double acRatio = Math.min(1.0, Math.exp(newGtVectorLogTreePriorRatio));
//			System.out.println("iter = " + st + " ac = " + acRatio); 
			double rr1 = _rng.nextDouble();
			if (rr1 <= acRatio){
				sampledGTVectorList.add(newGtVectorObj.gtVector);
			}
			else{
				sampledGTVectorList.add(currGtVector);
			}
		}
		return sampledGTVectorList.get(sampleIter-1);
	}
	
	/**
	 * Sample Genotype vector for a clone from its conditional posterior distribution
	 * Sample complete vector by Metropolis-Hastings algo
	 * @param OtherLeavesGenotypesList
	 * @param C
	 * @param t
	 * @param m
	 * @param cellGTVectorMap
	 * @param errorLikelihoodObj
	 * @param df
	 * @param nMut
	 * @return
	 * Created On: Oct 13, 2017
	 */
	public static Integer[] sampleCloneGenotypeVectorPosteriorMetropolis(ArrayList<GenotypeObservation> OtherLeavesGenotypesList, int nOtherClones, Clone C, 
																		 Tree t, ComplexEvolutionModel m, HashMap<Integer, Integer[]> cellGTVectorMap, 
																		 SCSAmplificationErrorLikelihood errorLikelihoodObj, int df, int nMut, int sampleIter){
		// Local lists for storing samples
		ArrayList<Integer[]> sampledGTVectorList = new ArrayList<>();
		ArrayList<Double> sampledLikelihoodList = new ArrayList<>();
		
		// Add current vector to list
		sampledGTVectorList.add(C.cloneGTVector);
		
		// Tree-Likelihood Object (Acting as prior here)
		TreeGenotypeLikelihood currTreeLikelihoodObj = new TreeGenotypeLikelihood(t, m);
		
		for (int st = 1; st < sampleIter; st++){
			Integer[] currGtVector = sampledGTVectorList.get(st - 1);
			CloneGtVector newGtVectorObj = proposeNewGtVector(currGtVector, nMut, df);
			// Calculate likelihood ratio
			double newVectorLogErrorLikelihoodRatio = errorLikelihoodObj.computeCloneNewGenotypeErrorLogLikelihoodRatio(newGtVectorObj.currGt, newGtVectorObj.newGt, newGtVectorObj.changeGtIndex, C, cellGTVectorMap);
			// Calculate prior ratio, comes from the tree
			GenotypeObservation otherClonesGtObservation = OtherLeavesGenotypesList.get(newGtVectorObj.changeGtIndex);
			GenotypeObservation currGtAllClonesGtObs = getFullGenotypeObs(otherClonesGtObservation, nOtherClones, C.cloneName, newGtVectorObj.currGt);
			GenotypeObservation newGtAllClonesGtObs = getFullGenotypeObs(otherClonesGtObservation, nOtherClones, C.cloneName, newGtVectorObj.newGt);
			double newGtVectorLogTreePriorRatio = currTreeLikelihoodObj.getNewGtObsLogLikelihoodRatio(newGtAllClonesGtObs, currGtAllClonesGtObs, df);
			double acRatio = Math.min(1.0, Math.exp(newVectorLogErrorLikelihoodRatio + newGtVectorLogTreePriorRatio));
//			System.out.println("iter = " + st + " ac = " + acRatio); 
			double rr1 = _rng.nextDouble();
			if (rr1 <= acRatio){
				sampledGTVectorList.add(newGtVectorObj.gtVector);
			}
			else{
				sampledGTVectorList.add(currGtVector);
			}
		}
		
		return sampledGTVectorList.get(sampleIter-1);
	}
	
	/**
	 * Return a full genotype observation from a genotype observation of other clones and genotype of thisClone
	 * @param otherLeaves
	 * @param nOtherClones
	 * @param thisClone
	 * @param genotype
	 * @return
	 * Created On: Oct 13, 2017
	 */
	private static GenotypeObservation getFullGenotypeObs(GenotypeObservation otherLeaves, int nOtherClones, String thisClone, int genotype){
		ArrayList<String> cloneNames = new ArrayList<>();
		Integer[] genotypes = new Integer[nOtherClones+1];
		int i = 0;
		for (Map.Entry<String, Integer> entry: otherLeaves.cellGenotypeMap.entrySet()){
			cloneNames.add(entry.getKey());
			genotypes[i] = entry.getValue();
			i++;
		}
		cloneNames.add(thisClone);
		genotypes[nOtherClones] = genotype;
		GenotypeObservation fullGtObs = new GenotypeObservation(cloneNames, genotypes);
		return fullGtObs;
	}
	
	/**
	 * Propose a new Genotype vector from a current genotype vector.
	 * Change a random position of the genotype
	 * @param currGtVector
	 * @param nMut
	 * @param df
	 * @return
	 * Created On: Oct 13, 2017
	 */
	private static CloneGtVector proposeNewGtVector(Integer[] currGtVector, int nMut, int df){
		Integer[] newGtVector = new Integer[nMut];
		for (int i = 0; i < nMut; i++){
			newGtVector[i] = currGtVector[i];
		}
		int changeGtId = _rng.nextInt(nMut);
		int currGt = currGtVector[changeGtId];
		int newGt = proposeNewgt(currGt, df);
		newGtVector[changeGtId] = newGt;
		CloneGtVector newGtVectorObj = new CloneGtVector(newGtVector, currGt, newGt, changeGtId);
		return newGtVectorObj;
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
	 * Sample the Genotype Vector of a clone using the Tree and Error Likelihood
	 * Considers the presence of doublets
	 * @param OtherLeavesGenotypesList
	 * @param nuClone
	 * @param listClone
	 * @param c2CellCloneList
	 * @param cellDoubletFlagArr
	 * @param t
	 * @param m
	 * @param cellGTVectorMap
	 * @param errorLikelihoodObj
	 * @param df
	 * @param nMut
	 * @param SCF
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static Integer[] sampleCloneGenotypeVectorFromTreeAndErrorLikelihoodDistDoublet(ArrayList<GenotypeObservation> OtherLeavesGenotypesList, Clone nuClone,
																						   ArrayList<Clone> listClone, ArrayList<Integer> c2CellCloneList, int[] cellDoubletFlagArr,
																						   Tree t, ComplexEvolutionModel m, HashMap<Integer, Integer[]> cellGTVectorMap,
																						   SCSAmplificationErrorLikelihood errorLikelihoodObj, int df, int nMut, SCFUtilityFunctions SCF){
		Integer[] cloneGenotypeVector = new Integer[nMut];
		// Construct the tree likelihood object
		TreeGenotypeLikelihood newTreeLikelihoodObj = new TreeGenotypeLikelihood(t, m);
		for (int i = 0; i < nMut; i++){
			ArrayList<Double> genotypeLikelihoodDist_locus_i = getCloneGenotypePosteriorDoubletList(OtherLeavesGenotypesList.get(i), listClone, c2CellCloneList, cellDoubletFlagArr, nuClone, t, i, cellGTVectorMap, newTreeLikelihoodObj, errorLikelihoodObj, SCF, df);
			cloneGenotypeVector[i] = sampleSingleMutCloneGenotypeFromTreeLikelihoodDist(genotypeLikelihoodDist_locus_i);
		}
		return cloneGenotypeVector;
	}

	/**
	 * Sample the genotype vector for the new clone
	 * @param OtherLeavesGenotypesList, List of Genotypes for other clones
	 * @param nuClone, new clone whose genotype vector is to be proposed
	 * @param t, the new tree with the new clone placed on a leaf
	 * @param m, evolution model
	 * @param df, data Flag
	 * @param nMut, number of mutation loci
	 * @return
	 * Created On: Aug 12, 2017
	 */
	public static Integer[] sampleCloneGenotypeVectorFromTreeLikelihoodDist(ArrayList<GenotypeObservation> OtherLeavesGenotypesList,
			Clone nuClone, Tree t, ComplexEvolutionModel m, int df, int nMut){
		Integer[] cloneGenotypeVector = new Integer[nMut];
		// Construct the tree likelihood object
		TreeGenotypeLikelihood newTreeLikelihoodObj = new TreeGenotypeLikelihood(t, m);
		for (int i = 0; i < nMut; i++){
			ArrayList<Double> genotypeLikelihoodDist_locus_i = getSingleMutCloneGenotypeLikelihoodDistFromTree(OtherLeavesGenotypesList.get(i), nuClone, newTreeLikelihoodObj, i, df);
//			System.out.println("locus="+i);
//			System.out.println(genotypeLikelihoodDist_locus_i);
			cloneGenotypeVector[i] = sampleSingleMutCloneGenotypeFromTreeLikelihoodDist(genotypeLikelihoodDist_locus_i);

		}
		return cloneGenotypeVector;
		
	}
	
//	public static Integer[] sampleCloneGenotypeVectorFromErrorLikelihood()
	/**
	 * @param args
	 * Created On: Aug 8, 2017
	 */
	public static void main(String[] args) {
		ArrayList<String> cs = new ArrayList<>();
		for (int i = 0; i < 5; i++){
			cs.add("C"+Integer.toString(i));
		}
		Integer[] gts = new Integer[]{1,0,0,1,1};
		GenotypeObservation oldG = new GenotypeObservation(cs, gts);
		System.out.println(oldG.cellGenotypeMap);
		
		GenotypeObservation newG = SampleCloneGenotype.getFullGenotypeObs(oldG, 5, "C6", 2);
		
		newG.cellGenotypeMap.put("C4", 2);
		System.out.println(newG.cellGenotypeMap);
		System.out.println(oldG.cellGenotypeMap);

	}

}
