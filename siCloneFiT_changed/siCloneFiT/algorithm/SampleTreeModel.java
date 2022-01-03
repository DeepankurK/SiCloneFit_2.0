/**
 * Sep 18, 2017
 */
package siCloneFiT.algorithm;

import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.distribution.BetaDistribution;

import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import siCloneFiT.likelihood.CompleteLikelihood;
import siCloneFiT.objects.Clone;
import siCloneFiT.proposal.ModelParamProposal;
import siCloneFiT.proposal.TreeProposal;

/**
 * @author hz22
 * Sep 18, 2017
 */
public class SampleTreeModel {
	
	public static Random _rng = new Random();
	
	/**
	 * Sample a new Tree and Model parameters from the conditional distribution
	 * Tree on the same set of leaves is sampled
	 * @param currTree - current Tree
	 * @param CLObj - Complete Likelihood Object
	 * @param delProb - current delProb
	 * @param LOHProb - current LOHProb
	 * @param sampleTreeIter - Number of iterations to run MCMC
	 * @param delProbPriorDist - prior distribution on delProb
	 * @param LOHProbPriorDist - prior distribution on LOHProb
	 * @param cloneGenotypeObsMatrix - Current genotypes of clones
	 * @param nMut - Number of mutations
	 * @param dataFlag 
	 * @param proposeTreeObj
	 * @param nClone
	 * @return
	 * Created On: Sep 18, 2017
	 */
	public static STITree<Clone> sampleCloneTreeMH(STITree<Clone> currTree, CompleteLikelihood CLObj, double delProb, double LOHProb,
										 int sampleTreeIter, BetaDistribution delProbPriorDist, BetaDistribution LOHProbPriorDist,
										 ArrayList<GenotypeObservation> cloneGenotypeObsMatrix, int nMut, int dataFlag, 
										 TreeProposal proposeTreeObj, int nClone){
		
		// Local lists for storing samples
		ArrayList<STITree<Clone>> sampledTreeList = new ArrayList<>();
		ArrayList<Double> sampledDelProbList = new ArrayList<>();
		ArrayList<Double> sampledLOHProbList = new ArrayList<>();
		ArrayList<Double> sampledLikelihoodList = new ArrayList<>();
		
		double currTreeLikelihood = CLObj.treeGenotypeLikelihoodObj.computeTreeCloneLogLikelihood(cloneGenotypeObsMatrix, currTree, nMut, dataFlag);
		
		// Add current tree, model parameters and likelihoods to lists
		sampledLikelihoodList.add(currTreeLikelihood);
		sampledTreeList.add(currTree);
		sampledDelProbList.add(delProb);
		sampledLOHProbList.add(LOHProb);
		
		double bestDelProb = delProb;
		double bestLOHProb = LOHProb;
		STITree<Clone> bestTree = currTree;
		double bestLikelihood = currTreeLikelihood;
		
		
		for (int st = 1; st < sampleTreeIter; st++){
			STITree<Clone> iter_stTree = sampledTreeList.get(st-1);
			double currDelProb = sampledDelProbList.get(st-1);
        	double currLOHProb = sampledLOHProbList.get(st-1);
        	double currLikelihood = sampledLikelihoodList.get(st-1);
        	
        	// Choosing the best
        	if (currLikelihood > bestLikelihood){
        		bestLikelihood = currLikelihood;
        		bestDelProb = currDelProb;
        		bestLOHProb = currLOHProb;
        		bestTree = iter_stTree;
        	}
        	double rr = _rng.nextDouble();
        	if (rr <= 0.5){
        		double rr2 = _rng.nextDouble();
        		// Propose a new value for delProb
        		if (rr2 <= 0.5){
        			double newDelProb = CLObj.treeGenotypeLikelihoodObj.model.proposeNewParamVal(currDelProb, currDelProb*0.1);
        			double newDelProbLogPriorRatio = ModelParamProposal.getLogPriorRatioModelParam(newDelProb, currDelProb, delProbPriorDist);
        			double newDelProbLogProposalRatio = ModelParamProposal.getProposalRatio(currDelProb, newDelProb, currDelProb*0.1);
        			double newDelProbLogLikelihoodRatio = CLObj.treeGenotypeLikelihoodObj.computeTreeCloneNewParamLogLikelihoodRatio(cloneGenotypeObsMatrix, currTree, newDelProb, currLikelihood, nMut, dataFlag, 0);
        			double newDelProbAcceptanceRatio = ModelParamProposal.computeAcceptanceRatioFn(newDelProbLogLikelihoodRatio, newDelProbLogPriorRatio, newDelProbLogProposalRatio);
        			boolean AcceptFlag = ModelParamProposal.getAcceptanceFlag(newDelProbAcceptanceRatio);
        			if (AcceptFlag == true){
        				CLObj.treeGenotypeLikelihoodObj.model.delProb = newDelProb;
        				sampledDelProbList.add(newDelProb);
        				sampledLikelihoodList.add(CLObj.treeGenotypeLikelihoodObj.newParamSameTreeLikelihood);
        			}
        			else{
        				sampledDelProbList.add(currDelProb);
        				sampledLikelihoodList.add(currLikelihood);
        			}
        			sampledLOHProbList.add(currLOHProb);
        			sampledTreeList.add(iter_stTree);
        		}
        		// Propose a new value for LOHProb
        		else{
        			double newLOHProb = CLObj.treeGenotypeLikelihoodObj.model.proposeNewParamVal(currLOHProb, currLOHProb*0.1);
        			double newLOHProbLogPriorRatio = ModelParamProposal.getLogPriorRatioModelParam(newLOHProb, currLOHProb, LOHProbPriorDist);
        			double newLOHProbLogProposalRatio = ModelParamProposal.getProposalRatio(currLOHProb, newLOHProb, currLOHProb*0.1);
        			double newLOHProbLogLikelihoodRatio = CLObj.treeGenotypeLikelihoodObj.computeTreeCloneNewParamLogLikelihoodRatio(cloneGenotypeObsMatrix, currTree, newLOHProb, currLikelihood, nMut, dataFlag, 1);
        			double newLOHProbAcceptanceRatio = ModelParamProposal.computeAcceptanceRatioFn(newLOHProbLogLikelihoodRatio, newLOHProbLogPriorRatio, newLOHProbLogProposalRatio);
        			boolean AcceptFlag = ModelParamProposal.getAcceptanceFlag(newLOHProbAcceptanceRatio);
        			if (AcceptFlag == true){
        				CLObj.treeGenotypeLikelihoodObj.model.LOHProb = newLOHProb;
        				sampledLOHProbList.add(newLOHProb);
        				sampledLikelihoodList.add(CLObj.treeGenotypeLikelihoodObj.newParamSameTreeLikelihood);
        			}
        			else{
        				sampledLOHProbList.add(currLOHProb);
        				sampledLikelihoodList.add(currLikelihood);
        			}
        			sampledDelProbList.add(currDelProb);
        			sampledTreeList.add(iter_stTree);
        		}
        	}
        	else{
        		STITree<Clone> newTree = proposeTreeObj.proposeTree(iter_stTree, nClone);
        		double newTreeLogLikelihoodRatio = CLObj.treeGenotypeLikelihoodObj.computeTreeCloneLogLikelihoodRatio(cloneGenotypeObsMatrix, currTreeLikelihood, newTree, nMut, dataFlag);
        		double newTreeLogPriorRatio = proposeTreeObj.getBranchLengthLogPriorRatio(newTree, iter_stTree);
//        		double newTreeLogLikelihoodRatio = CLObj.treeGenotypeLikelihoodObj.computeTreeCloneLikelihoodRatio(cloneGenotypeObsMatrix, cloneGenotypeObsMatrix, newTree, iter_stTree, nMut, dataFlag);
        		double treeAcceptanceRatio = ModelParamProposal.computeAcceptanceRatioTree(newTreeLogLikelihoodRatio, newTreeLogPriorRatio, proposeTreeObj.getProposalRatio());
        		boolean AcceptFlag = ModelParamProposal.getAcceptanceFlag(treeAcceptanceRatio);
        		if (newTreeLogLikelihoodRatio > 1){
//        		if (AcceptFlag == true){
        			sampledTreeList.add(newTree);
        			sampledLikelihoodList.add(CLObj.treeGenotypeLikelihoodObj.newTreeSameParamLikelihood);
        		}
        		else{
        			sampledTreeList.add(iter_stTree);
        			sampledLikelihoodList.add(currLikelihood);
        		}
        		sampledDelProbList.add(currDelProb);
        		sampledLOHProbList.add(currLOHProb);
        	}
//        	if (sampledLikelihoodList)
		}
		
		CLObj.treeGenotypeLikelihoodObj.model.delProb = bestDelProb;
		CLObj.treeGenotypeLikelihoodObj.model.LOHProb = bestLOHProb;
		return bestTree;
		
//		CLObj.treeGenotypeLikelihoodObj.model.delProb = sampledDelProbList.get(sampleTreeIter-1);
//		CLObj.treeGenotypeLikelihoodObj.model.LOHProb = sampledLOHProbList.get(sampleTreeIter-1);
//		return sampledTreeList.get(sampleTreeIter-1);
	}

	/**
	 * @param args
	 * Created On: Sep 18, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
