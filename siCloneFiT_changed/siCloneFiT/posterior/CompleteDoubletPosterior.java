/**
 * Oct 30, 2017
 */
package siCloneFiT.posterior;

import java.util.ArrayList;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;

import SiFit.model.ComplexEvolutionModel;
import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import siCloneFiT.likelihood.CompleteLikelihood;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.proposal.ModelParamProposal;
import siCloneFiT.proposal.TreeProposal;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Oct 30, 2017
 */
public class CompleteDoubletPosterior extends CompleteSingletPosterior {
	
	/**
	 * Compute the Full posterior probability of one sample of the doublet model
	 * @param listSingleCells
	 * @param cellDoubletFlagArr
	 * @param listClone
	 * @param c2CellCloneList
	 * @param cloneGenotypeObsMatrix
	 * @param CLObj
	 * @param proposeTreeObj
	 * @param tree
	 * @param m
	 * @param nMut
	 * @param nCell
	 * @param df
	 * @param Fn
	 * @param FnPriorDist
	 * @param Fp
	 * @param FpPriorDist
	 * @param _doublet
	 * @param doubletPriorDist
	 * @param delProb
	 * @param delProbPriorDist
	 * @param LOHProb
	 * @param LOHProbPriorDist
	 * @param alpha_0
	 * @param alpha_0PriorDist
	 * @param SCF
	 * @return
	 * Created On: Oct 30, 2017
	 */
	public static double computeFullModelLogPosterior(ArrayList<SingleCell> listSingleCells, int[] cellDoubletFlagArr,
													  ArrayList<Clone> listClone, ArrayList<Integer> c2CellCloneList, 
													  ArrayList<GenotypeObservation> cloneGenotypeObsMatrix,
													  CompleteLikelihood CLObj,
													  TreeProposal proposeTreeObj,
													  Tree tree, ComplexEvolutionModel m,
													  int nMut, int nCell, int df,
													  double Fn, BetaDistribution FnPriorDist,
													  double Fp, BetaDistribution FpPriorDist,
													  double _doublet, BetaDistribution doubletPriorDist,
													  double delProb, BetaDistribution delProbPriorDist,
													  double LOHProb, BetaDistribution LOHProbPriorDist,
													  double alpha_0, GammaDistribution alpha_0PriorDist,
													  SCFUtilityFunctions SCF){
		double logLikelihood = CLObj.computeFullLikelihoodDoublet(listSingleCells, cellDoubletFlagArr, listClone, c2CellCloneList, cloneGenotypeObsMatrix, tree, m, SCF, Fp, Fn, nMut, df);
		double partition = computePartitionLogProb(listClone, alpha_0, nCell);
		double treePrior = proposeTreeObj.getTreeLogPriorProb(tree);
		double paramValPrior = getParameterLogPrior(Fn, FnPriorDist, Fp, FpPriorDist, _doublet, doubletPriorDist, delProb, delProbPriorDist, LOHProb, LOHProbPriorDist, alpha_0, alpha_0PriorDist);
		return logLikelihood + treePrior + partition + paramValPrior;
	}
	/**
	 * Compute the log of prior for the parameter values
	 * @param Fn
	 * @param FnPriorDist
	 * @param Fp
	 * @param FpPriorDist
	 * @param _doublet
	 * @param doubletPriorDist
	 * @param delProb
	 * @param delProbPriorDist
	 * @param LOHProb
	 * @param LOHProbPriorDist
	 * @param alpha_0
	 * @param alpha_0PriorDist
	 * @return
	 * Created On: Oct 30, 2017
	 */
	public static double getParameterLogPrior(double Fn, BetaDistribution FnPriorDist,
												double Fp, BetaDistribution FpPriorDist,
												double _doublet, BetaDistribution doubletPriorDist,
												double delProb, BetaDistribution delProbPriorDist,
												double LOHProb, BetaDistribution LOHProbPriorDist,
												double alpha_0, GammaDistribution alpha_0PriorDist){
		return ModelParamProposal.logBetaPDF(Fn, FnPriorDist) + 
				ModelParamProposal.logBetaPDF(Fp, FpPriorDist) +
				ModelParamProposal.logBetaPDF(_doublet, doubletPriorDist) +
				ModelParamProposal.logBetaPDF(delProb, delProbPriorDist) +
				ModelParamProposal.logBetaPDF(LOHProb, LOHProbPriorDist) +
				Math.log(alpha_0PriorDist.density(alpha_0));
	}

	/**
	 * @param args
	 * Created On: Oct 30, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
