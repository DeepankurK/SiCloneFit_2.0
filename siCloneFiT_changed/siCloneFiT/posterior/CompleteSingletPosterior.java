/**
 * Sep 18, 2017
 */
package siCloneFiT.posterior;

import java.util.ArrayList;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.special.Gamma;

import SiFit.model.ComplexEvolutionModel;
import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import siCloneFiT.likelihood.CompleteLikelihood;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.proposal.ModelParamProposal;
import siCloneFiT.proposal.TreeProposal;

/**
 * @author hz22
 * Sep 18, 2017
 */
public class CompleteSingletPosterior {
	
	/**
	 * Compute the Full posterior probability of one sample of the model
	 * @param listSingleCells
	 * @param listClone
	 * @param cloneGenotypeObsMatrix
	 * @param CLObj
	 * @param proposeTreeObj
	 * @param tree
	 * @param m
	 * @param nMut
	 * @param nCell
	 * @param Fn
	 * @param FnPriorDist
	 * @param Fp
	 * @param FpPriorDist
	 * @param delProb
	 * @param delProbPriorDist
	 * @param LOHProb
	 * @param LOHProbPriorDist
	 * @param alpha_0
	 * @param alpha_0PriorDist
	 * @return
	 * Created On: Sep 18, 2017
	 */
	public static double computeFullModelLogPosterior(ArrayList<SingleCell> listSingleCells,ArrayList<SingleCell> listSingleCells_tot,												  ArrayList<Clone> listClone,
													  ArrayList<GenotypeObservation> cloneGenotypeObsMatrix,
													  CompleteLikelihood CLObj,
													  TreeProposal proposeTreeObj,
													  Tree tree, ComplexEvolutionModel m,
													  int nMut, int nCell, int df,
													  ArrayList<Double> Fn, BetaDistribution FnPriorDist,
													  double Fp, BetaDistribution FpPriorDist,
													  double delProb, BetaDistribution delProbPriorDist,
													  double LOHProb, BetaDistribution LOHProbPriorDist,
													  double alpha_0, GammaDistribution alpha_0PriorDist){
		double logLikelihood = CLObj.computeFullLikelihood(listSingleCells, listSingleCells_tot, listClone, tree, m, Fp, Fn, cloneGenotypeObsMatrix, nMut, df);
//		double logLikelihood = CLObj.computeFullLikelihood(listSingleCells, listClone, tree, m, Fp, Fn, cloneGenotypeObsMatrix, nMut);
		double partition = computePartitionLogProb(listClone, alpha_0, nCell);
//		double partition = 0;
		double treePrior = proposeTreeObj.getTreeLogPriorProb(tree);
//		double treePrior = proposeTreeObj.getTreeLogPriorProb(tree, nCell);
//		double treePrior = 0;
		double paramValPrior = getParameterLogPrior(Fn, FnPriorDist, Fp, FpPriorDist, delProb, delProbPriorDist, LOHProb, LOHProbPriorDist, alpha_0, alpha_0PriorDist);
		return logLikelihood + treePrior + partition + paramValPrior;
	}
	
	/**
	 * Compute the log of prior for the parameter values
	 * @param Fn
	 * @param FnPriorDist
	 * @param Fp
	 * @param FpPriorDist
	 * @param delProb
	 * @param delProbPriorDist
	 * @param LOHProb
	 * @param LOHProbPriorDist
	 * @param alpha_0
	 * @param alpha_0PriorDist
	 * @return
	 * Created On: Sep 18, 2017
	 */
	public static double getParameterLogPrior(ArrayList<Double> Fn, BetaDistribution FnPriorDist,
											  double Fp, BetaDistribution FpPriorDist,
											  double delProb, BetaDistribution delProbPriorDist,
											  double LOHProb, BetaDistribution LOHProbPriorDist,
											  double alpha_0, GammaDistribution alpha_0PriorDist){
		double sum= ModelParamProposal.logBetaPDF(Fp, FpPriorDist) + 
			   ModelParamProposal.logBetaPDF(delProb, delProbPriorDist) +
			   ModelParamProposal.logBetaPDF(LOHProb, LOHProbPriorDist) +
			   Math.log(alpha_0PriorDist.density(alpha_0));
		for(Double beta: Fn){
			sum+=ModelParamProposal.logBetaPDF(beta, FnPriorDist);
		}
		return sum;
	}
	
	/**
	 * Compute the probability of a partition of the cells into clones
	 * log P(C | alpha_0) = log(Gamma(alpha_0)) + K*log(alpha_0) + \sum_j = 1^K log (Gamma(m_j)) - log(Gamma(alpha_0 + m))
	 * @param listClone, current partition of cells
	 * @param alpha_0 - concentration parameter
	 * @param nCell - m in above formula, the number of cells
	 * @return
	 * Created On: Sep 18, 2017
	 */
	public static double computePartitionLogProb(ArrayList<Clone> listClone, double alpha_0, int nCell){
		int nClone = listClone.size();
		double logProb = Gamma.logGamma(alpha_0) + nClone * Math.log(alpha_0) - Gamma.logGamma(alpha_0 + nCell);
		for (Clone C: listClone){
			logProb += Gamma.logGamma(C.memberCellList.size());
		}
		return logProb;
		
	}
	
	public static double computePartitionLogProbTest(ArrayList<Integer> listClone, double alpha_0, int nCell){
		int nClone = listClone.size();
		double logProb = Gamma.logGamma(alpha_0) + nClone * Math.log(alpha_0) - Gamma.logGamma(alpha_0 + nCell);
		for (Integer i: listClone)
			logProb += Gamma.logGamma(i);
		return logProb;
	}

	/**
	 * @param args
	 * Created On: Sep 18, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ArrayList<Integer> lC = new ArrayList<>();
		lC.add(2);
		lC.add(3);
		lC.add(4);
		lC.add(1);
		System.out.println(CompleteSingletPosterior.computePartitionLogProbTest(lC, 3, 10));
		ArrayList<Integer> nC = new ArrayList<>();
		nC.add(3);
		nC.add(3);
		nC.add(3);
		nC.add(1);
		System.out.println(CompleteSingletPosterior.computePartitionLogProbTest(nC, 3, 10));
	}

}
