/**
 * Sep 14, 2017
 */
package proposal;

import java.util.Random;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * @author hz22
 * Sep 14, 2017
 */
public class ModelParamProposal {
	
	public static Random _rng = new Random();
	
	public static double logBetaPDF(double x, double a, double b){
		BetaDistribution betaObj = new BetaDistribution(a, b);
		return Math.log(betaObj.density(x));
	}
	
	public static double logBetaPDF(double x, BetaDistribution xPrior){
		return Math.log(xPrior.density(x));
	}
	
	public static double getLogPriorRatioModelParam(double newParamVal, double oldParamVal, BetaDistribution paramPrior){
		return logBetaPDF(newParamVal, paramPrior) - logBetaPDF(oldParamVal, paramPrior);
	}
	
	public static double getPriorRatioModelParam(double newParamVal, double oldParamVal, BetaDistribution paramPrior){
		return Math.exp(getLogPriorRatioModelParam(newParamVal, oldParamVal, paramPrior));
	}
	
	public static double getProposalRatio(double currFn, double newFn, double fnSD){
		NormalDistribution forwardNormalDist = new NormalDistribution(currFn, fnSD);
		NormalDistribution backwardNormalDist = new NormalDistribution(newFn, fnSD);
		return Math.log(backwardNormalDist.density(currFn)) - Math.log(forwardNormalDist.density(newFn));
	}
	
	public static double computeAcceptanceRatioFn(double loglikelihoodRatio, double priorRatio, double proposalRatio){
		double r = Math.exp(loglikelihoodRatio + priorRatio + proposalRatio);
		return Math.min(1,r);
	} 
	
	/**
	 * Compute M-H acceptance ratio for a new tree
	 * @param loglikelihoodRatio
	 * @param logPriorRatio
	 * @param proposalRatio
	 * @return
	 * Created On: Sep 18, 2017
	 */
	public static double computeAcceptanceRatioTree(double loglikelihoodRatio, double logPriorRatio, double proposalRatio){
		double r = Math.exp(loglikelihoodRatio + logPriorRatio) * proposalRatio;
		return Math.min(1,r);
	}

	public static boolean getAcceptanceFlag(double acceptance_ratio){
		if (acceptance_ratio == 1){
			return true;
		}
		else {
			double u = _rng.nextDouble();
			if (u < acceptance_ratio){
				return true;
			}
			else{
				return false;
			}
		}
	}
	/**
	 * @param args
	 * Created On: Sep 14, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
