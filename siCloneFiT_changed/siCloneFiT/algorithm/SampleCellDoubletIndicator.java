/**
 * Oct 9, 2017
 */
package siCloneFiT.algorithm;

import java.util.ArrayList;

import siCloneFiT.likelihood.SCSAmplificationErrorLikelihood;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.DoubletGenotype;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Oct 9, 2017
 */
public class SampleCellDoubletIndicator {
	
	/**
	 * Compute the conditional posterior of a cell being a doublet
	 * Y_j | X_j, c1_j, c2_j, G, \alpha, \beta, \delta  ~  P(X_j | Y_j, c1_j, c2_j, G, \alpha, \beta) * P(Y_j | \delta)
	 * if Y_j = 0, it is a singlet
	 * if Y_j = 1, it is a doublet
	 * @param doubletFlag
	 * @param s
	 * @param sClone
	 * @param doubletMateClone
	 * @param Fp
	 * @param Fn
	 * @param delta
	 * @param errorLikelihoodObj
	 * @param SCF
	 * @param df
	 * @param nMut
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static double getCellDoubletIndicatorPosterior(int doubletFlag, SingleCell s, Clone sClone, Clone doubletMateClone, double Fp, double Fn, double delta,
														  SCSAmplificationErrorLikelihood errorLikelihoodObj, int df, int nMut){
		double prior;
		double likelihood;
		// Cell is not doublet
		if (doubletFlag == 0){
			prior = 1 - delta;
			likelihood = errorLikelihoodObj.computeAmpErrorLogLikelihoodCellGTVector(sClone.cloneGTVector, s.observedGTVector, Fp, Fn, nMut);
		}
		// Cell is doublet
		else{
			prior = delta;
			Integer[] cellDoubletGtVector = DoubletGenotype.getDoubletGtVector(sClone.cloneGTVector, doubletMateClone.cloneGTVector, nMut, df);
			likelihood = errorLikelihoodObj.computeAmpErrorLogLikelihoodCellGTVector(cellDoubletGtVector, s.observedGTVector, Fp, Fn, nMut);
		}
		return prior * Math.exp(likelihood);
	}
	
	/**
	 * Sample a cellDoubletFlag from the posterior 
	 * @param s
	 * @param sClone
	 * @param doubletMateClone
	 * @param Fp
	 * @param Fn
	 * @param delta
	 * @param errorLikelihoodObj
	 * @param df
	 * @param nMut
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static int sampleCellDoubletIndicator(SingleCell s, Clone sClone, Clone doubletMateClone, double Fp, double Fn, double delta,
												SCSAmplificationErrorLikelihood errorLikelihoodObj, int df, int nMut){ 
		ArrayList<Double> cellDoubletIndicatorPosteriorList = new ArrayList<>();
		for (int i = 0; i < 2; i++){
			double val = getCellDoubletIndicatorPosterior(i, s, sClone, doubletMateClone, Fp, Fn, delta, errorLikelihoodObj, df, nMut);
			cellDoubletIndicatorPosteriorList.add(val);
		}
		ArrayList<Double> normalizedPosteriorList = SamplingAlgos.normalizeDist(cellDoubletIndicatorPosteriorList);
		return SamplingAlgos.sampleInverseCDFDiscrete(normalizedPosteriorList);
	}

	/**
	 * @param args
	 * Created On: Oct 9, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
