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
public class SampleIndicator2Doublet {
	
	/**
	 * Sample the clone of the doublet mate of the cell
	 * sample c2 for a cell j, currently only based on likelihood
	 * Samples only from the other clones (other than c1_j)
	 * c2_j = k (other than c1_j) | c1_j, X_j, G, \alpha, \beta, Y_j = 1  ~  P(X_j | G_c1_j + G_c2_j, \alpha, \beta, Y_j = 1)
	 * @param s_j
	 * @param C_j
	 * @param listClone
	 * @param currFp
	 * @param currFn
	 * @param nMut
	 * @param dataFlag
	 * @param errorLikelihoodObj
	 * @param SCF
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static Clone sampleDoubletMateClone(SingleCell s_j, Clone C_j, ArrayList<Clone> listClone, double currFp, double currFn, int nMut, 
											   int dataFlag, SCSAmplificationErrorLikelihood errorLikelihoodObj, SCFUtilityFunctions SCF){
		ArrayList<Clone> s_j_OtherCloneList = new ArrayList<>();
		ArrayList<Double> s_j_c2_unnormalizedDist = new ArrayList<>();
		for (Clone C: listClone){
			if (C.cloneID != C_j.cloneID){
				Clone C_copy = SCF.copyClone(C, nMut);
				s_j_OtherCloneList.add(C_copy);
				Integer[] C_copy_doubletGtVector = DoubletGenotype.getDoubletGtVector(C_j.cloneGTVector, C_copy.cloneGTVector, nMut, dataFlag);
				double val = errorLikelihoodObj.computeAmpErrorLogLikelihoodCellGTVector(C_copy_doubletGtVector, s_j.observedGTVector, currFp, currFn, nMut);
				s_j_c2_unnormalizedDist.add(val);
			}
		}
		Clone s_j_c2 = SamplingAlgos.sampleInverseCDFDiscrete(SamplingAlgos.normalizeDist(s_j_c2_unnormalizedDist), s_j_OtherCloneList);
		return s_j_c2;
	}

	/**
	 * @param args
	 * Created On: Oct 9, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
