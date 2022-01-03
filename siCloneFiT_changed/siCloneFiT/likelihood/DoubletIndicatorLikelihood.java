/**
 * Oct 9, 2017
 */
package siCloneFiT.likelihood;

/**
 * @author hz22
 * Oct 9, 2017
 */
public class DoubletIndicatorLikelihood {
	
	/**
	 * Compute the likelihood of observing a doublet indicator flag given delta
	 * @param cellDoubletFlagArr
	 * @param delta
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static double getDoubletFlagVectorLikelihood(int[] cellDoubletFlagArr, double delta){
		double logLikelihood = 0.0;
		for (int i = 0; i < cellDoubletFlagArr.length; i++){
			logLikelihood += getDoubletFlagLikelihood(cellDoubletFlagArr[i], delta);
		}
		return Math.exp(logLikelihood);
	}
	
	/**
	 * Compute pdf of Bernoulli Distribution on delta
	 * if 1, return delta, else return 1 - delta
	 * @param flag
	 * @param delta
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static double getDoubletFlagLikelihood(int flag, double delta){
		if (flag == 0)
			return Math.log(1 - delta);
		else
			return Math.log(delta);
	}

	/**
	 * @param args
	 * Created On: Oct 9, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
