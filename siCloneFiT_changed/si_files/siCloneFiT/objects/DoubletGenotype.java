/**
 * Oct 8, 2017
 */
package objects;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author hz22
 * Oct 8, 2017
 */
public class DoubletGenotype {
	
	/**
	 * Compute the merged genotype vector for a doublet
	 * @param gt1
	 * @param gt2
	 * @param nMut
	 * @param dataFlag
	 * @return
	 * Created On: Oct 8, 2017
	 */
	public static Integer[] getDoubletGtVector(Integer[] gt1, Integer[] gt2, int nMut, int dataFlag){
		Integer[] doubletGtVector = new Integer[nMut];
		if (dataFlag == 0){
			for (int i = 0; i < nMut; i++){
				doubletGtVector[i] = gt1[i] | gt2[i];
			}
		}
		else{
			for (int i = 0; i < nMut; i++){
				doubletGtVector[i] = getDoubletGt(gt1[i], gt2[i]);
			}
		}
		return doubletGtVector;
	}
	
	/**
	 * Return doublet merged genotype for a single position
	 * @param gt1
	 * @param gt2
	 * @param df
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static int getDoubletGt(Integer gt1, Integer gt2, int df){
		if (df == 0)
			return gt1 | gt2;
		else
			return getDoubletGt(gt1, gt2);
	}

	/**
	 * Return doublet merged genotype for ternary settings
	 * @param gt1
	 * @param gt2
	 * @return
	 * Created On: Oct 8, 2017
	 */
	private static Integer getDoubletGt(Integer gt1, Integer gt2) {
		if (gt1 == 0 & gt2 == 0)
			return 0;
		else if (gt1 == 2 & gt2 == 2)
			return 2;
		else
			return 1;
	}

	/**
	 * @param args
	 * Created On: Oct 8, 2017
	 */
	public static void main(String[] args) {
		Integer[] gt1 = new Integer[]{1,0,0,0,2};
		Integer[] gt2 = new Integer[]{2,1,1,0,2};
		Integer[] doublet = getDoubletGtVector(gt1, gt2, 5, 1);
		System.out.println(Arrays.toString(doublet));

		ArrayList<String> ss = new ArrayList<>();
		ss.add("s1");
		ss.add("s2");
		ss.add("s3");
		ss.set(1, "2s");
		System.out.println(ss);
		
	}

}
