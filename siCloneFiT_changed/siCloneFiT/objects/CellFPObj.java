/**
 * Oct 2, 2017
 */
package siCloneFiT.objects;

import java.util.HashMap;
import java.util.Set;

import SiFit.objects.AmpFPObj;

/**
 * @author hz22
 * Oct 2, 2017
 */
public class CellFPObj extends AmpFPObj {

	/**
	 * @param cellGenotypeArr
	 * @param cellName
	 */
	public CellFPObj(Integer[] cellGenotypeArr, String cellName) {
		super(cellGenotypeArr, cellName);
		// TODO Auto-generated constructor stub
	}
	
	public void getFPAffectedGenotypeArr(double FP_rate, Integer[] beforeFPGenotypeArr, HashMap<Integer, Set<Integer>> Genotype_flag_array, int df){
		this.afterFPGenotypeArr = new Integer[beforeFPGenotypeArr.length];
		if (df == 1){
			for (int i = 0; i < beforeFPGenotypeArr.length; i++){
				afterFPGenotypeArr[i] = beforeFPGenotypeArr[i];
				if (beforeFPGenotypeArr[i] == 0){
					double r_prob = _rng.nextDouble();
					if (r_prob < FP_rate){
						afterFPGenotypeArr[i] = 1;
						Genotype_flag_array.get(1).add(i);
						FPAffectedPosList.add(i);
					}
				}
				else if (beforeFPGenotypeArr[i] == 1){
					double r_prob = _rng.nextDouble();
					if (r_prob < FP_rate){
						afterFPGenotypeArr[i] = 2;
						Genotype_flag_array.get(2).add(i);
						FPAffectedPosList.add(i);
					}
				}
			}
		}
		else{
			for (int i = 0; i < beforeFPGenotypeArr.length; i++){
				afterFPGenotypeArr[i] = beforeFPGenotypeArr[i];
				if (beforeFPGenotypeArr[i] == 0){
					double r_prob = _rng.nextDouble();
					if (r_prob < FP_rate){
						afterFPGenotypeArr[i] = 1;
						Genotype_flag_array.get(1).add(i);
						FPAffectedPosList.add(i);
					}
				}
//				else if (beforeFPGenotypeArr[i] == 1){
//					double r_prob = _rng.nextDouble();
//					if (r_prob < FP_rate){
//						afterFPGenotypeArr[i] = 2;
//						Genotype_flag_array.get(2).add(i);
//						FPAffectedPosList.add(i);
//					}
//				}
			}
		}
	}

	/**
	 * @param args
	 * Created On: Oct 2, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
