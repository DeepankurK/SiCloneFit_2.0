/**
 * Oct 2, 2017
 */
package objects;

import SiFit.objects.AmpADOObj;

/**
 * @author hz22
 * Oct 2, 2017
 */
public class CellADOObj extends AmpADOObj {

	/**
	 * @param cellName
	 */
	public CellADOObj(String cellName) {
		super(cellName);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param cellGenotypeArr
	 * @param cellName
	 */
	public CellADOObj(Integer[] cellGenotypeArr, String cellName) {
		super(cellGenotypeArr, cellName);
		// TODO Auto-generated constructor stub
	}

	public void getADOAffectedGenotypeArr(double ADO_rate, Integer[] beforeADOGenotypeArr, int dataFlag){
		this.afterADOGenotypeArr =  new Integer[beforeADOGenotypeArr.length];
		
		if (dataFlag == 1){
			for (int i = 0; i < beforeADOGenotypeArr.length; i++){
				afterADOGenotypeArr[i] = beforeADOGenotypeArr[i];
				if (beforeADOGenotypeArr[i] == 1){
					double r_prob = _rng.nextDouble();
					// ADO happening here
					if (r_prob < ADO_rate){
						double refORalt = _rng.nextDouble();
						// Alt allele is dropped
						if (refORalt <= 0.5){
							afterADOGenotypeArr[i] = 0;
						}
						// Ref allele is dropped
						else{
							afterADOGenotypeArr[i] = 2;
						}
						ADOAffectedPosList.add(i);
					}
				}
			}
		}
		else{
			for (int i = 0; i < beforeADOGenotypeArr.length; i++){
				afterADOGenotypeArr[i] = beforeADOGenotypeArr[i];
				if (beforeADOGenotypeArr[i] == 1){
					double r_prob = _rng.nextDouble();
					// ADO happening here
					if (r_prob < ADO_rate){
						double refORalt = _rng.nextDouble();
						// Alt allele is dropped
						if (refORalt <= 0.5){
							afterADOGenotypeArr[i] = 0;
						}
						// Ref allele is dropped
						else{
							afterADOGenotypeArr[i] = 1;
						}
						ADOAffectedPosList.add(i);
					}
				}
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
