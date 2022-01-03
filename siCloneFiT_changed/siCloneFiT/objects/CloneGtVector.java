/**
 * Oct 13, 2017
 */
package siCloneFiT.objects;

import java.util.Arrays;

/**
 * @author hz22
 * Oct 13, 2017
 */
public class CloneGtVector {
	
	public Integer[] gtVector;
	public int currGt;
	public int newGt;
	public int changeGtIndex;
	
	public CloneGtVector(Integer[] newGtVector, int cg, int ng, int index){
		this.gtVector = newGtVector;
		this.currGt = cg;
		this.newGt = ng;
		this.changeGtIndex = index;
	}

	/**
	 * @param args
	 * Created On: Oct 13, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		Integer[] clone = new Integer[]{1,0,0,1,1};
		CloneGtVector gtv = new CloneGtVector(clone, 1, 0, 1);
		System.out.println(Arrays.toString(gtv.gtVector));
		System.out.println(gtv.gtVector[gtv.changeGtIndex] == gtv.newGt);
	}

}
