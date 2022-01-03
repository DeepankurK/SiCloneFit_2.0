/**
 * Oct 14, 2017
 */
package objects;

import java.util.ArrayList;

/**
 * @author hz22
 * Oct 14, 2017
 */
public class ChangeCloneGtMatrix {
	
	public ArrayList<Integer[]> newCloneGtMatrix;
	public int currGt;
	public int newGt;
	public int changeGtIndex;
	public String changeCloneName;
	public int changeCloneID;
	
	public ChangeCloneGtMatrix(ArrayList<Integer[]> newMat, int cg, int ng, int index, String cloneName, int cloneID){
		this.newCloneGtMatrix = newMat;
		this.currGt = cg;
		this.newGt = ng;
		this.changeGtIndex = index;
		this.changeCloneName = cloneName;
		this.changeCloneID = cloneID;
	}

	/**
	 * @param args
	 * Created On: Oct 14, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
