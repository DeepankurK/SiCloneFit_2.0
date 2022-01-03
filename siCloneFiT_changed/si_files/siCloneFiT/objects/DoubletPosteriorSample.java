/**
 * Oct 30, 2017
 */
package objects;

import java.util.ArrayList;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * @author hz22
 * Oct 30, 2017
 */
public class DoubletPosteriorSample {
	public ArrayList<Clone> listClone;
	public int[] cellCloneIDVector;
	public ArrayList<Integer> cell2ndCloneIDVector;
	public STITree<Clone> cloneTree;
	public STITree<Clone> cellCloneTree;
	public int[] cellDoubletFlagArr;
	public double logPosteriorProb;
	public double logLikelihood;
	public double sampleFp;
	public double sampleFn;
	public double sampleDoublet;
	
	public DoubletPosteriorSample(ArrayList<Clone> cloneL, int[] cellCloneID, ArrayList<Integer> cell2ndCloneVector, STITree<Clone> clonalTree, int[] doubletFlagArr, double posterior, double likelihood, double Fp, double Fn, double doublet){
		this.listClone = cloneL;
		this.cellCloneIDVector = cellCloneID;
		this.cell2ndCloneIDVector = cell2ndCloneVector;
		this.cloneTree = clonalTree;
		this.cellDoubletFlagArr = doubletFlagArr;
		this.logPosteriorProb = posterior;
		this.logLikelihood = likelihood;
		this.sampleFp = Fp;
		this.sampleFn = Fn;
		this.sampleDoublet = doublet;
	}
	
	public DoubletPosteriorSample(ArrayList<Clone> cloneL, ArrayList<Integer> cell2ndCloneVector, STITree<Clone> clonalTree, int[] doubletFlagArr, double likelihood, double Fp, double Fn, double doublet){
		this.listClone = cloneL;
		
		this.cell2ndCloneIDVector = cell2ndCloneVector;
		this.cloneTree = clonalTree;
		this.cellDoubletFlagArr = doubletFlagArr;
		
		this.logLikelihood = likelihood;
		this.sampleFp = Fp;
		this.sampleFn = Fn;
		this.sampleDoublet = doublet;
	}

	/**
	 * @param args
	 * Created On: Oct 30, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
