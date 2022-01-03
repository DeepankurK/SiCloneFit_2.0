/**
 * Oct 26, 2017
 */
package objects;

import java.util.ArrayList;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * @author hz22
 * Oct 26, 2017
 */
public class SingletPosteriorSample {
	public ArrayList<Clone> listClone;
	public int[] cellCloneIDVector;
	public STITree<Clone> cloneTree;
	public STITree<Clone> cellCloneTree;
	public double logPosteriorProb;
	public double logLikelihood;
	public double sampleFp;
	public ArrayList<Double> sampleFn;
	public double sampleDelProb;
	public double sampleLOHProb;

	public SingletPosteriorSample(ArrayList<Clone> cloneL, int[] cellCloneID, STITree<Clone> clonalTree, double posterior, double likelihood, double Fp, ArrayList<Double> Fn){
		this.listClone = cloneL;
		this.cellCloneIDVector = cellCloneID;
		this.cloneTree = clonalTree;
		this.logPosteriorProb = posterior;
		this.logLikelihood = likelihood;
		this.sampleFp = Fp;
		this.sampleFn = Fn;
	}
	
	public SingletPosteriorSample(ArrayList<Clone> cloneL, int[] cellCloneID, STITree<Clone> clonalTree, double posterior, double likelihood, double Fp, ArrayList<Double> Fn, double dProb, double LProb){
		this.listClone = cloneL;
		this.cellCloneIDVector = cellCloneID;
		this.cloneTree = clonalTree;
		this.logPosteriorProb = posterior;
		this.logLikelihood = likelihood;
		this.sampleFp = Fp;
		this.sampleFn = Fn;
		this.sampleDelProb = dProb;
		this.sampleLOHProb = LProb;
	}
	
	/**
	 * @param args
	 * Created On: Oct 26, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
