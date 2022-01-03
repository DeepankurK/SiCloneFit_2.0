/**
 * Aug 12, 2017
 */
package siCloneFiT.likelihood;

import java.util.ArrayList;

import org.apache.commons.math3.special.Gamma;

import SiFit.model.ComplexEvolutionModel;
import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Aug 12, 2017
 */
public class CompleteLikelihood {

	/**
	 * @param args
	 * Created On: Aug 12, 2017
	 */
	public TreeGenotypeLikelihood treeGenotypeLikelihoodObj;
	public SCSAmplificationErrorLikelihood ampErrLikelihoodObj;
	
	public CompleteLikelihood(Tree t, ComplexEvolutionModel m, double fp, ArrayList<Double> fn, int df){
		treeGenotypeLikelihoodObj = new TreeGenotypeLikelihood(t, m);
		ampErrLikelihoodObj = new SCSAmplificationErrorLikelihood(fp, fn, df);
	}
	
	/**
	 * Compute the full likelihood of the tree and error
	 * @param listCells
	 * @param listClone
	 * @param fp
	 * @param fn
	 * @param cloneGenotypeObsMatrix
	 * @param nMut
	 * @return
	 * Created On: Aug 12, 2017
	 */
	public double computeFullLikelihood(ArrayList<SingleCell> listCells, ArrayList<SingleCell> listCells_tot,ArrayList<Clone> listClone,
			double fp, ArrayList<Double> fn, ArrayList<GenotypeObservation> cloneGenotypeObsMatrix, int nMut){
		double t1 = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listCells, listCells_tot,listClone, fp, fn, nMut);
		double t2 = treeGenotypeLikelihoodObj.computeTreeCloneLogLikelihood(cloneGenotypeObsMatrix, nMut);
		return t1 + t2;
	}
	
	/**
	 * Compute the full likelihood of the tree and error
	 * @param listCells
	 * @param listClone
	 * @param tree
	 * @param m
	 * @param fp
	 * @param fn
	 * @param cloneGenotypeObsMatrix
	 * @param nMut
	 * @return
	 * Created On: Aug 15, 2017
	 */
	public double computeFullLikelihood(ArrayList<SingleCell> listCells, ArrayList<Clone> listClone, Tree tree,
			ComplexEvolutionModel m, double fp, ArrayList<Double> fn, ArrayList<GenotypeObservation> cloneGenotypeObsMatrix, int nMut){
		double t1 = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listCells, listClone, fp, fn, nMut);
		treeGenotypeLikelihoodObj = new TreeGenotypeLikelihood(tree, m);
		double t2 = treeGenotypeLikelihoodObj.computeTreeCloneLogLikelihood(cloneGenotypeObsMatrix, nMut);
		return t1 + t2;
	}
	
	/**
	 * Compute the full likelihood of the tree and error
	 * @param listCells
	 * @param listClone
	 * @param tree
	 * @param m
	 * @param fp
	 * @param fn
	 * @param cloneGenotypeObsMatrix
	 * @param nMut
	 * @param df
	 * @return
	 * Created On: Oct 11, 2017
	 */
	public double computeFullLikelihood(ArrayList<SingleCell> listCells,ArrayList<SingleCell> listCells_tot, ArrayList<Clone> listClone, Tree tree,
			ComplexEvolutionModel m, double fp, ArrayList<Double> fn, ArrayList<GenotypeObservation> cloneGenotypeObsMatrix, int nMut, int df){
		double t1 = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixFpFn(listCells,listCells_tot, listClone, fp, fn, nMut);
		treeGenotypeLikelihoodObj = new TreeGenotypeLikelihood(tree, m);
		double t2 = treeGenotypeLikelihoodObj.computeTreeCloneLogLikelihood(cloneGenotypeObsMatrix, tree, nMut, df);
		return t1 + t2;
	}
	
	/**
	 * Compute the full likelihood of the tree and error for the doublet model
	 * @param listCells
	 * @param cellDoubletFlagArr
	 * @param listClone
	 * @param c2CellCloneList
	 * @param cloneGenotypeObsMatrix
	 * @param tree
	 * @param m
	 * @param SCF
	 * @param fp
	 * @param fn
	 * @param nMut
	 * @param df
	 * @return
	 * Created On: Oct 30, 2017
	 */
	public double computeFullLikelihoodDoublet(ArrayList<SingleCell> listCells, int[] cellDoubletFlagArr,
											   ArrayList<Clone> listClone, ArrayList<Integer> c2CellCloneList, 
											   ArrayList<GenotypeObservation> cloneGenotypeObsMatrix, 
											   Tree tree, ComplexEvolutionModel m, SCFUtilityFunctions SCF,
											   double fp, double fn, int nMut, int df){
		double t1 = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listCells, listClone, c2CellCloneList, cellDoubletFlagArr, fp, fn, nMut, SCF, df);
		treeGenotypeLikelihoodObj = new TreeGenotypeLikelihood(tree, m);
		double t2 = treeGenotypeLikelihoodObj.computeTreeCloneLogLikelihood(cloneGenotypeObsMatrix, tree, nMut, df);
		return t1 + t2;
	}
	

	
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
