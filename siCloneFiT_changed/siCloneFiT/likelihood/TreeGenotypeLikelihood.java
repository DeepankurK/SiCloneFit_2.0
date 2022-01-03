/**
 * Aug 7, 2017
 */
package siCloneFiT.likelihood;

import static edu.rice.hj.Module1.forallChunked;

import java.util.ArrayList;

import SiFit.BasicUtilityFunctions;
import SiFit.model.ComplexEvolutionModel;
import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import jeigen.DenseMatrix;
import siCloneFiT.objects.Clone;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Aug 7, 2017
 */
public class TreeGenotypeLikelihood {
	
	public Tree tree;
	public ComplexEvolutionModel model;
	public static final String genotypes = "012";
	public double newParamSameTreeLikelihood;
	public double newTreeSameParamLikelihood;
	
	public TreeGenotypeLikelihood(Tree t, ComplexEvolutionModel m){
		this.tree = t;
		this.model = m;
	}
	
	/**
	 * Create the base case leaf matrix for the given genotype
	 * @param genotype
	 * @return
	 * Created On: Aug 7, 2017
	 */
    private static DenseMatrix createLeafMatrix(Integer genotype){
    	/* Handle missing data. Missing data will be 
    	 * encoded by '-' or '3' and for that
    	 * each genotype is 1.
    	 */
    	if (genotype == 3){
            double [][] t = {{1},{1},{1}};
            DenseMatrix D = new DenseMatrix(t);
            return D;
    	}
    	double [][] temp = {
                {0},
                {0},
                {0}            
            };
    	if (genotypes.indexOf(Integer.toString(genotype)) == -1)
    		throw new RuntimeException(genotype +" is not a valid genotype");

    	// Evolution model check
    	temp[genotype][0] = 1;
    	return new DenseMatrix(temp);    		    	   	
    }
    
    private static DenseMatrix createLeafMatrixBinary(Integer genotype){
    	/* Handle missing data. Missing data will be 
    	 * encoded by '-' or '3' and for that
    	 * each genotype is 1.
    	 */
    	if (genotype == 3){
            double [][] t = {{1},{1}};
            DenseMatrix D = new DenseMatrix(t);
            return D;
    	}
    	double [][] temp = {
                {0},
                {0}                          
            };
    	if (genotypes.indexOf(Integer.toString(genotype)) == -1)
    		throw new RuntimeException(genotype +" is not a valid genotype");

    	// Evolution model check
    	temp[genotype][0] = 1;
    	return new DenseMatrix(temp);    		    	   	
    }
	
    /**
     * This returns the likelihood of a node having each genotype
     */
    private DenseMatrix probabilityOfNode(TNode node, GenotypeObservation obs_genotypes){
    	if (node.isLeaf()){    	    		
    		return createLeafMatrix(obs_genotypes.cellGenotypeMap.get(node.getName()));
    	}
    	else{
    		double[][] probs = new double[][]{{1.0}, {1.0}, {1.0}};
            for (TNode child : node.getChildren()) {
                double dist = child.getParentDistance();
                if (Double.isNaN(dist) || Double.isInfinite(dist))
                    throw new IllegalArgumentException("Node : " + child.getName() + " has non-finite dist: "+ dist);
                double[] edgeP = new double[3];
                DenseMatrix childProbs = probabilityOfNode(child, obs_genotypes);
             
                double[][] TransMat = new double[3][3];
                // Now just testing the new complex evolution model
                
                TransMat = model.getTransitionMatrixExp(dist);
                for (int i = 0; i < 3; i++){
                	edgeP[0] += TransMat[0][i] * childProbs.get(i, 0);
                	edgeP[1] += TransMat[1][i] * childProbs.get(i, 0);
                	edgeP[2] += TransMat[2][i] * childProbs.get(i, 0);
                	
                }
                probs[0][0] *= edgeP[0]; probs[1][0] *= edgeP[1]; probs[2][0] *= edgeP[2];
            }
            DenseMatrix probMat = new DenseMatrix(probs);
            return probMat;
    	}
    }
    
    private DenseMatrix probabilityOfNodeBinary(TNode node, GenotypeObservation obs_genotypes){
    	if (node.isLeaf()){    	    		
    		return createLeafMatrixBinary(obs_genotypes.cellGenotypeMap.get(node.getName()));
    	}
    	else{
    		double[][] probs = new double[][]{{1.0}, {1.0}};
            for (TNode child : node.getChildren()) {
                double dist = child.getParentDistance();
                if (Double.isNaN(dist) || Double.isInfinite(dist))
                    throw new IllegalArgumentException("Node : " + child.getName() + " has non-finite dist: "+dist);
                double[] edgeP = new double[2];
                DenseMatrix childProbs = probabilityOfNodeBinary(child, obs_genotypes);
             
                double[][] TransMat = new double[2][2];
                // Now just testing the new complex evolution model
                
                TransMat = model.getTransitionMatrixBinary(dist);
                for (int i = 0; i < 2; i++){
                	edgeP[0] += TransMat[0][i] * childProbs.get(i, 0);
                	edgeP[1] += TransMat[1][i] * childProbs.get(i, 0);
                	
                	
                }
                probs[0][0] *= edgeP[0]; probs[1][0] *= edgeP[1]; 
            }
            DenseMatrix probMat = new DenseMatrix(probs);
            return probMat;
    	}
    }
    
    /**
     * Function that returns the likelihood of genotype observation at given locus
     * @param arg0
     * @return
     * Created On: Aug 7, 2017
     */
	public double getGenotypeProbability(GenotypeObservation arg0) {		
		return probabilityOfNode(tree.getRoot(), arg0).mul(model.getEquilibriumVector()).sum().s(); // returns densematrix.get(0,0);
	}
	
	public double getGenotypeProbability(GenotypeObservation arg0, Tree nuTree) {		
		return probabilityOfNode(nuTree.getRoot(), arg0).mul(model.getEquilibriumVector()).sum().s(); // returns densematrix.get(0,0);
	}
	
	/**
	 * Function that returns the likelihood of genotype observation at given locus for binary data
	 * @param arg0
	 * @param nuTree
	 * @return
	 * Created On: Aug 16, 2017
	 */
	public double getGenotypeProbabilityBinary(GenotypeObservation arg0, Tree nuTree) {		
		return probabilityOfNodeBinary(nuTree.getRoot(), arg0).mul(model.getEquilibriumVectorBinary()).sum().s(); // returns densematrix.get(0,0);
	}
	
	public double getGenotypeProbabilityBinary(GenotypeObservation arg0) {		
		return probabilityOfNodeBinary(tree.getRoot(), arg0).mul(model.getEquilibriumVectorBinary()).sum().s(); // returns densematrix.get(0,0);
	}
	
	/**
	 * Compute the Log likelihood of a genotype observation given tree and model
	 * Also applicable when genotype observation corresponds to a single genotype of a single leaf
	 * This is the prior probability of that single genotype of that leaf leaf given other leaves genotypes 
	 * @param arg0
	 * @param df
	 * @return
	 * Created On: Oct 13, 2017
	 */
	public double getLogLikelihoodGtObs(GenotypeObservation arg0, int df){
		if (df == 0)
			return Math.log(getGenotypeProbabilityBinary(arg0));
		else
			return Math.log(getGenotypeProbability(arg0));
	}
	
	/**
	 * Compute the log-likelihood ratio of two genotype vectors (differs only in one position) of the same clone
	 * @param newGtObs
	 * @param currGtObs
	 * @param df
	 * @return
	 * Created On: Oct 13, 2017
	 */
	public double getNewGtObsLogLikelihoodRatio(GenotypeObservation newGtObs, GenotypeObservation currGtObs, int df){
		double newGtObsLogLikelihood = getLogLikelihoodGtObs(newGtObs, df);
		double currGtObsLogLikelihood = getLogLikelihoodGtObs(currGtObs, df);
		return newGtObsLogLikelihood - currGtObsLogLikelihood;
	}
	
	/**
	 * Return the likelihood of a tree
	 * @param cloneGenotypeObsMatrix
	 * @param nMut
	 * @return
	 * Created On: Aug 12, 2017
	 */
	public double computeTreeCloneLogLikelihood(ArrayList<GenotypeObservation> cloneGenotypeObsMatrix,
												int nMut){
		double sum = 0;
		double[] val = new double[nMut];
		for (int i = 0; i < nMut; i++){
			val[i] = getGenotypeProbability(cloneGenotypeObsMatrix.get(i));
		}
		for (double d: val){
//			System.out.println(d);
			double log_d;
			if (d == 0)
				log_d = -10000;
			else
				log_d = Math.log(d);
			sum += log_d;
		}
		return sum;
//		forallChunked(0, nMut -1, (i) -> {
//			val[i] = getGenotypeProbability(cloneGenotypeObsMatrix.get(i));
//		});
	}
	
	/**
	 * Compute the logLikelihood of a new parameter value for the model of evolution
	 * @param cloneGenotypeObsMatrix
	 * @param tree
	 * @param newParamVal
	 * @param nMut
	 * @param df
	 * @param modelParamFlag
	 * @return
	 * Created On: Sep 14, 2017
	 */
	public double computeTreeCloneNewParamLogLikelihood(ArrayList<GenotypeObservation> cloneGenotypeObsMatrix,
														Tree tree, double newParamVal, int nMut, int df, int modelParamFlag){
		double currParamVal = 0;
		// delProb is being used
		if (modelParamFlag == 0){
			currParamVal = model.delProb;
			model.delProb = newParamVal;
		}
		// LOHProb is being used
		else if (modelParamFlag == 1){
			currParamVal = model.LOHProb;
			model.LOHProb = newParamVal;
		}
		double logLikelihood = computeTreeCloneLogLikelihood(cloneGenotypeObsMatrix, tree, nMut, df);
		if (modelParamFlag == 0){
			model.delProb = currParamVal;
		}
		else{
			model.LOHProb = currParamVal;
		}
		return logLikelihood;
		
	}
	public double computeTreeCloneLogLikelihood(ArrayList<GenotypeObservation> cloneGenotypeObsMatrix,
			                                    Tree nuTree, int nMut, int df){
		double sum = 0;
		double[] val = new double[nMut];
		if (df == 1){
			for (int i = 0; i < nMut; i++){
				val[i] = getGenotypeProbability(cloneGenotypeObsMatrix.get(i), nuTree);
			}
		}
		else{
			for (int i = 0; i < nMut; i++){
				val[i] = getGenotypeProbabilityBinary(cloneGenotypeObsMatrix.get(i), nuTree);
			}
		}
		for (double d: val){
			//System.out.println(d);
			double log_d;
			if (d == 0)
				log_d = -10000;
			else
				log_d = Math.log(d);
			sum += log_d;
		}
		return sum;
		//forallChunked(0, nMut -1, (i) -> {
		//val[i] = getGenotypeProbability(cloneGenotypeObsMatrix.get(i));
		//});
	}
	
	/**
	 * Compute the Likelihood ratio for the nu tree obtained by placing the new clone as a leaf
	 * @param nuCloneGenotypeObsMatrix
	 * @param oldCloneGenotypeObsMatrix
	 * @param nuTree
	 * @param oldTree
	 * @param nMut
	 * @return
	 * Created On: Aug 14, 2017
	 */
	public double computeTreeCloneLikelihoodRatio(ArrayList<GenotypeObservation> nuCloneGenotypeObsMatrix, 
												  ArrayList<GenotypeObservation> oldCloneGenotypeObsMatrix,
												  Tree nuTree, Tree oldTree, int nMut, int df){
		double nuCloneTreeLogLikelihood = this.computeTreeCloneLogLikelihood(nuCloneGenotypeObsMatrix, nuTree, nMut, df);
//		this.newTreeSameParamLikelihood = nuCloneTreeLogLikelihood;
		double oldCloneTreeLogLikelihood = this.computeTreeCloneLogLikelihood(oldCloneGenotypeObsMatrix, oldTree, nMut, df);
//		System.out.printf("nuClone = %f\toldClone = %f\n", nuCloneTreeLogLikelihood, oldCloneTreeLogLikelihood);
		return Math.exp(nuCloneTreeLogLikelihood - oldCloneTreeLogLikelihood);
	}
	
	/**
	 * Compute Log-likelihood ratio when nu tree on same number of leaves is proposed. 
	 * Used for sampling tree from conditional posterior
	 * log( P(G | T_nu, M) / P(G | T_old, M))
	 * @param cloneGenotypeObsMatrix, clone Genotype observation matrix
	 * @param nuTree, new tree being proposed
	 * @param oldTree, old tree
	 * @param nMut
	 * @param df
	 * @return
	 * Created On: Sep 18, 2017
	 */
	public double computeTreeCloneLogLikelihoodRatio(ArrayList<GenotypeObservation> cloneGenotypeObsMatrix,
													 Tree nuTree, Tree oldTree, int nMut, int df){
		double nuCloneTreeLogLikelihood = this.computeTreeCloneLogLikelihood(cloneGenotypeObsMatrix, nuTree, nMut, df);
		this.newTreeSameParamLikelihood = nuCloneTreeLogLikelihood;
		double oldCloneTreeLogLikelihood = this.computeTreeCloneLogLikelihood(cloneGenotypeObsMatrix, oldTree, nMut, df);
		return nuCloneTreeLogLikelihood - oldCloneTreeLogLikelihood;
	}
	
	/**
	 * Compute Log-likelihood ratio when nu tree on same number of leaves is proposed. 
	 * Used for sampling tree from conditional posterior
	 * log( P(G | T_nu, M) / P(G | T_old, M))
	 * @param cloneGenotypeObsMatrix
	 * @param currLikelihood
	 * @param nuTree
	 * @param nMut
	 * @param df
	 * @return
	 * Created On: Sep 18, 2017
	 */
	public double computeTreeCloneLogLikelihoodRatio(ArrayList<GenotypeObservation> cloneGenotypeObsMatrix,
			 										 double currLikelihood, Tree nuTree, int nMut, int df){
		double nuCloneTreeLogLikelihood = this.computeTreeCloneLogLikelihood(cloneGenotypeObsMatrix, nuTree, nMut, df);
		this.newTreeSameParamLikelihood = nuCloneTreeLogLikelihood;
		return nuCloneTreeLogLikelihood - currLikelihood;
	}
	
	/**
	 * Compute log-likelihood ratio for proposing a new value of model parameter
	 * @param CloneGenotypeObsMatrix
	 * @param currTree
	 * @param newParamVal
	 * @param currLikelihood
	 * @param nMut
	 * @param df
	 * @param paramFlag
	 * @return
	 * Created On: Sep 18, 2017
	 */
	public double computeTreeCloneNewParamLogLikelihoodRatio(ArrayList<GenotypeObservation> CloneGenotypeObsMatrix,
														  Tree currTree, double newParamVal, double currLikelihood, int nMut, int df, int paramFlag){
		double newParamTreeCloneLogLikelihood = this.computeTreeCloneNewParamLogLikelihood(CloneGenotypeObsMatrix, currTree, newParamVal, nMut, df, paramFlag);
		this.newParamSameTreeLikelihood = newParamTreeCloneLogLikelihood;
		return newParamTreeCloneLogLikelihood - currLikelihood;
	}
	
	/**
	 * Computes the probability of a leaf having a genotype (single locus) given the tree and genotypes at the other leaves
	 * The tree here is the same as used to instantiate the class, the tree should have the leafset consisting of
	 * the current clusters as well as the new cluster
	 * @param OtherLeavesGenotypes, Genotypes at the other leaves
	 * @param leafName, name of the leaf (String), whose genotype is sampled 
	 * @param genotype, the genotype (0/1/2)
	 * @return
	 * Created On: Aug 7, 2017
	 */
	public double getSingleLeafGenotypeProbGivenTreeOtherGenotypes(GenotypeObservation OtherLeavesGenotypes,
																   String leafName, Integer genotype){
		GenotypeObservation allLeavesGenotypes = OtherLeavesGenotypes;
		allLeavesGenotypes.cellGenotypeMap.put(leafName, genotype);
		double val = getGenotypeProbability(allLeavesGenotypes);
//		System.out.printf("clone = %s, genotype = %d, likelihood = %f\n", leafName, genotype, val);
		return val;
	}
	
	/**
	 * Computes the probability of a leaf having a genotype (single locus) given the tree and genotypes at the other leaves
	 * The tree here is the same as used to instantiate the class, the tree should have the leafset consisting of
	 * the current clusters as well as the new cluster
	 * @param OtherLeavesGenotypes
	 * @param leafName
	 * @param genotype
	 * @param df
	 * @return
	 * Created On: Oct 2, 2017
	 */
	public double getSingleLeafGenotypeProbGivenTreeOtherGenotypes(GenotypeObservation OtherLeavesGenotypes,
			   String leafName, Integer genotype, int df){
//		System.out.println(leafName);
		GenotypeObservation allLeavesGenotypes = OtherLeavesGenotypes;
//		System.out.println("before " + allLeavesGenotypes.cellGenotypeMap);
		allLeavesGenotypes.cellGenotypeMap.put(leafName, genotype);
//		System.out.println("after " + allLeavesGenotypes.cellGenotypeMap);
		double val;
		if (df == 1)
			val = getGenotypeProbability(allLeavesGenotypes);
		else
			val = getGenotypeProbabilityBinary(allLeavesGenotypes);
		//System.out.printf("clone = %s, genotype = %d, likelihood = %f\n", leafName, genotype, val);
		return val;
}
	
	
	

	/**
	 * @param args
	 * Created On: Aug 7, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String t = "((C0:0.1,C1:0.05):0.02, (C2:0.04, C3:0.06):0.07);";
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
		STITree<Double> tree = SCF.getTree(t);
		ComplexEvolutionModel m = new ComplexEvolutionModel(0.01, 0.01);
		TreeGenotypeLikelihood TGL = new TreeGenotypeLikelihood(tree, m);
		
		ArrayList<String> leaves = new ArrayList<>();
		for (int i = 0; i < 3; i++){
			leaves.add("C"+Integer.toString(i));
		}
		Integer[] gt = new Integer[]{0,1,1};
		GenotypeObservation G = new GenotypeObservation(leaves, gt);
		for (int i = 0; i < 3; i++){
			System.out.println(TGL.getSingleLeafGenotypeProbGivenTreeOtherGenotypes(G, "C3", i));
		}
	}

}
