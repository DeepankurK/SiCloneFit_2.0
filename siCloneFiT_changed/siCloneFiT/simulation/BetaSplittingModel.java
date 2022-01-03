/**
 * Oct 24, 2017
 */
package siCloneFiT.simulation;

import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * @author hz22
 * Oct 24, 2017
 */
public class BetaSplittingModel {
	
    protected static final int FIRST = 0;
    protected static final int SECOND = 1;
    protected static final double RIGHT = 1;
    protected static final double LEFT = 0;
    protected static final int BRANCH_LENGTH = 2;
    protected static final int CHILD_SIDE = 3;
    public ExponentialDistribution branchLengthDist;

    double Alpha, Beta, Lambda;
    int Taxa;
    double[] Speciation_rates;
    double[] Speciation_events;
    double[] Branches;
    
    STITree<double []> Beta_Tree;
    
    public BetaSplittingModel(double alpha_param, double beta_param, double exponential_mean, int taxa_num) {
        Alpha = alpha_param;
        Beta = beta_param;
        Taxa = taxa_num;
        Lambda = exponential_mean;
        Speciation_rates = Gen_Beta_Variables(Taxa - 1);
        Speciation_events = Gen_Uniform_Variables(Taxa - 1);
        Branches = Gen_Exponential_Variables(Taxa - 1);
        branchLengthDist = new ExponentialDistribution(exponential_mean);
    }
    
    private double[] Gen_Beta_Variables(int taxa_number) {
//        double random;
        double x,b;
        double[] output = new double[taxa_number];
        BetaDistribution beta = new BetaDistribution(Alpha+1.0, Beta+1.0);
        for (int i = 0; i < taxa_number; i++) {
            x = Math.random();
            b = beta.inverseCumulativeProbability(x);
//            b = beta.sample();
            output[i]=b;
        }
//        System.out.println(Arrays.toString(output));
        return output;
    }

    private double[] Gen_Uniform_Variables(int taxa_number) {
        double x;
        double[] output = new double[taxa_number];
        for (int i = 0; i < taxa_number; i++) {
            x = Math.random();
            output[i] = x;
        }
        return output;
    }
    
    private double[] Gen_Exponential_Variables(int taxa_number) {
        double random;
        double[] output = new double[taxa_number];
        ExponentialDistribution exponential = new ExponentialDistribution(Lambda);
        for (int i = 0; i < taxa_number; i++) {
            random = exponential.sample();
            output[i] = random;
        }
        /**
         * Normalization of output elements
         */
        double sum = 0;
        for (int j = 0; j < output.length; j++) {
            sum = output[j] + sum;
        }
        for (int k = 0; k < output.length; k++) {
            output[k] = (double) (output[k] / sum);
        }
        return output;
    }

    private boolean is_in(double a, double b, double c) {
        if (a >= b && a <= c) {
            return true;
        } else {
            return false;
        }
    }
    
    public String Construct_Tree() {
        int Speciation_number = 0;
        /**
         * Initialization of the tree
         */
        Beta_Tree = new STITree<>("root", true);
        Beta_Tree.getRoot().setData(new double[]{0, 1});
        STINode<double[]> tempNode;
        tempNode = Beta_Tree.getRoot().createChildWithUniqueName();

        /**
         * Data in this case: [a,b,c,d] in which [a,b] is the speciation interval, c is the branch length, and d is the left or right label.
         */
        tempNode.setParentDistance(branchLengthDist.sample());
        tempNode.setData(new double[]{0, Speciation_rates[Speciation_number], Branches[Speciation_number], LEFT});
        tempNode = Beta_Tree.getRoot().createChildWithUniqueName();
        tempNode.setData(new double[]{Speciation_rates[Speciation_number], 1, Branches[Speciation_number], RIGHT});
        tempNode.setParentDistance(branchLengthDist.sample());
        Speciation_number++;
        double middle;
        STINode<double[]> iterativeNode;
        /**
         * While loop for constructing the rest of the tree
         */
        while (Speciation_number < (Taxa - 1)) {
            Iterator<STINode<double[]>> i = Beta_Tree.getNodes().iterator();
            while (i.hasNext()) {
                iterativeNode = i.next();
                if (iterativeNode.isLeaf() && is_in(Speciation_events[Speciation_number], iterativeNode.getData()[FIRST], iterativeNode.getData()[SECOND])) {
                    /**
                     * Two new children are born here; the left child is labeled by 0 and 1 is for the right one.
                     */
                    middle = (double)Speciation_rates[Speciation_number] * ((double)(iterativeNode.getData()[SECOND] - iterativeNode.getData()[FIRST]))
                            + (double)iterativeNode.getData()[FIRST];
                    tempNode = iterativeNode.createChildWithUniqueName();
                    tempNode.setParentDistance(branchLengthDist.sample());
                    tempNode.setData(new double[]{iterativeNode.getData()[FIRST], middle, Branches[Speciation_number], LEFT});
                    tempNode = iterativeNode.createChildWithUniqueName();
                    tempNode.setParentDistance(branchLengthDist.sample());
                    tempNode.setData(new double[]{middle, iterativeNode.getData()[SECOND], Branches[Speciation_number], RIGHT});
                    Speciation_number++;
                    break;
                }
            }
        }
        return Beta_Tree.toNewick();
    }

    /**
     * Finds all the descendants of a given node
     * @param node
     * @return
     */

    private int internal_descendants_Count(STINode<double []> node){
        int internal_descendants=0;
        if (node.isLeaf()){
            return 0;
        }
        else{
            for (STINode<double []> pass_node : node.getChildren()){
                internal_descendants =  internal_descendants + internal_descendants_Count(pass_node);
            }
            return internal_descendants+1;
        }
    }

    public double Probability(STINode<double []> Root){
        int left_internal_nodes=0;
        int right_internal_nodes=0;
        double logP=0;
        if (!Root.isLeaf()){
            for (STINode<double []> node : Root.getChildren()){
                logP = logP + Probability(node);
                if (node.getData()[CHILD_SIDE]==LEFT){
                    left_internal_nodes = internal_descendants_Count(node);
                }
                if (node.getData()[CHILD_SIDE]==RIGHT){
                    right_internal_nodes = internal_descendants_Count(node);
                }
            }
            logP = logP + (org.apache.commons.math3.special.Beta.logBeta(Alpha+1.0+(double)(left_internal_nodes), Beta+1.0+(double)(right_internal_nodes))-
            org.apache.commons.math3.special.Beta.logBeta(Alpha+1.0,Beta+1.0));
            return logP;
        }
        return logP;
    }
    
    public static void test(int Ts) throws IOException, ParseException {
        BetaSplittingModel test = new BetaSplittingModel(1000000000, 1000000000, 1,Ts);
        test.Construct_Tree();

    }

	/**
	 * @param args
	 * Created On: Oct 24, 2017
	 * @throws ParseException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, ParseException {
		// TODO Auto-generated method stub
		test(10);
	}

}
