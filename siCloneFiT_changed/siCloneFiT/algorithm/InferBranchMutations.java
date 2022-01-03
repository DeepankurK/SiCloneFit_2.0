/**
 * Mar 9, 2018
 */
package siCloneFiT.algorithm;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;

import org.apache.commons.lang3.StringUtils;

import SiFit.BasicUtilityFunctions;
//import SiFit.algorithm.InferAncestralStates.AncestralInfo;
import SiFit.io.VariantMatrixReader;
import SiFit.model.ComplexEvolutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * @author hz22
 * Mar 9, 2018
 */
public class InferBranchMutations {
	
	protected static class AncestralInfo{
		
		// fields
		/**
		 * The sequence at this node.  Generated for internal nodes, set in leaves.
		 */
		public Integer[] sequence;
		
		/**
		 * The node in the input tree that this node corresponds to.
		 */
		public TMutableNode node;
		
		/**
		 * The list of best genotypes for each mutation for the node.
		 * These will be filled during bottom up phase of the algorithm
		 * and used during the top down phase of the algorithm.
		 */
		public ArrayList<ArrayList<Integer>> setGTNode = new ArrayList<>();
		
		/**
		 * The list of likelihoods for each genotype for each mutation for the node.
		 * These will be filled during bottom up phase of the algorithm
		 * and used during the top down phase of the algorithm.
		 */
		public ArrayList<ArrayList<Double>> likelihoodGTNode = new ArrayList<>();
		
		public AncestralInfo(TMutableNode node, int len) {
			this.node = node;
			for (int i = 0; i < len; i++){
				setGTNode.add(new ArrayList<Integer>());
				likelihoodGTNode.add(new ArrayList<Double>());
			}
			
		}
	}
	
	public Tree tree;
	public ComplexEvolutionModel modelC;
	public static int dataFlag = 0;
	public static String genotypes = "01";
	public ArrayList<String> geneNamesList;			// List of genes on which mutations occur
    public STITree<AncestralInfo> nuPtree;		// Copy of the original tree but with branches that correspond to number of mutations
    
    /**
     * Constructor
     * @param delProb
     * @param LOHProb
     * @param recurProb
     * @param df
     */
    public InferBranchMutations(double delProb, double LOHProb, double recurProb,
    							int df){
    	this.modelC = new ComplexEvolutionModel(delProb, LOHProb, recurProb);
        dataFlag = df;
    }
    
    /**
     * 
     * @param tree
     * @param cellGTSeqMap
     * @param seqLen
     * @return
     * Created On: Mar 9, 2018
     */
    protected int computeAncestralStates(MutableTree tree, 
										 Hashtable<String, Integer[]> cellGTSeqMap,
										 int seqLen){
    	STITree<AncestralInfo> ptree = new STITree<>();
    	
    	// copy the tree
    	copyNode(ptree.getRoot(), tree.getRoot(), cellGTSeqMap, seqLen);
    	
    	for (STINode<AncestralInfo> n : ptree.getNodes()){
    		n.setName(n.getData().node.getName());
    	}

    	// Perform bottom up and top down phases
    	for(int pos = 0; pos < seqLen; pos++) {
    		performBottomUp(pos, ptree.getRoot());
    		performTopDown(pos, ptree.getRoot());
    	}
    	
    	this.nuPtree = ptree;
    	
    	return computePScore(ptree.getRoot());
    }
    
    /**
     * Top Down phase
     * @param pos
     * @param pnode
     * Created On: Mar 9, 2018
     */
    protected void performBottomUp(int pos, STINode<AncestralInfo> pnode){
    	AncestralInfo pi = pnode.getData();
    	// Leaf
    	if (pnode.isLeaf()){
    		// Get the transition matrix with branch length to the parent of leaf
    		double[][] transMat = modelC.getTransitionMatrixBinary(pi.node.getParentDistance());
    		// Invalid Genotype
    		if (genotypes.indexOf(Integer.toString(pi.sequence[pos])) == -1){
    			throw new RuntimeException(pi.sequence[pos] +" is not a valid genotype");
    		}
    		else{
    			for (int i = 0; i < genotypes.length(); i++){
    				pi.setGTNode.get(pos).add(pi.sequence[pos]);
    				double ll_i = Math.log(transMat[i][pi.sequence[pos]]);
    				pi.likelihoodGTNode.get(pos).add(ll_i);
    			}
    		}
    	}
    	else if (pnode.isRoot()){
    		pi.setGTNode.get(pos).add(0);
    		for (STINode<AncestralInfo> child : pnode.getChildren()) {
    			performBottomUp(pos, child);
    		}
    	}
    	else{
    		ArrayList<STINode<AncestralInfo>> p_children = new ArrayList<>();
    		for(STINode<AncestralInfo> child : pnode.getChildren()) {				
				performBottomUp(pos, child);
				p_children.add(child);
			}
    		double[][] JCTransMat = modelC.getTransitionMatrixBinary(pi.node.getParentDistance());
    		for (int i = 0; i < genotypes.length(); i++){
    			double bestL = Double.NEGATIVE_INFINITY;
				int bestGT = 0;
				double ll_i;
				for (int j = 0; j < genotypes.length(); j++){
					ll_i = Math.log(JCTransMat[i][j]) + p_children.get(0).getData().likelihoodGTNode.get(pos).get(j) + p_children.get(1).getData().likelihoodGTNode.get(pos).get(j);
					if (ll_i > bestL){
						bestL = ll_i;
						bestGT = j;
					}
				}
				pi.likelihoodGTNode.get(pos).add(bestL);
				pi.setGTNode.get(pos).add(bestGT);	
    		}
    	}
    }
    
    /**
     * 
     * @param pos
     * @param pnode
     * Created On: Mar 9, 2018
     */
    protected void performTopDown(int pos, STINode<AncestralInfo> pnode){
    	AncestralInfo pi = pnode.getData();
    	if (pnode.isRoot()){
    		pi.sequence[pos] = 0;
    	}
    	else{
    		int index = pnode.getParent().getData().sequence[pos];
    		pi.sequence[pos] = pi.setGTNode.get(pos).get(index);	
    	}
//    	else if (pnode.isLeaf() == false){
//    		int index = pnode.getParent().getData().sequence[pos];
//    		pi.sequence[pos] = pi.setGTNode.get(pos).get(index);
//    	}
//    	else{
//    		pi.sequence[pos] = pnode.getParent().getData().sequence[pos];
//    	}
    	for (STINode<AncestralInfo> child : pnode.getChildren()) {
			performTopDown(pos, child);
		}
    }
    
    /**
     * Compute the number of mutations in a tree
     * @param pnode
     * @return
     * Created On: Mar 9, 2018
     */
	protected int computePScore(STINode<AncestralInfo> pnode) {
		
		int pscore = 0;
		if(pnode.isRoot()) {
			pnode.setParentDistance(0);
		} 
		else {
			Integer[] ps = pnode.getParent().getData().sequence;
			Integer[] seq = pnode.getData().sequence;
			
			int num_diffs = 0;
			
			for(int i = 0; i < ps.length; i++) {
				if(ps[i] != seq[i]) {
					num_diffs++;
				}
			}

			pnode.setParentDistance(num_diffs);
		}
		
		// set the mirror node's distance
//		pnode.getData().node.setParentDistance(pnode.getParentDistance());
		
		// compute all children's pscores
		for(STINode<AncestralInfo> child : pnode.getChildren()) {
			pscore += computePScore(child) + child.getParentDistance();
		}
		
		return pscore;		
	}
    
    /**
     * Copy the input tree
     * @param pnode
     * @param n
     * @param cellGTSeqMap
     * @param seqLen
     * Created On: Mar 9, 2018
     */
    protected void copyNode(STINode<AncestralInfo> pnode, TMutableNode n, Hashtable<String, Integer[]> cellGTSeqMap, int seqLen){
		AncestralInfo pi = new AncestralInfo(n, seqLen);
		pnode.setData(pi);
		if(n.isLeaf()) {
			pi.sequence = cellGTSeqMap.get(n.getName());
			if(pi.sequence == null) {
				throw new RuntimeException("No sequence provided for leaf node " + n.getName());
			}
		}
		else{
			pi.sequence = new Integer[seqLen];
			for(TMutableNode c : n.getChildren()) {
				STINode<AncestralInfo> pchild = pnode.createChild();
				copyNode(pchild, c, cellGTSeqMap, seqLen);
			}
		}
	}
    
    /**
     * 
     * @param treeFile
     * @param varMatFile
     * @param singleCellFile
     * @param geneNamesFile
     * @throws IOException
     * Created On: Mar 9, 2018
     */
    public void annotateMutations(String treeFile, String varMatFile, String singleCellFile, String geneNamesFile, String expectedMatrixFile) throws IOException{
    	// Classes to use
    	BasicUtilityFunctions BUF = new BasicUtilityFunctions();

    	// Get the input tree
    	String intreeNewick = BUF.readNewickString(treeFile);
    	STITree<Double> intree = BUF.getTree(intreeNewick);

    	// Name the internal nodes of tree
    	intree.getRoot().setName("R");
    	int i = 1;
    	for (STINode<Double> node: intree.getNodes()){			
    		if (node.getName() == ""){
    			node.setName("in"+ Integer.toString(i));
    			i++;
    		}
    	}
    	System.out.println(intree.toNewick());
    	// Create the single cell genomes to use as input to parsimony score calculator
    	int nCell = intree.getLeafCount();
    	VariantMatrixReader vr = new VariantMatrixReader(varMatFile);
    	vr.populateVarGTMatrix(varMatFile, nCell);
    	ArrayList<Integer[]> obsGenotypeMatrix = vr.obsVarGTMatrix;

    	int nMut = obsGenotypeMatrix.size();
    	ArrayList<Integer[]> seqGT = BUF.getCellGenomes(obsGenotypeMatrix, nMut);

    	// Get the single cell names
    	String singleCellNameString = BUF.readNewickString(singleCellFile);
    	String[] cellNameArr = singleCellNameString.split(" ");

    	Hashtable<String, Integer[]> cellGTSeqMap = new Hashtable<>();
    	for (int j = 0; j < cellNameArr.length; j++){
    		cellGTSeqMap.put(cellNameArr[j], seqGT.get(j));
    	}

    	// Get the Gene Names
    	ArrayList<String> geneList = BUF.readGeneNames(geneNamesFile);
    	this.geneNamesList = geneList;
    	this.tree = intree;
    	int score = this.computeAncestralStates(intree, cellGTSeqMap, nMut);
    	System.out.println(score);
    	
    	for (STINode<AncestralInfo> node : this.nuPtree.getNodes()){
    		if (node.isRoot() == false){
    			if (node.isLeaf())
					continue;
				ArrayList<String> mutList = this.getMutationsP2C(node.getParent(), node, geneList, nMut);
				if (mutList.size()>0){
					this.printMutationsP2C(node.getParent().getName(), node.getName(), geneList, nMut);
				}
    		}
    	}
    	
    	HashMap<String, ArrayList<String>> cellGenomeMap = new HashMap<>();
		for (STINode<AncestralInfo> node : this.nuPtree.getNodes()){
			if (node.isLeaf()){
				ArrayList<String> cellGenome = new ArrayList<>();
				cellGenome.add(node.getData().node.getName());
				for (Integer g : node.getParent().getData().sequence)
					cellGenome.add(Integer.toString(g));
				cellGenomeMap.put(node.getData().node.getName(), cellGenome);
			}
		}
		
		try{
			PrintWriter writer = new PrintWriter(expectedMatrixFile, "UTF-8");
			for (String s : cellGenomeMap.keySet()){
				String cell = StringUtils.join(cellGenomeMap.get(s), ' ');
				writer.println(cell);
			}
			
			writer.close();
		} catch(IOException e){
			e.printStackTrace();
		};
    }
    
    /**
     * 
     * @param par
     * @param child
     * @param geneNames
     * @param nMut
     * @return
     * Created On: Mar 9, 2018
     */
	public ArrayList<String> getMutationsP2C(STINode<AncestralInfo> par, 
			 STINode<AncestralInfo> child,
			 ArrayList<String> geneNames,
			 int nMut){
		ArrayList<String> mutatedGenes = new ArrayList<>();
		for (int i = 0; i < nMut; i++){
			if (par.getData().sequence[i] != child.getData().sequence[i]){
				mutatedGenes.add(geneNames.get(i));
			}
		}
		return mutatedGenes;
	}

	/**
	 * 
	 * @param parname
	 * @param childname
	 * @param geneNames
	 * @param nMut
	 * @param mutPosFile
	 * @throws IOException
	 * Created On: Mar 9, 2018
	 */
	public void printMutationsP2C(String parname, String childname, ArrayList<String> geneNames, int nMut, String mutPosFile) throws IOException{
		ArrayList<String> mutPosList = new BasicUtilityFunctions().readGeneNames(mutPosFile);
		//ArrayList<String> synoList = new BasicUtilityFunctions().readGeneNames(synoFile);
		STINode<AncestralInfo> par = this.nuPtree.getNode(parname);
		STINode<AncestralInfo> child = this.nuPtree.getNode(childname);
		System.out.printf("%s\t%s%n", parname, childname);
		for (int i = 0; i < nMut; i++){
			if (par.getData().sequence[i] != child.getData().sequence[i]){
				System.out.printf("%s\t%s\t%d\t%d\n", mutPosList.get(i), geneNames.get(i), par.getData().sequence[i], child.getData().sequence[i]);
			}
		}
	}
	
	/**
	 * 
	 * @param parname
	 * @param childname
	 * @param geneNames
	 * @param nMut
	 * @throws IOException
	 * Created On: Mar 9, 2018
	 */
	public void printMutationsP2C(String parname, String childname, ArrayList<String> geneNames, int nMut) throws IOException{
		STINode<AncestralInfo> par = this.nuPtree.getNode(parname);
		STINode<AncestralInfo> child = this.nuPtree.getNode(childname);
		System.out.printf("%s\t%s%n", parname, childname);
		for (int i = 0; i < nMut; i++){
			if (par.getData().sequence[i] != child.getData().sequence[i]){
				System.out.printf("%s\t%d\t%d\n", geneNames.get(i), par.getData().sequence[i], child.getData().sequence[i]);
			}
		}
	}
    
	/**
	 * @param args
	 * Created On: Mar 9, 2018
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		runCO5();
	}
	
	public static void runCO8() throws IOException{
		BasicUtilityFunctions BUF = new BasicUtilityFunctions();
		String varFileName = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO8/Complete/best_posterior_genotype.txt";
		String treefile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO8/Complete/MPEAR_mcc_tree.txt";
		String cellfile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO8/CO8_cellnames.txt";
		String genefile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO8/CO8.genenames";
		String expectedMatFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO8/Complete/MPEAR_mcc_expected_Mat.txt";
		double d = 0.101696;
		double l = 0.105672;
		double r = 0.05;
		InferBranchMutations IB = new InferBranchMutations(d, l, r, 0);
		IB.annotateMutations(treefile, varFileName, cellfile, genefile, expectedMatFile);
	}
	
	public static void runCO5() throws IOException{
		BasicUtilityFunctions BUF = new BasicUtilityFunctions();
		String varFileName = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO5/samples/posterior_genotypes.txt";
		String treefile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO5/samples/best_mcc_tree.txt";
		String cellfile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO5/CO5_cellnames.txt";
		String genefile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO5/CO5.genenames";
		String expectedMatFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO5/samples/MPEAR_mcc_expected_Mat.txt";
		double d = 0.095601;
		double l = 0.090232;
		double r = 0.05;
		InferBranchMutations IB = new InferBranchMutations(d, l, r, 0);
		IB.annotateMutations(treefile, varFileName, cellfile, genefile, expectedMatFile);
	}

}
