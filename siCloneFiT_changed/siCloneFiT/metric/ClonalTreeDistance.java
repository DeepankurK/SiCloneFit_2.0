/**
 * Nov 13, 2017
 */
package siCloneFiT.metric;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import SiFit.BasicUtilityFunctions;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeBipartition;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.PostTraversal;
import jeigen.DenseMatrix;
import siCloneFiT.simulation.SimulateClonalTree;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Nov 13, 2017
 */
public class ClonalTreeDistance {

	/**
	 * Find the most recent common ancestor of two nodes in a tree
	 * @param tree
	 * @param leaf1
	 * @param leaf2
	 * @return
	 * Created On: Nov 13, 2017
	 */
	public static STINode<Double> findMRCA(STITree<Double> tree, STINode<Double> leaf1, STINode<Double> leaf2){
		if (leaf1.isAncestor(leaf2))
			return leaf1;
		else if (leaf2.isAncestor(leaf1))
			return leaf2;
		else{
		int OUT = 1;
		int IN = 2;
		
		Map<STINode<Double>, Integer> marks = new HashMap<>();
		for (STINode<Double> node: tree.getNodes())
			marks.put(node, OUT);
		marks.put(tree.getRoot(), IN);
		STINode<Double> node = leaf1;
		while(!node.isRoot()) {
			marks.put(node, IN);
			node = node.getParent();
		}
		
		node = leaf2;
		while(!node.isRoot()) {
			marks.put(node, IN);
			node = node.getParent();
		}
		
		STINode<Double> mrcaRoot = tree.getRoot();
		while(!mrcaRoot.isLeaf()) {
			 int count = 0;
			 STINode<Double> inChild = null;
			 for (STINode<Double> child: mrcaRoot.getChildren()){
				 if(marks.get(child) == IN) {
					 count++;
		             inChild = child;
				 }
			 }
			 if(count == 1)
				 mrcaRoot = inChild;
			 else
				 break;
		}
//		System.out.println("mrca = " + mrcaRoot.getID());
		return mrcaRoot;
		}
	}
	
	/**
	 * Return the common ancestor for a set of leaves
	 * @param tree
	 * @param leaves
	 * @return
	 * Created On: Dec 18, 2018
	 */
	public static STINode<Double> findMRCA(STITree<Double> tree, ArrayList<STINode<Double>> leaves){
		if (leaves.size() == 1)
			return leaves.get(0);
		else{
			STINode<Double> mrca = findMRCA(tree, leaves.get(0), leaves.get(1));
			//		System.out.println("start "+mrca.getName());
			for (int i = 2; i < leaves.size(); i++){
				STINode<Double> i_mrca = findMRCA(tree, mrca, leaves.get(i));
				mrca = i_mrca;
				//			System.out.println(i +" " +mrca.getName());
			}
			return mrca;
		}
	}
	
	/**
	 * Find the distance between two nodes in a tree
	 * @param tree
	 * @param leaf1
	 * @param leaf2
	 * @return
	 * Created On: Nov 13, 2017
	 */
	public static double findPairwiseDist(STITree<Double> tree, STINode<Double> leaf1, STINode<Double> leaf2){
		STINode<Double> MRCA = findMRCA(tree, leaf1, leaf2);
		double dist = 0;
		
		STINode<Double> node = leaf1;
		
//		System.out.println(MRCA.getID());
		while(node.getID() != MRCA.getID()) {
			dist += node.getParentDistance();
//			System.out.println(dist);
//			System.out.println(node.isRoot());
			node = node.getParent();
//			System.out.println("hehe " + node.getID());
		}
		node = leaf2;
		
//		System.out.println("mrca = "+ MRCA.getID());
		while(node.getID() != MRCA.getID()) {
			dist += node.getParentDistance();
//			System.out.println(node.isRoot());
//			if (node.isRoot()){
//				System.out.println(leaf1.getID());
//				System.out.println(leaf2.getID());
//				for (STINode<Double> nc: node.getChildren())
//					System.out.println(nc.getID());
//			}
			node = node.getParent();
		}
		return dist;
	}
	
	/** 
	 * Replace internal branch lengths by 1 and leaf branch lengths by 0
	 * @param tree
	 * @return
	 * Created On: Nov 13, 2017
	 */
	public static STITree<Double> getIntegerBLTree(STITree<Double> tree){
		STITree<Double> nuTree = new STITree<>(tree);
		for (STINode<Double> n: nuTree.getNodes()){
			if (n.isLeaf())
				n.setParentDistance(0);
			else
				n.setParentDistance(1);
		}
		nuTree.getRoot().setParentDistance(0);
		return nuTree;
	}
	
	/**
	 * Find the pairwise cell shortest path distance
	 * @param trueTree
	 * @param inferredTree
	 * @param nCell
	 * @return
	 * Created On: Nov 13, 2017
	 */
	public static double getPairwiseCellSPTreeDist(STITree<Double> trueTree, STITree<Double> inferredTree, int nCell){
		ArrayList<String> cellList = new ArrayList<>();
		for (int i = 0; i < nCell; i++){
			cellList.add("sc" + i);
		}
		double dist = 0;
		for (int i = 0; i < nCell; i++){
			for (int j = i + 1; j < nCell; j++){
//				System.out.println("before tt " + i + " "+ j);
//				double trueTreeDist = findPairwiseDist(trueTree, trueTree.getNode(cellList.get(i)).getParent(), trueTree.getNode(cellList.get(j)).getParent());
				double trueTreeDist = findPairwiseDist(trueTree, trueTree.getNode(cellList.get(i)), trueTree.getNode(cellList.get(j)));
//				System.out.printf("cell %s and cell %s have truedist = %f\n", cellList.get(i), cellList.get(j), trueTreeDist);
//				System.out.println("cell i = " + cellList.get(i) + " cell j = "+ cellList.get(j));
//				double iTreeDist = findPairwiseDist(inferredTree, inferredTree.getNode(cellList.get(i)).getParent(), inferredTree.getNode(cellList.get(j)).getParent());
				double iTreeDist = findPairwiseDist(inferredTree, inferredTree.getNode(cellList.get(i)), inferredTree.getNode(cellList.get(j)));
//				System.out.printf("cell %s and cell %s have inferreddist = %f\n", cellList.get(i), cellList.get(j), iTreeDist);
				//				System.out.println(iTreeDist);
//				if (iTreeDist == Double.NEGATIVE_INFINITY){
//					System.out.println("cell i = " + cellList.get(i) + " cell j = "+ cellList.get(j));
//				}
//				System.out.println("after it " + i + " "+ j);
//				double trueTreeDist = findPairwiseDist(trueTree, trueTree.getNode(cellList.get(i)), trueTree.getNode(cellList.get(j)));
//				double iTreeDist = findPairwiseDist(inferredTree, inferredTree.getNode(cellList.get(i)), inferredTree.getNode(cellList.get(j)));

				dist += Math.abs(trueTreeDist - iTreeDist); 
			}
		}
		return dist;
	}
	
	/**
	 * Find the pairwise cell shortest path distance when the distance matrix is already computed for a clonal tree
	 * Applicable to SCITE
	 * @param trueTree
	 * @param clonalTreeDistMat
	 * @param nCell
	 * @return
	 * Created On: Dec 21, 2018
	 */
	public static double getPairwiseCellSPTreeDist(STITree<Double> trueTree, double[][] clonalTreeDistMat, int nCell){
		ArrayList<String> cellList = new ArrayList<>();
		for (int i = 0; i < nCell; i++){
			cellList.add("sc" + i);
		}
		double dist = 0;
		for (int i = 0; i < nCell; i++){
			for (int j = i + 1; j < nCell; j++){
				double trueTreeDist = findPairwiseDist(trueTree, trueTree.getNode(cellList.get(i)), trueTree.getNode(cellList.get(j)));
				double iTreeDist = clonalTreeDistMat[i][j];
				dist += Math.abs(trueTreeDist - iTreeDist);
			}
		}
		return dist;
	}
	
	/**
	 * Get the pairwise cell shortest path distance for SiFit tree
	 * @param trueTree
	 * @param inferredTree
	 * @param nCell
	 * @param cloneDist
	 * @param cellToCloneMap
	 * @return
	 * Created On: Dec 18, 2018
	 */
	public static double getPairwiseCellSPTreeDistSiFit(STITree<Double> trueTree, STITree<Double> inferredTree, int nCell, 
														double[][] cloneDist, HashMap<String, Integer> cellToCloneMap){
		ArrayList<String> cellList = new ArrayList<>();
		for (int i = 0; i < nCell; i++){
			cellList.add("sc" + i);
		}
		double dist = 0;
		for (int i = 0; i < nCell; i++){
			for (int j = i + 1; j < nCell; j++){
				double trueTreeDist = findPairwiseDist(trueTree, trueTree.getNode(cellList.get(i)), trueTree.getNode(cellList.get(j)));
				Integer cell_i_clone = cellToCloneMap.get(cellList.get(i));
				Integer cell_j_clone = cellToCloneMap.get(cellList.get(j));
				double iTreeDist = cloneDist[cell_i_clone-1][cell_j_clone-1];
				dist += Math.abs(trueTreeDist - iTreeDist); 
			}
		}
		return dist;
	}
	
	/**
	 * Find the pairwise cell shortest path distance, when SCG reports only two clones
	 * @param trueTree
	 * @param cellCloneArr
	 * @param nCell
	 * @return
	 * Created On: Mar 17, 2018
	 */
	public static double getPairwiseCellSPTreeDist(STITree<Double> trueTree, String[] cellCloneArr, int nCell){
		ArrayList<String> cellList = new ArrayList<>();
		for (int i = 0; i < nCell; i++){
			cellList.add("sc" + i);
		}
		double dist = 0;
		for (int i = 0; i < nCell; i++){
			for (int j = i + 1; j < nCell; j++){
				double trueTreeDist = findPairwiseDist(trueTree, trueTree.getNode(cellList.get(i)).getParent(), trueTree.getNode(cellList.get(j)).getParent());
				double iTreeDist = 0;
				if (!cellCloneArr[i].equals(cellCloneArr[j])){
					iTreeDist = 2;
//					System.out.println(iTreeDist);
				}
				dist += Math.abs(trueTreeDist - iTreeDist); 
			}
		}
		return dist;
	}
	
	/**
	 * Find the pairwise cell shortest path distance for the whole tree
	 * Doublets are not considered
	 * @param trueTree
	 * @param inferredTree
	 * @param nCell
	 * @param doubletFlagList
	 * @return
	 * Created On: Mar 15, 2018
	 */
	public static double getPairwiseCellSPTreeDist(STITree<Double> trueTree, STITree<Double> inferredTree, int nCell, ArrayList<Integer> doubletFlagList){
		ArrayList<String> cellList = new ArrayList<>();
		for (int i = 0; i < nCell; i++){
			if (doubletFlagList.get(i) == 0)
				cellList.add("sc" + i);
		}
		double dist = 0;
		for (int i = 0; i < cellList.size(); i++){
			for (int j = i + 1; j < cellList.size(); j++){
//				double trueTreeDist = findPairwiseDist(trueTree, trueTree.getNode(cellList.get(i)).getParent(), trueTree.getNode(cellList.get(j)).getParent());
				double trueTreeDist = findPairwiseDist(trueTree, trueTree.getNode(cellList.get(i)), trueTree.getNode(cellList.get(j)));
//				double iTreeDist = findPairwiseDist(inferredTree, inferredTree.getNode(cellList.get(i)).getParent(), inferredTree.getNode(cellList.get(j)).getParent());
				double iTreeDist = findPairwiseDist(inferredTree, inferredTree.getNode(cellList.get(i)), inferredTree.getNode(cellList.get(j)));
				dist += Math.abs(trueTreeDist - iTreeDist); 
			}
		}
		return dist;
	}
	
	/**
	 * Find the pairwise cell shortest path distance when cell names are given
	 * @param trueTree
	 * @param inferredTree
	 * @param cellNamesFile
	 * @param nCell
	 * @return
	 * @throws IOException
	 * Created On: Feb 6, 2018
	 */
	public static double getPairwiseCellSPTreeDist(STITree<Double> trueTree, STITree<Double> inferredTree, String cellNamesFile, int nCell) throws IOException{
		ArrayList<String> cellList = new ArrayList<>();
		BufferedReader brCells = new BufferedReader(new FileReader(cellNamesFile));
		String[] singleCellNameArray = brCells.readLine().split(" ");
		for (String s: singleCellNameArray){
			cellList.add(s);
        }
		brCells.close();
		double dist = 0;
		for (int i = 0; i < nCell; i++){
			for (int j = i + 1; j < nCell; j++){
				double trueTreeDist = findPairwiseDist(trueTree, trueTree.getNode(cellList.get(i)).getParent(), trueTree.getNode(cellList.get(j)).getParent());
				double iTreeDist = findPairwiseDist(inferredTree, inferredTree.getNode(cellList.get(i)).getParent(), inferredTree.getNode(cellList.get(j)).getParent());
//				double trueTreeDist = findPairwiseDist(trueTree, trueTree.getNode(cellList.get(i)), trueTree.getNode(cellList.get(j)));
//				double iTreeDist = findPairwiseDist(inferredTree, inferredTree.getNode(cellList.get(i)), inferredTree.getNode(cellList.get(j)));

				dist += Math.abs(trueTreeDist - iTreeDist); 
			}
		}
		return dist;
	}
	
	/**
	 * Get pairwise SP dist between SiFit tree and true tree
	 * @param trueFile
	 * @param inferredFile
	 * @param nCell
	 * @param BUF
	 * @return
	 * @throws IOException
	 * Created On: Dec 17, 2018
	 */
	public static double getPairwiseCellSPSiFitTreeDist(String trueFile, String inferredFile, String SiFitCloneFile, int nCell, BasicUtilityFunctions BUF) throws IOException{
		STITree<Double> trueTree = BUF.getTree(BUF.readNewickString(trueFile));
//		System.out.println(Arrays.toString(trueTree.getLeaves()));
		trueTree.rerootTreeAtEdge("sc0");
		STITree<Double> itrueTree = getIntegerBLTree(trueTree);
		STINode<Double> rT = itrueTree.getRoot();
		ArrayList<STINode<Double>> Trueroot_chileLeaves = new ArrayList<>();
		
		for (STINode<Double> rc: rT.getChildren()){
			if (rc.isLeaf())
				Trueroot_chileLeaves.add(rc);
		}
		if (Trueroot_chileLeaves.size()>0){
//			System.out.println("hello");
			rT.createChild("rlpT");
			STINode<Double> rlpT = itrueTree.getNode("rlpT");
			rlpT.setParentDistance(1);
			for (STINode<Double> rc: Trueroot_chileLeaves){
				rlpT.adoptChild(rc);
			}
		}
		
		// Analyze the clone file
		String SiFitCellClone = BUF.readNewickString(SiFitCloneFile);
		String[] SiFitCellCloneArr = SiFitCellClone.split(" ");
		HashMap<Integer, ArrayList<String>> cloneToCellMap = new HashMap<>();
		HashMap<String, Integer> cellToCloneMap = new HashMap<>();
		for (int i = 0; i < nCell; i++){
			Integer cloneID = Integer.parseInt(SiFitCellCloneArr[i]);
			cellToCloneMap.put("sc" + Integer.toString(i), cloneID);
			if (cloneToCellMap.containsKey(cloneID))
				cloneToCellMap.get(cloneID).add("sc" + Integer.toString(i));
			else{
				ArrayList<String> cloneCells = new ArrayList<>();
				cloneCells.add("sc" + Integer.toString(i));
				cloneToCellMap.put(cloneID, cloneCells);
			}
		}
//		System.out.println(cloneToCellMap);
//		System.out.println(cellToCloneMap);
		
//		double[][] cloneDist = new double[cloneToCellMap.size()][cloneToCellMap.size()];
		
		
		STITree<Double> inferredTree = BUF.getTree(BUF.readNewickString(inferredFile));
		// Rename the leaf nodes (sc0 to scm-1)
		for (int i = 1; i <= nCell; i++){
			STINode<Double> i_leaf = inferredTree.getNode("sc" + Integer.toString(i));
			i_leaf.setName("sc" + Integer.toString(i-1));
		}
//		System.out.println(Arrays.toString(inferredTree.getLeaves()));
		// Set the parent distance to 1 for whatever branch has mutations
		for (STINode<Double> node: inferredTree.getNodes()){
			if (node.getParentDistance() > 0.0){
				node.setParentDistance(1);
			}
//			if (!node.isLeaf()){
//				System.out.println(Arrays.toString(node.getTree().getLeaves()));
//			}
		}
		
		
		
		inferredTree.rerootTreeAtEdge("sc0");
		STINode<Double> r = inferredTree.getRoot();
		ArrayList<STINode<Double>> root_chileLeaves = new ArrayList<>();
		
		for (STINode<Double> rc: r.getChildren()){
			if (rc.isLeaf())
				root_chileLeaves.add(rc);
		}
		if (root_chileLeaves.size()>0){
			r.createChild("rlp");
			STINode<Double> rlp = inferredTree.getNode("rlp");
			rlp.setParentDistance(1);
			for (STINode<Double> rc: root_chileLeaves){
				rlp.adoptChild(rc);
			}
		}
		for (STINode<Double> node: inferredTree.getNodes()){
			if (node.isLeaf() == false){
				node.setName(Integer.toString(node.getID()));
			}
		}
		
		double[][] cloneDist = getClonalDistanceMatrix(cloneToCellMap, inferredTree);
//		double dist1 = getPairwiseCellSPTreeDist(itrueTree, inferredTree, nCell);
//		System.out.println(2*dist1/(nCell*(nCell-1)));
		
		double dist = getPairwiseCellSPTreeDistSiFit(itrueTree, inferredTree, nCell, cloneDist, cellToCloneMap);
		
		// checking if rooting at 0 works
//		nCell = nCell - 1;
		return (2*dist/(nCell*(nCell-1)));
	}
	
	public static double getPairwiseCellSPSiFitTreeDistNotRooted(String trueFile, String inferredFile, String SiFitCloneFile, int nCell, BasicUtilityFunctions BUF) throws IOException{
		STITree<Double> trueTree = BUF.getTree(BUF.readNewickString(trueFile));
		STITree<Double> itrueTree = getIntegerBLTree(trueTree);
//		System.out.println(itrueTree.toNewick());
		// Analyze the clone file
		String SiFitCellClone = BUF.readNewickString(SiFitCloneFile);
		String[] SiFitCellCloneArr = SiFitCellClone.split(" ");
		HashMap<Integer, ArrayList<String>> cloneToCellMap = new HashMap<>();
		HashMap<String, Integer> cellToCloneMap = new HashMap<>();
		for (int i = 0; i < nCell; i++){
			Integer cloneID = Integer.parseInt(SiFitCellCloneArr[i]);
			cellToCloneMap.put("sc" + Integer.toString(i), cloneID);
			if (cloneToCellMap.containsKey(cloneID))
				cloneToCellMap.get(cloneID).add("sc" + Integer.toString(i));
			else{
				ArrayList<String> cloneCells = new ArrayList<>();
				cloneCells.add("sc" + Integer.toString(i));
				cloneToCellMap.put(cloneID, cloneCells);
			}
		}
//		System.out.println(cloneToCellMap);
//		System.out.println(cellToCloneMap);
		
		STITree<Double> inferredTree = BUF.getTree(BUF.readNewickString(inferredFile));
		// Rename the leaf nodes (sc0 to scm-1)
		for (int i = 1; i <= nCell; i++){
			STINode<Double> i_leaf = inferredTree.getNode("sc" + Integer.toString(i));
			i_leaf.setName("sc" + Integer.toString(i-1));
		}
		
		// Set the parent distance to 1 for whatever branch has mutations
		for (STINode<Double> node: inferredTree.getNodes()){
			if (node.getParentDistance() > 0.0){
				node.setParentDistance(1);
			}			
		}
		double[][] cloneDist = getClonalDistanceMatrix(cloneToCellMap, inferredTree);
		
//		System.out.println(inferredTree.toNewick());
//		System.out.println(new DenseMatrix(cloneDist));
//		double dist1 = getPairwiseCellSPTreeDist(itrueTree, inferredTree, nCell);
//		System.out.println(2*dist1/(nCell*(nCell-1)));
		double dist = getPairwiseCellSPTreeDistSiFit(itrueTree, inferredTree, nCell, cloneDist, cellToCloneMap);
		
		return (2*dist/(nCell*(nCell-1)));
	}
	
	/**
	 * Compute the distance matrix for all the clones
	 * @param cloneToCellMap
	 * @param tree
	 * @return
	 * Created On: Dec 18, 2018
	 */
	public static double[][] getClonalDistanceMatrix(HashMap<Integer, ArrayList<String>> cloneToCellMap, STITree<Double> tree){
		for (Integer c: cloneToCellMap.keySet()){
			ArrayList<String> clone = cloneToCellMap.get(c);
			resizeCloneInternalBranches(clone, tree);
		}
		
		double[][] cloneDist = new double[cloneToCellMap.size()][cloneToCellMap.size()];
		int nClone = cloneToCellMap.size();
		for (int i = 1; i <= nClone; i++){
			for (int j = i+1; j <= nClone; j++){
				double dist_i_j = getMinClonalDistance(cloneToCellMap.get(i), cloneToCellMap.get(j), tree);
				cloneDist[i-1][j-1] = dist_i_j;
				cloneDist[j-1][i-1] = dist_i_j;
			}
		}
		return cloneDist;
	}
	
	/**
	 * set the branch length of branches that connect all the cells in a clone to 0
	 * @param clone
	 * @param tree
	 * Created On: Dec 19, 2018
	 */
	public static void resizeCloneInternalBranches(ArrayList<String> clone, STITree<Double> tree){
		ArrayList<STINode<Double>> cloneLeaves = new ArrayList<>();
		for (String leaf: clone)
			cloneLeaves.add(tree.getNode(leaf));
		STINode<Double> cloneRoot = findMRCA(tree, cloneLeaves);
		for (STINode<Double> node: cloneLeaves){
			while(node.getID() != cloneRoot.getID()) {
				node.setParentDistance(0);
				node = node.getParent();
			}
		}
	}
	
	/**
	 * Compute the minimum tree distance between two clones
	 * Minimum of the distances between the member cells 
	 * @param clone1
	 * @param clone2
	 * @param tree
	 * @return
	 * Created On: Dec 18, 2018
	 */
	public static double getMinClonalDistance(ArrayList<String> clone1, ArrayList<String> clone2, STITree<Double> tree){
		double minDist = Double.MAX_VALUE;
		for (String ci: clone1){
			for (String cj: clone2){
				double ci_cj_dist = findPairwiseDist(tree, tree.getNode(ci), tree.getNode(cj));
				if (ci_cj_dist == Double.NEGATIVE_INFINITY)
					ci_cj_dist = Double.POSITIVE_INFINITY;
				if (ci_cj_dist < minDist)
					minDist = ci_cj_dist;
			}
		}
		return minDist;
	}
	
	/**
	 * Find the pairwise SP distance for SCITE tree
	 * @param trueFile
	 * @param SCITE_cloneTreeDistFile
	 * @param nCell
	 * @param BUF
	 * @return
	 * @throws IOException
	 * Created On: Dec 21, 2018
	 */
	public static double getSCITETreePairwiseCellSPDist(String trueFile, String SCITE_cloneTreeDistFile, int nCell, BasicUtilityFunctions BUF) throws IOException{
		STITree<Double> trueTree = BUF.getTree(BUF.readNewickString(trueFile));
		STITree<Double> itrueTree = getIntegerBLTree(trueTree);
		double[][] SCITECloneTreeDist = new double[nCell][nCell];
		try (BufferedReader br = new BufferedReader(new FileReader(SCITE_cloneTreeDistFile))) {
			String line;
			int c = 0;
		    while ((line = br.readLine()) != null) {
		    	String[] row = line.split(" ");
		    	for (int i = 0; i < nCell; i++){
		    		SCITECloneTreeDist[c][i] = Double.parseDouble(row[i]);
		    	}
		    	c++;
		    }
		}
//		System.out.println(new DenseMatrix(SCITECloneTreeDist));
		double dist = getPairwiseCellSPTreeDist(itrueTree, SCITECloneTreeDist, nCell);
		return (2*dist/(nCell*(nCell-1)));
	}
	
	
	
	
	/**
	 * Find the pairwise cell shortest path distance
	 * @param trueFile
	 * @param inferredFile
	 * @param nCell
	 * @param BUF
	 * @return
	 * @throws IOException
	 * Created On: Nov 13, 2017
	 */
	public static double getPairwiseCellSPTreeDist(String trueFile, String inferredFile, int nCell, BasicUtilityFunctions BUF) throws IOException{
		STITree<Double> trueTree = BUF.getTree(BUF.readNewickString(trueFile));
		trueTree.rerootTreeAtEdge("sc0");
//		System.out.println(trueTree.toNewick());
//		STINode<Double> root = trueTree.getBipartitionClusters(arg0, arg1)
		STITree<Double> itrueTree = getIntegerBLTree(trueTree);
		
		STINode<Double> rT = itrueTree.getRoot();
		ArrayList<STINode<Double>> Trueroot_chileLeaves = new ArrayList<>();
		
		for (STINode<Double> rc: rT.getChildren()){
			if (rc.isLeaf())
				Trueroot_chileLeaves.add(rc);
		}
//		System.out.println(root_chileLeaves.size());
		
		if (Trueroot_chileLeaves.size()>0){
//			System.out.println("hello");
			rT.createChild("rlpT");
			STINode<Double> rlpT = itrueTree.getNode("rlpT");
			rlpT.setParentDistance(1);
			for (STINode<Double> rc: Trueroot_chileLeaves){
				rlpT.adoptChild(rc);
			}
		}
//		STINode<Double> rlpT = itrueTree.getNode("rlpT");
//		System.out.println(rlpT.getChildCount());
//		System.out.println(itrueTree.toNewick());
		
		STITree<Double> inferredTree = BUF.getTree(BUF.readNewickString(inferredFile));
		inferredTree.rerootTreeAtEdge("sc0");
//		System.out.println(inferredTree.toNewick());
		STITree<Double> iinferredTree = getIntegerBLTree(inferredTree);
		STINode<Double> r = iinferredTree.getRoot();
		ArrayList<STINode<Double>> root_chileLeaves = new ArrayList<>();
		
		for (STINode<Double> rc: r.getChildren()){
			if (rc.isLeaf())
				root_chileLeaves.add(rc);
		}
//		System.out.println(root_chileLeaves.size());
		
		if (root_chileLeaves.size()>0){
			r.createChild("rlp");
			STINode<Double> rlp = iinferredTree.getNode("rlp");
			rlp.setParentDistance(1);
			for (STINode<Double> rc: root_chileLeaves){
				rlp.adoptChild(rc);
			}
		}
		STINode<Double> rlp = iinferredTree.getNode("rlp");
		for (STINode<Double> node: iinferredTree.getNodes()){
			if (node.isLeaf() == false){
				node.setName(Integer.toString(node.getID()));
			}
		}
		
//		System.out.println(rlp.getChildCount());
//		STITree<Double> iinferredTree2 = BUF.getTree(iinferredTree.toNewick()); 
//		System.out.println(iinferredTree.toNewick());
//		System.out.println(iinferredTree.isRooted());
//		System.out.println(iinferredTree.getLeafCount());
		double dist = getPairwiseCellSPTreeDist(itrueTree, iinferredTree, nCell);
		
		// checking if rooting at 0 works
//		nCell = nCell - 1;
		return (2*dist/(nCell*(nCell-1)));
	}
	
	
	/**
	 * Find the pairwise cell shortest path distance for the whole inferred Tree
	 * without considering the doublets, only singlets are considered
	 * @param trueFile
	 * @param inferredFile
	 * @param nCell
	 * @param inferredDoubletFile
	 * @param BUF
	 * @return
	 * @throws IOException
	 * Created On: Mar 15, 2018
	 */
	public static double getPairwiseCellSPTreeDistSinglets(String trueFile, String inferredFile, int nCell, String inferredDoubletFile, BasicUtilityFunctions BUF) throws IOException{
		STITree<Double> trueTree = BUF.getTree(BUF.readNewickString(trueFile));
		trueTree.rerootTreeAtEdge("sc0");
		STITree<Double> itrueTree = getIntegerBLTree(trueTree);
		STITree<Double> inferredTree = BUF.getTree(BUF.readNewickString(inferredFile));
		inferredTree.rerootTreeAtEdge("sc0");
		STITree<Double> iinferredTree = getIntegerBLTree(inferredTree);
		String[] doubletFlagArr = BUF.readNewickString(inferredDoubletFile).split(" ");
		ArrayList<Integer> doubletFlagList = new ArrayList<>();
		int singletCount = 0;
		for (String f : doubletFlagArr){
			doubletFlagList.add(Integer.parseInt(f));
			if (f.equals("0"))
				singletCount++;
		}
//		System.out.println(singletCount);
		double dist = getPairwiseCellSPTreeDist(itrueTree, iinferredTree, nCell, doubletFlagList);
		return (2*dist/(singletCount*(singletCount-1)));
	}
	
	/**
	 * Find the pairwise cell shortest path distance
	 * @param trueFile
	 * @param inferredTree
	 * @param nCell
	 * @param SCF
	 * @return
	 * @throws IOException
	 * Created On: Feb 14, 2018
	 */
	public static double getPairwiseCellSPTreeDist(String trueFile, STITree<Double> inferredTree, int nCell, SCFUtilityFunctions SCF) throws IOException{
		STITree<Double> trueTree = SCF.getTree(SCF.readNewickString(trueFile));
		STITree<Double> itrueTree = getIntegerBLTree(trueTree);
		STITree<Double> iinferredTree = getIntegerBLTree(inferredTree);
		double dist = getPairwiseCellSPTreeDist(itrueTree, iinferredTree, nCell);
		return (2*dist/(nCell*(nCell-1)));
	}
	/**
	 * Get the SCG tree
	 * @param cloneTreeFile
	 * @param cellCloneFile
	 * @param BUF
	 * @return
	 * @throws IOException
	 * Created On: Nov 13, 2017
	 */
    public static STITree<Double> getSCGCellLineageTree(String cloneTreeFile, String cellCloneFile, BasicUtilityFunctions BUF) throws IOException{
    	STITree<Double> SCGTree = BUF.getTree(BUF.readNewickString(cloneTreeFile));
    	String cellCloneS = BUF.readNewickString(cellCloneFile);
    	String[] cellClone = cellCloneS.split(" ");
//    	System.out.println(Arrays.toString(cellClone));
    	for (int i = 0; i < cellClone.length; i++){
    		String cell = "sc" + i;
    		STINode<Double> cloneNode = SCGTree.getNode(cellClone[i]);
//    		System.out.println(cloneNode.getName());
    		cloneNode.createChild(cell);
    	}
    	return SCGTree;
    }
	
    /**
     * Find pairwise cell shortest path distance for SCG tree
     * @param trueFile
     * @param cloneTreeFile
     * @param cellCloneFile
     * @param nCell
     * @param BUF
     * @return
     * @throws IOException
     * Created On: Nov 13, 2017
     */
    public static double getSCGTreePairwiseCellSPDist(String trueFile, String cloneTreeFile, String cellCloneFile, int nCell, BasicUtilityFunctions BUF) throws IOException{
    	STITree<Double> trueTree = BUF.getTree(BUF.readNewickString(trueFile));
    	trueTree.rerootTreeAtEdge("sc0");
		STITree<Double> itrueTree = getIntegerBLTree(trueTree);
		
    	STITree<Double> SCGTree = getSCGCellLineageTree(cloneTreeFile, cellCloneFile, BUF);
    	SCGTree.rerootTreeAtEdge("sc0");
    	STITree<Double> iinferredTree = getIntegerBLTree(SCGTree);
//    	iinferredTree.rerootTreeAtEdge("sc0");
//    	System.out.println(iinferredTree.toNewick());
    	double dist = getPairwiseCellSPTreeDist(itrueTree, iinferredTree, nCell);
		return (2*dist/(nCell*(nCell-1)));
//    	return getPairwiseCellSPTreeDist(itrueTree, iinferredTree, nCell);
    }
    
    public static void writeSCGTree(String cloneTreeFile, String cellCloneFile, String cellTreeFile, BasicUtilityFunctions BUF) throws IOException{
    	STITree<Double> SCGTree = getSCGCellLineageTree(cloneTreeFile, cellCloneFile, BUF);
    	STITree<Double> iinferredTree = getIntegerBLTree(SCGTree);
    	BUF.writeNewickTree(cellTreeFile, iinferredTree.toNewick());
    }
    
    /**
     * Find pairwise cell shortest path distance for SCG tree when only two clones are reported
     * @param trueFile
     * @param cellCloneFile
     * @param nCell
     * @param BUF
     * @return
     * @throws IOException
     * Created On: Mar 17, 2018
     */
    public static double getSCGTreePairwiseCellSPDist(String trueFile, String cellCloneFile, int nCell, BasicUtilityFunctions BUF) throws IOException{
    	STITree<Double> trueTree = BUF.getTree(BUF.readNewickString(trueFile));
		STITree<Double> itrueTree = getIntegerBLTree(trueTree);
		String[] cellCloneArr = BUF.readNewickString(cellCloneFile).split(" ");
		double dist = getPairwiseCellSPTreeDist(itrueTree, cellCloneArr, nCell);
		return (2*dist/(nCell*(nCell-1)));
    }
    
    /**
     * Find pairwise cell shortest path distance for SCG tree (only for singlets)
     * Doublets are not considered
     * @param trueFile
     * @param cloneTreeFile
     * @param cellCloneFile
     * @param nCell
     * @param inferredDoubletFile
     * @param BUF
     * @return
     * @throws IOException
     * Created On: Mar 15, 2018
     */
    public static double getSCGTreePairwiseCellSPDistSinglets(String trueFile, String cloneTreeFile, String cellCloneFile, int nCell, String inferredDoubletFile, BasicUtilityFunctions BUF) throws IOException{
    	STITree<Double> trueTree = BUF.getTree(BUF.readNewickString(trueFile));
    	trueTree.rerootTreeAtEdge("sc0");
		STITree<Double> itrueTree = getIntegerBLTree(trueTree);
    	STITree<Double> SCGTree = getSCGCellLineageTree(cloneTreeFile, cellCloneFile, BUF);
    	SCGTree.rerootTreeAtEdge("sc0");
    	STITree<Double> iinferredTree = getIntegerBLTree(SCGTree);
    	String[] doubletFlagArr = BUF.readNewickString(inferredDoubletFile).split(" ");
		ArrayList<Integer> doubletFlagList = new ArrayList<>();
		int singletCount = 0;
		for (String f : doubletFlagArr){
			doubletFlagList.add(Integer.parseInt(f));
			if (f.equals("0"))
				singletCount++;
		}
//		System.out.println(singletCount);
		double dist = getPairwiseCellSPTreeDist(itrueTree, iinferredTree, nCell, doubletFlagList);
		return (2*dist/(singletCount*(singletCount-1)));
    }
    
    /**
     * 
     * @param oncoNEMcloneTree
     * @param cellCloneFile
     * @param BUF
     * @return
     * @throws IOException
     * Created On: Nov 22, 2017
     */
    public static STITree<Double> getOncoNEMCellTree(STITree<Double> oncoNEMcloneTree, String cellCloneFile, BasicUtilityFunctions BUF) throws IOException{
    	String cellCloneS = BUF.readNewickString(cellCloneFile);
    	String[] cellClone = cellCloneS.split(" ");
    	for (int i = 0; i < cellClone.length; i++){
    		String cell = "sc" + i;
    		STINode<Double> cloneNode = oncoNEMcloneTree.getNode(cellClone[i]);
    		cloneNode.createChild(cell);
    	}
    	// Remove empty clones (now present as leaves)
    	for (String nodeName: oncoNEMcloneTree.getLeaves()){
    		STINode<Double> node = oncoNEMcloneTree.getNode(nodeName);
    		if (node.getName().startsWith("sc") == false){
    			// This is an empty clone leaf
    			STINode<Double> nodePar = node.getParent();
    			STINode<Double> nodeGrandPar = nodePar.getParent();
    			STINode<Double> nodeSib = BUF.getOtherChild(nodePar, node);
    			nodeGrandPar.adoptChild(nodeSib);
    			nodePar.removeAllChildren();
    			nodeGrandPar.removeChild(nodePar, true);
    		}
    	}
    	return oncoNEMcloneTree;
    }
    
    /**
     * Find pairwise cell shortest path distance for OncoNEM tree
     * @param trueFile
     * @param cloneTreeFile
     * @param cellCloneFile
     * @param nCell
     * @param BUF
     * @return
     * @throws IOException
     * Created On: Nov 22, 2017
     */
    public static double getOncoNEMTreePairwiseCellSPDist(String trueFile, String cloneTreeFile, String cellCloneFile, int nCell, BasicUtilityFunctions BUF) throws IOException{
    	STITree<Double> trueTree = BUF.getTree(BUF.readNewickString(trueFile));
    	trueTree.rerootTreeAtEdge("sc0");
		STITree<Double> itrueTree = getIntegerBLTree(trueTree);
		STITree<Double> oncoNEMcloneTree = BUF.getTree(BUF.readNewickString(cloneTreeFile));
    	STITree<Double> OncoNEMtree = getOncoNEMCellTree(oncoNEMcloneTree, cellCloneFile, BUF);
    	OncoNEMtree.rerootTreeAtEdge("sc0");
    	STITree<Double> iinferredTree = getIntegerBLTree(OncoNEMtree);
//    	System.out.println(iinferredTree.toNewick());
    	double dist = getPairwiseCellSPTreeDist(itrueTree, iinferredTree, nCell);
		return (2*dist/(nCell*(nCell-1)));
    }
    
    public static void writeOncoNEMTree(String cloneTreeFile, String cellCloneFile, String cellTreeFile, BasicUtilityFunctions BUF) throws IOException{
    	STITree<Double> oncoNEMcloneTree = BUF.getTree(BUF.readNewickString(cloneTreeFile));
    	STITree<Double> OncoNEMtree = getOncoNEMCellTree(oncoNEMcloneTree, cellCloneFile, BUF);
    	STITree<Double> iinferredTree = getIntegerBLTree(OncoNEMtree);
    	BUF.writeNewickTree(cellTreeFile, iinferredTree.toNewick());
    }
    
	public static ArrayList<STINode<Double>> findAncestors(STITree<Double> tree, STINode<Double> leaf1){
		ArrayList<STINode<Double>> leaf1Ancestors = new ArrayList<>();
		STINode<Double> leaf1Par = leaf1.getParent();
		leaf1Ancestors.add(leaf1Par);
		while(leaf1Par.getID() != tree.getRoot().getID()){
			STINode<Double> node = leaf1Par.getParent();
			leaf1Par = node;
		}
		return leaf1Ancestors;
	}
	
	public static STITreeBipartition getBipartition(String trueFile, BasicUtilityFunctions BUF) throws IOException{
		STITree<Double> trueTree = BUF.getTree(BUF.readNewickString(trueFile));
		STINode<Double> root = trueTree.getRoot();
		root.setName("R");
		PostTraversal<Double> traversal = new PostTraversal<Double>(root);
		ArrayList<STINode<Double>> rootC = new ArrayList<>();
		for (STINode<Double> node: root.getChildren())
			rootC.add(node);
		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		List<STITreeBipartition> bipartitions = new LinkedList<STITreeBipartition>();
		String[] leaves = new String[trueTree.getLeafCount()];
		int j = 0;
		for (TNode node : trueTree.getNodes()) {
			if (node.isLeaf()) {
				leaves[j] = node.getName();
				j++;
			}
		}
		for (TNode node : traversal) {
			BitSet bs = new BitSet();
			
			if (node.isLeaf()) {	// Create a bipartition with a single leaf.
				// Find the index of this leaf.
				int i = 0;
				
				for (i = 0; i < leaves.length; i++) {
					if (node.getName().equals(leaves[i])) {
						break;
					}
				}
				
				// The leaf must always be found.
				assert(i < leaves.length);
				
				bs.set(i);				
				map.put(node, bs);
			}
			else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = map.get(child);
					bs.or(childCluster);
				}
				
				map.put(node, bs);
			}
			
			// Record only nontrivial, non-duplicate bipartitions.			
			if (bs.cardinality() > 1 && bs.cardinality() < leaves.length - 1) {			
				boolean duplicate = false;
				
				STITreeBipartition tb = new STITreeBipartition(leaves);
				tb.setBipartition(bs, null);
				
				for (int i = 0; i < bipartitions.size(); i++) {
					if (tb.isEqual(bipartitions.get(i))) {
						duplicate = true;
						break;
					}
				}
				
				if (!duplicate) {
					bipartitions.add(tb);
				}
			}
		}
		BitSet rootbs = map.get(root);
		STITreeBipartition tbR = new STITreeBipartition(leaves);
		tbR.setBipartition(rootbs, rootbs);
		return tbR;
//		System.out.println(Arrays.toString(rootC.get(0).getLeaves()));
//		for (TNode node : traversal){
//			System.out.println(node.getName() + " "+ node.getID());
//		}
//		return null;
	}
	/**
	 * @param args
	 * Created On: Nov 13, 2017
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		BasicUtilityFunctions BUF = new BasicUtilityFunctions();
		
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
		ArrayList<String> cloneNames = new ArrayList<>();
		for (int i = 0; i < 16; i++)
			cloneNames.add("C" + Integer.toString(i));
//		STITree<Double> lt = getLinearTree(8, 1, 1.1, SCF);
//		STITree<Double> lt = SimulateClonalTree.getBranchingTree(cloneNames, 10, 1, 1.1, SCF);
//		for (STINode<Double> node: lt.getNodes()){
////			System.out.println(node.getName());
//			if (node.getName() == "")
//				node.setName("i"+Integer.toString(node.getID()));
//		}
		
//		String tt = "(C9:0.021963236979083638,((((C8:0.07115438758940615,C7:0.017521777979975343)i5:0.056694850156808065,C6:0.11966369790021499)i4:0.03510439805850386,((C5:0.06362539598768709,(C4:0.09391876969591087,C3:0.007371840039728682)i12:0.08830497801439058)i10:0.01186608972770352,C2:0.11966379723932288)i9:0.1385332139337697)i3:0.004841670886448254,(C1:0.026467998174073158,C0:0.09560100386196492)i16:0.1003336183637787)i2:0.027369275411229802)i0;";
//		STITree<Double> lt = SCF.getTree(tt);
//		
//		for (STINode<Double> node: lt.getNodes()){
////			System.out.println(node.getName());
//			if (node.getName() == "")
//				node.setName("i"+Integer.toString(node.getID()));
//		}
//		System.out.println(lt.toNewick());
//		ArrayList<STINode<Double>> leaves = new ArrayList<>();
//		leaves.add(lt.getNode("C2"));
//		leaves.add(lt.getNode("C3"));
//		leaves.add(lt.getNode("C9"));
//		leaves.add(lt.getNode("C6"));
//		
//		STINode<Double> mrca = findMRCA(lt, leaves);
//		System.out.println(mrca.getName());
		
//		STINode<Double> mm = findMRCA(lt, lt.getNode("C4"), lt.getNode("i9"));
//		System.out.println("see " + mm.getName());
		
//		String tt = "((1,2),(3,((4,5,6),(7,8))));";
		
//		String tt = "(C0,(C1,(C2,C3)));";
//		STITree<Double> ttr = BUF.getTree(tt);
//		ttr.getNode("C0").createChild("sc0");
//		ttr.getNode("C0").createChild("sc1");
//		ttr.getNode("C1").createChild("sc2");
//		ttr.getNode("C2").createChild("sc3");
//		ttr.getNode("C2").createChild("sc4");
//		ttr.getNode("C2").createChild("sc5");
//		ttr.getNode("C3").createChild("sc6");
//		ttr.getNode("C3").createChild("sc7");
//		
//		STITree<Double> ittr = ClonalTreeDistance.getIntegerBLTree(ttr);
//		System.out.println(ittr.toNewick());
//		
//		String ti = "(C0,(C1,C2));";
//		STITree<Double> tti = BUF.getTree(ti);
//		tti.getNode("C0").createChild("sc0");
//		tti.getNode("C0").createChild("sc1");
//		tti.getNode("C1").createChild("sc2");
//		tti.getNode("C1").createChild("sc3");
//		tti.getNode("C1").createChild("sc4");
//		tti.getNode("C1").createChild("sc5");
//		tti.getNode("C2").createChild("sc6");
//		tti.getNode("C2").createChild("sc7");
//		
//		STITree<Double> itti = ClonalTreeDistance.getIntegerBLTree(tti);
//		System.out.println(itti.toNewick());
//		
//		System.out.println(ClonalTreeDistance.getPairwiseCellSPTreeDist(ittr, itti, 8));
//		
//		String trueFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/set2/dataset7/Orig_tree_dataset7.txt";
//		String SCFTree = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/set2/dataset7/SCF_tree.txt";
//		
//		double dist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueFile, SCFTree, 100, BUF);
//		System.out.println("dist = " + dist);
		
//		String scgTreeFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/set2/dataset7/SCG_results/scg_clone_tree.txt";
//		String scgCellCloneFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/set2/dataset7/SCG_results/scg_cell_clusters.txt";
//		double scgDist = ClonalTreeDistance.getSCGTreePairwiseCellSPDist(trueFile, scgTreeFile, scgCellCloneFile, 100, BUF);
//		System.out.println("SCG dist = " + scgDist);
		
		/*
		 * Directory Information
		 */
//		String dir = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/fnRateExpt/fn4/SiFit_n_SCITE/";
//		String dir = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/set5/";
		String dir = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/100_cells/400_sites-1/10_clones/";
//		String dir = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/fsTest/isa/100_sites/";
		
		
		
		
		// Doublet
//		String dir = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/10pDoublet/500_cells/100_sites/10_clones/";
		
		
		/*
		 * Singlet Datasets
		 */
		
		String[] dists = new String[10];
		int nCell = 100;
		
		/**
		 * SCITE results
		 */
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//			String SCITE_cloneTreeDistFile = dir + dataset + "/SCITE_rslts/SCITE_" + dataset + "_cloneTree_dist_mat.txt";
//			double dist = ClonalTreeDistance.getSCITETreePairwiseCellSPDist(trueFile, SCITE_cloneTreeDistFile, nCell, BUF);
//			System.out.println(dist);
//		}
		
		/**
		 * SiFit results
		 */
		
//		String[] nonRootedDists = new String[10];
//		System.out.println("not rooted at sc0");
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//			String SiFitTree = dir + dataset + "/SiFit_" + dataset + "_mutation_tree.newick";
//			String SiFitCloneFile = dir + dataset + "/SiFit_" + dataset + "_mutationTree_clustering.txt";
//			double dist = ClonalTreeDistance.getPairwiseCellSPSiFitTreeDistNotRooted(trueFile, SiFitTree, SiFitCloneFile, nCell, BUF);
//			System.out.println(dist);
//			nonRootedDists[i-1] = Double.toString(dist);
//		}
//		System.out.println(" ");
//		String SiFitDist = String.join(" ", nonRootedDists);
//		System.out.println(SiFitDist);
//		
//		System.out.println(" ");
//		System.out.println("rooted at sc0");
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//			String SiFitTree = dir + dataset + "/SiFit_" + dataset + "_mutation_tree.newick";
//			String SiFitCloneFile = dir + dataset + "/SiFit_" + dataset + "_mutationTree_clustering.txt";
//			double dist = ClonalTreeDistance.getPairwiseCellSPSiFitTreeDist(trueFile, SiFitTree, SiFitCloneFile, nCell, BUF);
//			System.out.println(dist);
//			dists[i-1] = Double.toString(dist);
//		}
//		System.out.println(" ");
//		String SCFDist = String.join(" ", dists);
//		System.out.println(SCFDist);
		
		
		
		
		/**
		 * SCG tree rooted results
		 */
		
		
//		System.out.println("Root at sc0");
//		System.out.println("SCG Results");
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//			
//			// No missing data
//			String scgCloneTreeFile = dir + dataset + "/SCG_results/scg_clone_tree.txt";
//			String scgCellTreeFile = dir + dataset + "/SCG_results/scg_cell_tree.txt";
//			String scgCellCloneFile = dir + dataset + "/SCG_results/scg_cell_clusters.txt";
//			
//			// 15% missing data
////			String scgCloneTreeFile = dir + dataset + "/15p_SCG_results/scg_clone_tree.txt";
////			String scgCellTreeFile = dir + dataset + "/15p_SCG_results/scg_cell_tree.txt";
////			String scgCellCloneFile = dir + dataset + "/15p_SCG_results/scg_cell_clusters.txt";
//			
//			// 30% missing data
////			String scgCloneTreeFile = dir + dataset + "/30p_SCG_results/scg_clone_tree.txt";
////			String scgCellTreeFile = dir + dataset + "/30p_SCG_results/scg_cell_tree.txt";
////			String scgCellCloneFile = dir + dataset + "/30p_SCG_results/scg_cell_clusters.txt";
//			
//			ClonalTreeDistance.writeSCGTree(scgCloneTreeFile, scgCellCloneFile, scgCellTreeFile, BUF);
//			double scgDist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueFile, scgCellTreeFile, nCell, BUF);
////			double scgDist = ClonalTreeDistance.getSCGTreePairwiseCellSPDist(trueFile, scgCloneTreeFile, scgCellCloneFile, nCell, BUF);;
//			System.out.println(scgDist);
//			dists[i-1] = Double.toString(scgDist);
//		}
//		System.out.println(" ");
//		String SCGDist = String.join(" ", dists);
//		System.out.println(SCGDist);
//		
		/**
		 * SCF tree rooted results
		 */
//		/**
		System.out.println(" ");
		System.out.println("MCC tree (root at sc0)");
		System.out.println("SCF Results");
		for (int i = 10; i <= 10; i++){
			String dataset = "dataset" + i;
			System.out.println(dataset);
			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
			
			// No missing data
			String SCFTree = dir + dataset + "/samples/best/best_mcc_tree.txt";
			
			// 15% missing data
//			String SCFTree = dir + dataset + "/15p_missing_samples/best/best_mcc_tree.txt";
			
			// 15% missing data
//			String SCFTree = dir + dataset + "/30p_missing_samples/best/best_mcc_tree.txt";
			
			double dist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueFile, SCFTree, nCell, BUF);
			System.out.println(dist);
			dists[i-1] = Double.toString(dist);
		}
		System.out.println(" ");
		String SCFDist = String.join(" ", dists);
		System.out.println(SCFDist);
		
		/*
		 * Doublet datasets
		 */
		
//		/**
//		 * SCG tree rooted results
//		 */
//		
//		String[] dists = new String[10];
//		int nCell = 500;
//		System.out.println("Root at sc0");
//		System.out.println("SCG Results");
//		for (int i = 5; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//			
//			// No missing data
////			String scgCloneTreeFile = dir + dataset + "/SCG_results/scg_clone_tree.txt";
////			String scgCellTreeFile = dir + dataset + "/SCG_results/scg_cell_tree.txt";
////			String scgCellCloneFile = dir + dataset + "/SCG_results/scg_cell_clusters.txt";
////			String doubletFile = dir + dataset + "/SCG_results/scg_inferred_doubletFlags.txt";
//			
//			// 15% missing data
////			String scgCloneTreeFile = dir + dataset + "/15p_SCG_results/scg_clone_tree.txt";
////			String scgCellTreeFile = dir + dataset + "/15p_SCG_results/scg_cell_tree.txt";
////			String scgCellCloneFile = dir + dataset + "/15p_SCG_results/scg_cell_clusters.txt";
////			String doubletFile = dir + dataset + "/15p_SCG_results/scg_inferred_doubletFlags.txt";
//			
//			// 30% missing data
//			String scgCloneTreeFile = dir + dataset + "/30p_SCG_results/scg_clone_tree.txt";
//			String scgCellTreeFile = dir + dataset + "/30p_SCG_results/scg_cell_tree.txt";
//			String scgCellCloneFile = dir + dataset + "/30p_SCG_results/scg_cell_clusters.txt";
//			String doubletFile = dir + dataset + "/30p_SCG_results/scg_inferred_doubletFlags.txt";
//			
//			ClonalTreeDistance.writeSCGTree(scgCloneTreeFile, scgCellCloneFile, scgCellTreeFile, BUF);
////			double scgDist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueFile, scgCellTreeFile, nCell, BUF);
////			double scgDist = ClonalTreeDistance.getSCGTreePairwiseCellSPDist(trueFile, scgCloneTreeFile, scgCellCloneFile, nCell, BUF);
//			double scgDist = ClonalTreeDistance.getSCGTreePairwiseCellSPDistSinglets(trueFile, scgCloneTreeFile, scgCellCloneFile, nCell, doubletFile, BUF);
//			System.out.println(scgDist);
//			dists[i-1] = Double.toString(scgDist);
//		}
//		System.out.println(" ");
//		String SCGDist = String.join(" ", dists);
//		System.out.println(SCGDist);
		
//		/**
//		 * SCF tree rooted results
//		 */
//		System.out.println(" ");
//		System.out.println("MCC tree (root at sc0)");
//		System.out.println("SCF Results");
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//			
//			// No missing data
////			String SCFTree = dir + dataset + "/best/best_mcc_tree.txt";
////			String doubletFile = dir + dataset + "/best/best_posterior_doublets.txt";
//			
//			// 15% missing data
////			String SCFTree = dir + dataset + "/15p_missing_samples/best/best_mcc_tree.txt";
////			String doubletFile = dir + dataset + "/15p_missing_samples/best/best_posterior_doublets.txt";
//			
//			// 30% missing data
//			String SCFTree = dir + dataset + "/30p_missing_samples/best/best_mcc_tree.txt";
//			String doubletFile = dir + dataset + "/30p_missing_samples/best/best_posterior_doublets.txt";
//			
////			double dist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueFile, SCFTree, nCell, BUF);
//			double dist = ClonalTreeDistance.getPairwiseCellSPTreeDistSinglets(trueFile, SCFTree, nCell, doubletFile, BUF);
//			System.out.println(dist);
//			dists[i-1] = Double.toString(dist);
//		}
//		System.out.println(" ");
//		String SCFDist = String.join(" ", dists);
//		System.out.println(SCFDist);
		
//		System.out.println(" ");
//		System.out.println("Root at sc0");
//		System.out.println("OncoNEM Results");
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//			
//			// No Missing Data
//			String oncoNEMCloneTreeFile = dir + "OncoNEM_rslts/" + dataset + "/oncoNEM_clonal_tree_noisy_genotype_mat_" + dataset + ".txt";
//			String oncoCellCloneFile = dir + "OncoNEM_rslts/" + dataset + "/oncoNEM_cell_clones_" + dataset + ".txt";
//			double oncoDist = ClonalTreeDistance.getOncoNEMTreePairwiseCellSPDist(trueFile, oncoNEMCloneTreeFile, oncoCellCloneFile, nCell, BUF);
//			System.out.println(oncoDist);
//			dists[i-1] = Double.toString(oncoDist);
//		}
//		System.out.println(" ");
//		String OncoNEMDist = String.join(" ", dists);
//		System.out.println(OncoNEMDist);
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
//		String trueFile = dir + "dataset1" + "/Orig_tree_" + "dataset1" + ".txt";
//		String SCFTree = dir + "dataset1" + "/best/best_MAP_tree.txt";
//		String doubletFile = dir + "dataset1" + "/best/best_posterior_doublets.txt";
//		double dist = ClonalTreeDistance.getPairwiseCellSPTreeDistSinglets(trueFile, SCFTree, 100, doubletFile, BUF);
//		System.out.println(dist);
		
		/*
		 * No Missing Data SCF results, doublets removed
		 */
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
////			String trueFile = dir + dataset + "/t.txt";
//			String SCFTree = dir + dataset + "/best/best_mcc_tree_r.txt";
//			String doubletFile = dir + dataset + "/best/best_posterior_doublets.txt";
//			
////			String SCFTree = dir + dataset + "/15p_missing_samples/best/best_MAP_tree.txt";
////			String doubletFile = dir + dataset + "/15p_missing_samples/best/best_posterior_doublets.txt";
//			
////			String SCFTree = dir + dataset + "/30p_missing_samples/best/best_MAP_tree.txt";
////			String doubletFile = dir + dataset + "/30p_missing_samples/best/best_posterior_doublets.txt";
//			double dist = ClonalTreeDistance.getPairwiseCellSPTreeDistSinglets(trueFile, SCFTree, 100, doubletFile, BUF);
//			System.out.println(dist);
//		}
		
		/*
		 * Bad SCG results
		 */
//		int i = 4;
//		String dataset = "dataset" + i;
//		System.out.println(dataset);
//		String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//		String SCFTree = dir + dataset + "/best/d_tree.txt";
////		String SCFTree = dir + dataset + "/15p_missing_samples/best/best_MAP_tree.txt";
//		double dist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueFile, SCFTree, 100, BUF);
//		System.out.println(dist);
//		String scgCellCloneFile = dir + dataset + "/30p_SCG_results/scg_cell_clusters.txt";
//		double scgDist = ClonalTreeDistance.getSCGTreePairwiseCellSPDist(trueFile, scgCellCloneFile, 500, BUF);
//		System.out.println(scgDist);
		
		/*
		 * No Missing Data SCG results 
		 */
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//			String scgTreeFile = dir + dataset + "/SCG_results/scg_clone_tree.txt";
////			String scgCellTreeFile = dir + dataset + "/SCG_results/scg_cell_tree.txt";
//			String scgCellCloneFile = dir + dataset + "/SCG_results/scg_cell_clusters.txt";
////			ClonalTreeDistance.writeSCGTree(scgTreeFile, scgCellCloneFile, scgCellTreeFile, BUF);
//			
////			double scgDist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueFile, scgTreeFile, 100, BUF);
////			String scgTreeFile = dir + dataset + "/30p_SCG_results/scg_clone_tree.txt";
////			String scgCellCloneFile = dir + dataset + "/30p_SCG_results/scg_cell_clusters.txt";
//			double scgDist = ClonalTreeDistance.getSCGTreePairwiseCellSPDist(trueFile, scgTreeFile, scgCellCloneFile, 100, BUF);
//			System.out.println(scgDist);
//		}
//		
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//			String scgCellTreeFile = dir + dataset + "/SCG_results/scg_cell_tree_r.txt";
//			double scgDist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueFile, scgCellTreeFile, 100, BUF);
//			System.out.println(scgDist);
//		}
		/*
		 * SCG results, doublets
		 */
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//			String scgTreeFile = dir + dataset + "/SCG_results/scg_clone_tree.txt";
//			String scgCellCloneFile = dir + dataset + "/SCG_results/scg_cell_clusters.txt";
//			String doubletFile = dir + dataset + "/SCG_results/scg_inferred_doubletFlags.txt";
//			
////			String scgTreeFile = dir + dataset + "/15p_SCG_results/scg_clone_tree.txt";
////			String scgCellCloneFile = dir + dataset + "/15p_SCG_results/scg_cell_clusters.txt";
////			String doubletFile = dir + dataset + "/15p_SCG_results/scg_inferred_doubletFlags.txt";
//			
////			String scgTreeFile = dir + dataset + "/30p_SCG_results/scg_clone_tree.txt";
////			String scgCellCloneFile = dir + dataset + "/30p_SCG_results/scg_cell_clusters.txt";
////			String doubletFile = dir + dataset + "/30p_SCG_results/scg_inferred_doubletFlags.txt";
//			double scgDist = ClonalTreeDistance.getSCGTreePairwiseCellSPDistSinglets(trueFile, scgTreeFile, scgCellCloneFile, 500, doubletFile, BUF);
//			System.out.println(scgDist);
//		}
		
		/*
		 * 15% Missing Data SCG results 
		 */
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
////			String scgTreeFile = dir + dataset + "/15p_SCG_results/scg_clone_tree.txt";
////			String scgCellCloneFile = dir + dataset + "/15p_SCG_results/scg_cell_clusters.txt";
//			
//			String scgTreeFile = dir + dataset + "/30p_SCG_results/scg_clone_tree.txt";
//			String scgCellCloneFile = dir + dataset + "/30p_SCG_results/scg_cell_clusters.txt";
//			double scgDist = ClonalTreeDistance.getSCGTreePairwiseCellSPDist(trueFile, scgTreeFile, scgCellCloneFile, 500, BUF);
//			System.out.println(scgDist);
//		}
		
		/*
		 * No Missing Data SCF results
		 */
		
//		String dataset = "dataset" + 3;
//		String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//		STITreeBipartition sb = ClonalTreeDistance.getBipartition(trueFile, BUF);
//		System.out.println(sb);
		
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
////			String trueFile = dir + dataset + "/t.txt";
////			String SCFTree = dir + dataset + "/best/best_MAP_tree.txt";
//			String SCFTree = dir + dataset + "/samples/best/best_mcc_tree_r.txt";
////			String SCFTree = dir + dataset + "/30p_missing_samples/best/nu_MAP_tree.txt";
////			String SCFTree = dir + dataset + "/15p_missing_samples/best/best_mcc_tree_r.txt";
//			double dist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueFile, SCFTree, 100, BUF);
//			System.out.println(dist);
//		}
		
		/*
		 * 30% Missing Data SCF results
		 */
//		for (int i = 1; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
////			String SCFTree = dir + dataset + "/samples/best/best_MAP_tree.txt";
//			String SCFTree = dir + dataset + "/30p_missing_samples/best_MAP_tree.txt";
//			double dist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueFile, SCFTree, 500, BUF);
//			System.out.println(dist);
//		}
		
		
		/*
		 * OncoNEM results 
		 */
//		for (int i = 5; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//			
//			// No Missing Data
//			String oncoNEMCloneTreeFile = dir + "OncoNEM_rslts/" + dataset + "/oncoNEM_clonal_tree_noisy_genotype_mat_" + dataset + ".txt";
//			String oncoCellCloneFile = dir + "OncoNEM_rslts/" + dataset + "/oncoNEM_cell_clones_" + dataset + ".txt";
//			
//			String oncoCellTreeFile = dir + "OncoNEM_rslts/" + dataset + "/oncoNEM_cell_tree.txt";
//			
//			ClonalTreeDistance.writeOncoNEMTree(oncoNEMCloneTreeFile, oncoCellCloneFile, oncoCellTreeFile, BUF);
////			double dist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueFile, oncoCellTreeFile, 100, BUF);
////			System.out.println(dist);
//			// 15% Missing data
////			String oncoNEMCloneTreeFile = dir + "OncoNEM_rslts/" + dataset + "/oncoNEM_15p_missing_clonal_tree_noisy_genotype_mat_" + dataset + ".txt";
////			String oncoCellCloneFile = dir + "OncoNEM_rslts/" + dataset + "/oncoNEM_15p_missing_cell_clones_" + dataset + ".txt";
//						
////			double oncoDist = ClonalTreeDistance.getOncoNEMTreePairwiseCellSPDist(trueFile, oncoNEMCloneTreeFile, oncoCellCloneFile, 100, BUF);
////			System.out.println(oncoDist);
//		}
		
//		for (int i = 5; i <= 10; i++){
//			String dataset = "dataset" + i;
//			System.out.println(dataset);
//			String trueFile = dir + dataset + "/Orig_tree_" + dataset + ".txt";
//			String oncoCellTreeFile = dir + "OncoNEM_rslts/" + dataset + "/oncoNEM_cell_tree_r.txt";
//			double dist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueFile, oncoCellTreeFile, 100, BUF);
//			System.out.println(dist);
//		}
		
		
//		ttr.getRoot().setName("R");
//		// Name each internal node
//		for (STINode<Double> n : ttr.getNodes()){
//			if (n.getName().equals("")){
//				n.setName("IN"+Integer.toString(n.getID()));
//			}			
//		}
//		for (STINode<Double> n : ttr.getNodes()){
//			if (!n.isLeaf())
//				n.setParentDistance(1);
//			else 
//				n.setParentDistance(0);
//		}
//		ttr.getRoot().setParentDistance(0);
//		
//		System.out.println(ttr.toNewick());
//		
//		STINode<Double> test = ClonalTreeDistance.findMRCA(ttr, ttr.getNode("1"), ttr.getNode("2"));
//		System.out.println(test.getName());
//		
//		double dist = ClonalTreeDistance.findPairwiseDist(ttr, ttr.getNode("5"), ttr.getNode("2"));
//		System.out.println(dist);
	}

}
