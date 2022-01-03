package siCloneFiT.objects;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;

import cc.mallet.types.DenseMatrix;
import java.util.Arrays;
//import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeBipartition;
import siCloneFiT.utils.SCFUtilityFunctions;

public class HasseDiGraph {
	
	protected HashMap<String, HasseGraphNode> _name2Node = new HashMap<String, HasseGraphNode>();
	protected HasseGraphNode root;
	protected Hashtable<String, HasseGraphNode> _name2NodeInDegree = new Hashtable<>();
	protected Hashtable<STITreeCluster, HasseGraphNode> cluster2Node = new Hashtable<>();
	protected ArrayList<HasseGraphNode> leaves = new ArrayList<>();
	protected ArrayList<HasseGraphNode> clusterNodes = new ArrayList<>();
	protected String[] leafSet;
	
	/**
	 * Constructor, constructs Hasse graph from a tree
	 * @param tree
	 */
	public HasseDiGraph (STITree<Double> tree){
		// Create the nodes for the leaves of the tree
		this.leafSet = tree.getLeaves();
		for (String leaf: leafSet){
//			leafSet.add(leaf);
			HasseGraphNode node = new HasseGraphNode(leaf);	
			node.nodeLeaves.add(leaf);
			_name2Node.put(leaf, node);
			leaves.add(node);
		}
		Arrays.sort(leafSet);
		
		List<STITreeCluster> clusters = tree.getClusters(tree.getLeaves(), false);
		// Clusters obtained from tree are created
		for (STITreeCluster s: clusters){
			String[] leavesS = s.getClusterLeaves();
			String name = String.join("_", leavesS);
			HasseGraphNode node = new HasseGraphNode(name);	
			for (String leaf: s.getClusterLeaves())
				node.nodeLeaves.add(leaf);
			_name2Node.put(name, node);
			cluster2Node.put(s, node);
			clusterNodes.add(node);
		}
		// Partial order between the clusters are established
		for (STITreeCluster s: clusters){
			if (getNumSubClusters(s, clusters) == 1){
				String[] leavesS = s.getClusterLeaves();
				for (String leaf: leavesS){
					HasseGraphNode node = _name2Node.get(leaf);
					node.outDegree++;
					node.adjacentNextNodes.add(cluster2Node.get(s).name);
				}
			}
			else{
				for (STITreeCluster sc: getSubClusters(s, clusters)){
					HasseGraphNode node = cluster2Node.get(sc);
					node.outDegree++;
					node.adjacentNextNodes.add(cluster2Node.get(s).name);
				}
			}
		}
	}
	
	/**
	 * Return number of subclusters
	 * @param s
	 * @param clusters
	 * @return
	 * Created On: Feb 8, 2018
	 */
	public int getNumSubClusters(STITreeCluster s, List<STITreeCluster> clusters){
		int num = 0;
		for (STITreeCluster c: clusters){
			if (s.containsCluster(c))
				num++;
		}
		return num;
	}
	
	/**
	 * Get the list of subclusters
	 * @param s
	 * @param clusters
	 * @return
	 * Created On: Feb 8, 2018
	 */
	public ArrayList<STITreeCluster> getSubClusters(STITreeCluster s, List<STITreeCluster> clusters){
		ArrayList<STITreeCluster> subClusters = new ArrayList<>();
		for (STITreeCluster c: clusters){
			if (s.containsCluster(c))
				if (s != c)
					subClusters.add(c);
		}
		return subClusters;
	}
	
	public void augmentHasseDiGraph(String newTreeS, SCFUtilityFunctions SCF){
		STITree<Double> newTree = SCF.getTree(newTreeS);
		String[] newTreeLeaves = newTree.getLeaves();
		if (! compareArrays(leafSet, newTreeLeaves))
			throw new IllegalArgumentException("The tree has leaf-set different from the leaf-set of the Hasse Graph");
		// Increase the count of the leaves
		for (String leaf: newTreeLeaves){
			this._name2Node.get(leaf).incrCount();
		}
		
		List<STITreeCluster> newTreeClusters = newTree.getClusters(newTreeLeaves, false);
		
		// Add the clusters
		for (STITreeCluster sc: newTreeClusters){
			HasseGraphNode nodeCluster = checkNodeExistence(sc);
			// Cluster exists, increment occurrence
			if (nodeCluster != null)
				nodeCluster.incrCount();
			// Cluster does not exist, add a new node in graph
			else{
				String[] scLeaves = sc.getClusterLeaves();
				String name = String.join("_", scLeaves);
				HasseGraphNode node = new HasseGraphNode(name);
				for (String leaf: scLeaves)
					node.nodeLeaves.add(leaf);
				_name2Node.put(name, node);
				cluster2Node.put(sc, node);
				clusterNodes.add(node);
			}
		}
		
		// Partial order between the clusters are established
		for (STITreeCluster s: newTreeClusters){
			String[] leavesS = s.getClusterLeaves();
			String name = String.join("_", leavesS);
			HasseGraphNode sNode = checkNodeExistence(s);
//			System.out.println(s);
			if (getNumSubClusters(s, newTreeClusters) == 1){
				
				
				for (String leaf: leavesS){
//					System.out.println(leaf);
					HasseGraphNode node = _name2Node.get(leaf);
//					System.out.println(node.adjacentNextNodes);
//					System.out.println(cluster2Node.get(s));
//					boolean bb = node.adjacentNextNodes.contains(cluster2Node.get(s).name);
//					System.out.println(bb);
					
					if (! node.adjacentNextNodes.contains(sNode.name)){
						node.outDegree++;
						node.adjacentNextNodes.add(sNode.name);
					}
				}
			}
			else{
				for (STITreeCluster sc: getSubClusters(s, newTreeClusters)){
					String scName = String.join("_", sc.getClusterLeaves());
					HasseGraphNode node = checkNodeExistence(sc);
//					System.out.println("sc= " + sc);
//					System.out.println(node.name);
//					System.out.println(node.adjacentNextNodes);
//					System.out.println("here = "+sNode.name);
					if (! node.adjacentNextNodes.contains(sNode.name)){
						node.outDegree++;
						node.adjacentNextNodes.add(sNode.name);
					}
				}
			}
		}
		
			
	}
	
	/**
	 * Check if a cluster already exists in the graph
	 * @param sc
	 * @return
	 * Created On: Feb 8, 2018
	 */
	public HasseGraphNode checkNodeExistence(STITreeCluster sc){
		String[] scLeaves = sc.getClusterLeaves();
		for (HasseGraphNode cluster: this.clusterNodes){
			if (compareArrNodeLeaves(scLeaves, cluster)){
				return cluster;
			}
		}
		return null;
	}
	
	/**
	 * Check if the the array of leaves is same as that of the node
	 * @param scLeaves
	 * @param node
	 * @return
	 * Created On: Feb 8, 2018
	 */
	public boolean compareArrNodeLeaves(String[] scLeaves, HasseGraphNode node){
		String[] nodeLeaves = node.nodeLeaves.toArray(new String[node.nodeLeaves.size()]);
		Arrays.sort(nodeLeaves);
		Arrays.sort(scLeaves);
		return Arrays.equals(nodeLeaves, scLeaves);
	}
	
	/**
	 * Compare if two arrays are equal
	 * @param arr1
	 * @param arr2
	 * @return
	 * Created On: Feb 8, 2018
	 */
	public boolean compareArrays(String[] arr1, String[] arr2){
		Arrays.sort(arr2);
		return Arrays.equals(arr1, arr2);
	}
	
	/**
	 * Print the adjacency matrix
	 * 
	 * Created On: Feb 8, 2018
	 */
	public void printAdjMatrix(){
		ArrayList<HasseGraphNode> allNodes = new ArrayList<>();
		for (HasseGraphNode node: this.leaves){
			allNodes.add(node);
		}
		for (HasseGraphNode node: this.clusterNodes){
			allNodes.add(node);
		}
		String nodeNames = "";
		for (HasseGraphNode node: allNodes)
			nodeNames += " " + node.name;
		System.out.println(nodeNames);
		
		String[][] adjMat = new String[allNodes.size()][allNodes.size()] ;
		for (int i = 0; i< allNodes.size(); i++){
			HasseGraphNode iNode = allNodes.get(i);
			for (int j = 0; j < allNodes.size(); j++){
				String jNodeName = allNodes.get(j).name;
				if (iNode.adjacentNextNodes.contains(jNodeName))
					adjMat[i][j] = "TRUE";
				else
					adjMat[i][j] = "FALSE";
				
			}
		}
		for (String[] s : adjMat){
			System.out.println(Arrays.toString(s));
		}
	}

	public static void main(String[] args) {
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
		// TODO Auto-generated method stub
		String t = "((1,2,3,4),((5,6),(7,8,9,10)));";
		STITree<Double> tree = SCF.getTree(t);
		List<STITreeCluster> clusters = tree.getClusters(tree.getLeaves(), false);
		
		HasseDiGraph HD = new HasseDiGraph(tree);

		String t2 = "(((1,2),(3,4)),(((5,6,7),8),(9,10)));";
		HD.augmentHasseDiGraph(t2, SCF);
		
		for (HasseGraphNode n : HD._name2Node.values()){
			System.out.println(n.name);
//			System.out.println(n.outDegree);
			System.out.println(n.countOccurrence);
		}
		
		HD.printAdjMatrix();
//		for (STITreeCluster s: clusters){
//			System.out.println(s);
//			System.out.println(HasseDiGraph.getNumSubClusters(s, clusters));			
//		}
	}

}
