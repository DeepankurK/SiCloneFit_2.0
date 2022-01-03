/**
 * Aug 10, 2017
 */
package siCloneFiT.proposal;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.distribution.ExponentialDistribution;

import cern.jet.random.tdouble.Normal;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import cern.jet.random.tdouble.engine.DoubleRandomEngine;
import cern.jet.stat.tdouble.Gamma;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import siCloneFiT.objects.Clone;
import siCloneFiT.utils.MultifurcatedPrior;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Aug 10, 2017
 */
public class TreeProposal {
	
	private Random _rng;
	protected double _w1;
	protected double _w2;
	protected double[] _targetBranchPartitionLengths = new double[2];
	protected double _nuNodeBranchLength;
	protected double _sigma = 0.5;
	protected double _jacobian;
	protected double _branchLengthPriorRatio;
	public ExponentialDistribution blPriorDist = new ExponentialDistribution(1);
	public SCFUtilityFunctions SCF = new SCFUtilityFunctions();
	
//	private Random _rng;
	private boolean _topologyHasChanged = false;
	private double _extension = 0.5;
	private double _tuning = 2.0 * Math.log(1.6);
	private double _tuningNNI = 2.0 * Math.log(1.20);
	private double _local = 2.0 * Math.log(1.05);
	private double _sigmaN = 0.03;
	private double _proposalRatio = 0.0;
	public double _maxV = 20.0;

	
	public TreeProposal(){
		this._rng = new Random();
		
	}
	
	/**
	 * Add a new Clone to the tree, adds to the branch connected to oldClone
	 * Adds the new clone as a sibling of the old clone
	 * @param intree
	 * @param nuClone
	 * @param oldClone
	 * @return
	 * Created On: Oct 3, 2017
	 */
	public STITree<Clone> addNode(STITree<Clone> intree, Clone nuClone, Clone oldClone){
		//Create copy of the original tree.
		STITree<Clone> nuTree = new STITree<>(intree);
		STINode<Clone> u = nuTree.getNode(oldClone.cloneName);
		double targetBranchLength = u.getParentDistance();
		// Create the intermediate and new node to add in the tree
		STINode<Clone> uParent = u.getParent();
		STINode<Clone> intermediateNode = uParent.createChild();
		STINode<Clone> nuCloneNode = intermediateNode.createChild();
		// Populate the new clone node
		nuCloneNode.setData(nuClone);
		nuCloneNode.setName(nuClone.cloneName);

		// Intermediate node adopts u as child
		intermediateNode.adoptChild(u);

		// Compute the branch Lengths of the new branches
		_w1 = _rng.nextDouble();		
		_targetBranchPartitionLengths[0] = targetBranchLength * _w1;
		_targetBranchPartitionLengths[1] = targetBranchLength * (1 - _w1);

		_w2 = _rng.nextDouble();		
		
		_nuNodeBranchLength = -Math.log(1 - _w2);  // New formulation for branch length, -ln(1-w2)/delta
		
		_jacobian = targetBranchLength/(1 - _w2); // new jacobian

		//		System.out.println("_jacobian = "+ _jacobian);

		// Calculation of the Branch Length Prior Ratio
		double newTreeBLPrior = blPriorDist.density(_targetBranchPartitionLengths[0])*blPriorDist.density(_targetBranchPartitionLengths[1])*blPriorDist.density(_nuNodeBranchLength);
		double oldTreeBLPrior = blPriorDist.density(targetBranchLength);
		this._branchLengthPriorRatio = newTreeBLPrior/oldTreeBLPrior;
		//		System.out.println("_branchLengthPriorRatio = "+_branchLengthPriorRatio);

		// Set the branch lengths of the new branches
		u.setParentDistance(_targetBranchPartitionLengths[1]);
		intermediateNode.setParentDistance(_targetBranchPartitionLengths[0]);
		nuCloneNode.setParentDistance(_nuNodeBranchLength);
		
		String nuTreeNewick = nuTree.toNewick();
		return SCF.getCloneTree(nuTreeNewick);
	}
	
	/**
	 * Add a new Clone to the tree, adds to a random branch
	 * @param intree, input tree
	 * @param nuClone, nu clone to add as a leaf
	 * @return The new tree
	 * Created On: Aug 10, 2017
	 */
	public STITree<Clone> addNode(STITree<Clone> intree, Clone nuClone){
		//Create copy of the original tree.
		STITree<Clone> nuTree = new STITree<>(intree);
		
		if (intree.getLeafCount() > 1){
			//		STINode<Clone> u = nuTree.getNode(oldClone.cloneName);
			// Select a random node other than root
			STINode<Clone> u = nuTree.selectRandomNode(true, false);
			//		STINode<Clone> u = selectRandomNode(nuTree, true, false, _rng);
			double targetBranchLength = u.getParentDistance();


			// Create the intermediate and new node to add in the tree
			STINode<Clone> uParent = u.getParent();
			STINode<Clone> intermediateNode = uParent.createChild();
			STINode<Clone> nuCloneNode = intermediateNode.createChild();

			// Populate the new clone node
			nuCloneNode.setData(nuClone);
			nuCloneNode.setName(nuClone.cloneName);

			// Intermediate node adopts u as child
			intermediateNode.adoptChild(u);

			// Compute the branch Lengths of the new branches
			_w1 = _rng.nextDouble();		
			_targetBranchPartitionLengths[0] = targetBranchLength * _w1;
			_targetBranchPartitionLengths[1] = targetBranchLength * (1 - _w1);

			_w2 = _rng.nextDouble();		
//			_nuNodeBranchLength = targetBranchLength * Math.exp(-_sigma*(_w2 - 0.5)); // old formulation for branch length
					_nuNodeBranchLength = -Math.log(1 - _w2);  // New formulation for branch length, -ln(1-w2)/delta
//			_jacobian = _sigma * _nuNodeBranchLength * targetBranchLength; // old jacobian
					_jacobian = targetBranchLength/(1 - _w2); // new jacobian

			//		System.out.println("_jacobian = "+ _jacobian);

			// Calculation of the Branch Length Prior Ratio
			double newTreeBLPrior = blPriorDist.density(_targetBranchPartitionLengths[0])*blPriorDist.density(_targetBranchPartitionLengths[1])*blPriorDist.density(_nuNodeBranchLength);
			double oldTreeBLPrior = blPriorDist.density(targetBranchLength);
			this._branchLengthPriorRatio = newTreeBLPrior/oldTreeBLPrior;
			//		System.out.println("_branchLengthPriorRatio = "+_branchLengthPriorRatio);

			// Set the branch lengths of the new branches
			u.setParentDistance(_targetBranchPartitionLengths[1]);
			intermediateNode.setParentDistance(_targetBranchPartitionLengths[0]);
			nuCloneNode.setParentDistance(_nuNodeBranchLength);
		}
		// Current tree contains only one node, the new clone 
		// will be a sibling of this node, new root to create
		else{
			STINode<Clone> u = nuTree.getNode(nuTree.getLeaves()[0]);
			double targetBranchLength = u.getParentDistance();
//			System.out.println(u.getName());
//			System.out.println(u.getParentDistance());
			
			STINode<Clone> newRoot = u.createChild("R");
			
			STINode<Clone> nuCloneNode = newRoot.createChild();
			STINode<Clone> uCopy = newRoot.createChild();
			uCopy.setParentDistance(targetBranchLength);

			// Populate the new clone node
			nuCloneNode.setData(nuClone);
			nuCloneNode.setName(nuClone.cloneName);
			
			// Root the tree at new Root
			nuTree.rerootTreeAtNode(newRoot);
						
			// Remove u, name new copy of u
			nuTree.removeNode(u.getName());
			uCopy.setName(u.getName());
			
			_w2 = _rng.nextDouble();		
//			_nuNodeBranchLength = targetBranchLength * Math.exp(-_sigma*(_w2 - 0.5));
			_nuNodeBranchLength = -Math.log(1 - _w2);  // New formulation for branch length, -ln(1-w2)/delta
			nuCloneNode.setParentDistance(_nuNodeBranchLength);
//			_jacobian = _sigma * _nuNodeBranchLength * targetBranchLength; // old jacobian
			_jacobian = targetBranchLength/(1 - _w2); // new jacobian
			this._branchLengthPriorRatio = blPriorDist.density(_nuNodeBranchLength);
			
		}
		String nuTreeNewick = nuTree.toNewick();
		return SCF.getCloneTree(nuTreeNewick);
		
	}
	
	/**
	 * Remove a clone from the tree
	 * @param intree
	 * @param oldClone
	 * @return
	 * Created On: Aug 10, 2017
	 */
	public STITree<Clone> removeNode(STITree<Clone> intree, Clone oldClone){
		//Create copy of the original tree.
		STITree<Clone> nuTree = new STITree<>(intree);
		
		// Get the leaf that is to be removed from the tree
		STINode<Clone> cloneToRemove = nuTree.getNode(oldClone.cloneName);
		STINode<Clone> cloneToRemoveParent = cloneToRemove.getParent();
		
		// Obtain the other child of the intermediate node, assumes binary tree
		ArrayList<STINode<Clone>> children = new ArrayList<>();
		for (STINode<Clone> c : cloneToRemoveParent.getChildren()){
			if (c.getID() != cloneToRemove.getID())
				children.add(c);
		}
		
		STINode<Clone> otherChild = children.get(0);
		
//		if (intree.getLeafCount() == 2)
//			System.out.println("intree = " + cloneToRemoveParent.isRoot());
//		System.out.println(otherChild.getID());
		
		if (cloneToRemoveParent.isRoot() == false){
		
			// Obtain the branch lengths for the branches to be affected
			double bl1 = cloneToRemoveParent.getParentDistance();
			double bl2 = otherChild.getParentDistance();
			double bl3 = cloneToRemove.getParentDistance();

			// Calculate Jacobian
//			_jacobian = 1/(_sigma * bl3 * (bl1 + bl2)); // old jacobian
			
			_jacobian = Math.exp(-bl3)/(bl1 + bl2);  // new jacobian

			// Calculation of the Branch Length Prior Ratio
			double newTreeBLPrior = blPriorDist.density(bl1 + bl2);
			double oldTreeBLPrior = blPriorDist.density(bl1) * blPriorDist.density(bl2) * blPriorDist.density(bl3);
			this._branchLengthPriorRatio = newTreeBLPrior/oldTreeBLPrior;

			// Node removal and adoption
			
			cloneToRemoveParent.getParent().adoptChild(otherChild);
			otherChild.setParentDistance(bl1 + bl2);
			cloneToRemoveParent.removeChild(cloneToRemove, true);
			cloneToRemoveParent.getParent().removeChild(cloneToRemoveParent, true);
		}
		else{
//			System.out.println("root has ID = "+ cloneToRemoveParent.getID());
			// Obtain the branch lengths for the branches to be affected
			
			double bl2 = otherChild.getParentDistance();
			double bl3 = cloneToRemove.getParentDistance();
			// Calculate Jacobian
//			_jacobian = 1/(_sigma * bl3 * (bl2));
			
			_jacobian = Math.exp(-bl3)/(bl2);

			// Calculation of the Branch Length Prior Ratio
			double newTreeBLPrior = 1;
			double oldTreeBLPrior = blPriorDist.density(bl2) * blPriorDist.density(bl3);
			this._branchLengthPriorRatio = newTreeBLPrior/oldTreeBLPrior;
			
			// Node removal and adoption
			cloneToRemoveParent.removeChild(cloneToRemove, true);
			
			nuTree.rerootTreeAtNode(otherChild);
			
			// The current tree has only two leaves.
			// After removing the clone, the other leaf will be the tree
			if (intree.getLeafCount() == 2){
				
//				System.out.println(otherChild.getID());
//				System.out.println(nuTree.toNewick());
//				System.out.println(nuTree.getRoot().getID());
				return SCF.getCloneTree(otherChild.toNewick());
			}
//			for (STINode<Clone> n: otherChild.getChildren())
//				System.out.println(n.getID());
//			nuTree.rerootTreeAtEdge(otherChild);
			otherChild.removeChild(cloneToRemoveParent, true);
//			System.out.println("previous root removed");
//			for (STINode<Clone> n: otherChild.getChildren())
//				System.out.println(n.getID());
		}
		String nuTreeNewick = nuTree.toNewick();
		return SCF.getCloneTree(nuTreeNewick);
//		return nuTree;
		
	}
	
	
	/**
	 * Randomly selects a node in the tree 
	 * @param tree
	 * @param include_leaves
	 * @param include_root
	 * @param rng
	 * @return
	 * Created On: Aug 10, 2017
	 */
	private STINode<Clone> selectRandomNode(STITree<Clone> tree, boolean include_leaves, boolean include_root, Random rng) {
		ArrayList<Integer> idList = new ArrayList<Integer> ();
		
		for (int i = 0; i < tree.getNodeCount(); i++) {
			STINode<Clone> curNode = tree.getNode(i);
			if (curNode != null) {
				if (curNode.isLeaf()) {
					if (include_leaves) {
						idList.add(i);
					}
				}
				else if (curNode.isRoot()) {
					if (include_root) {
						idList.add(i);
					}
				}
				else {
					idList.add(i);
				}
			}
			else {
				System.out.println("ERROR: Random node selection: Null Node.");
				return null;
			}
		}
		
		if (idList.size() == 0) {
			System.out.println("ERROR: Random node selection: Empty Set.");
			return null;
		}		
		double randDouble = rng.nextDouble();
		int randIndex = (int) (randDouble * idList.size());
		return tree.getNode(idList.get(randIndex));
	}
	
	/**
	 * Compute the prior ratio for two trees,
	 * product of the topologyRatio and branchLengthRatio
	 * topologyRatio = 
	 * @param nuTreeLeafCount, number of leaves in the new tree
	 * @param moveFlag, 1 if addEdge move generated the new tree, 0 else
	 * @return
	 * Created On: Aug 14, 2017
	 */
	public double getTreePriorRatio(int nuTreeLeafCount, int moveFlag){
		double topologyPriorRatio;
		if (moveFlag == 1)
			topologyPriorRatio = (double)(nuTreeLeafCount - 1)/((nuTreeLeafCount - 2)*(2*nuTreeLeafCount -3));
		else
			topologyPriorRatio = (double)((nuTreeLeafCount - 1)*(2*nuTreeLeafCount -1))/(nuTreeLeafCount);
//		System.out.println("tpr="+Math.log(topologyPriorRatio));
		return topologyPriorRatio*_branchLengthPriorRatio;
	}
	
	/**
	 * Compute tree prior ratio for multifurcated topology
	 * @param nuTreeLeafCount
	 * @param nCell
	 * @param moveFlag
	 * @return
	 * Created On: Feb 8, 2018
	 */
	public double getTreePriorRatio(int nuTreeLeafCount, int nCell, int moveFlag){
		double topologyPriorRatio;
		double currTreeProb = MultifurcatedPrior.computeNumberTrees(nCell,nuTreeLeafCount);
		double oldTreeProb = MultifurcatedPrior.computeNumberTrees(nCell,nuTreeLeafCount-1);
		if (moveFlag == 1)
			topologyPriorRatio = oldTreeProb/currTreeProb;
		else
			topologyPriorRatio = currTreeProb/oldTreeProb;
		return topologyPriorRatio*_branchLengthPriorRatio;
	}
	
	/**
	 * Compute prior ratio for the number of clones
	 * @param nuListClone
	 * @param oldListClone
	 * @param alpha_0
	 * @return
	 * Created On: Sep 19, 2017
	 */
	public double getClusteringLogPriorRatio(ArrayList<Clone> nuListClone, ArrayList<Clone> oldListClone, double alpha_0){
		return getClusteringLogPriorProb(nuListClone, alpha_0) - getClusteringLogPriorProb(oldListClone, alpha_0);
	}
	
	public double getClusteringLogPriorProb(ArrayList<Clone> listClone, double alpha_0){
		int nClone = listClone.size();
		double logProb = nClone * Math.log(alpha_0);
		for (Clone C: listClone){
			logProb += Gamma.logGamma(C.memberCellList.size());
		}
		return logProb;
	}
	
	/**
	 * Compute Hastings ratio for adding a new clone
	 * @param alpha_0
	 * @param nCloneNu
	 * @param nCloneOld
	 * @return
	 * Created On: Aug 14, 2017
	 */
	public double getHastingsRatioAddClone(double alpha_0, int nCloneNu, int nCloneOld){
		double nEdge = 2*nCloneOld - 2;
//		return nEdge/(alpha_0 * nCloneNu);
		return (double) nEdge/nCloneNu;
	}
	
	public double getHastingsRatioRemoveClone(int nCloneNu, int nCloneOld){
		double nEdge = 2*nCloneNu - 2;
		return nCloneOld/nEdge;
	}

	public double getJacobian(){
		return this._jacobian;
	}
	
	/**
	 * Bunch of setter functions follow:
	 */
	
	/*
	 * maxBranchLength: Maximum allowed branch lengths for consideration.
	 */
	public void setMaxBranchLength(double maxBranchLength) {
		this._maxV = maxBranchLength;
	}
	
	/*
	 * LOCAL tuning: Used in LOCAL(Larget and Simon)
	 */	
	public void setLocalParameter(double local) {
		this._local = local;
	}

	/*
	 * tuningParameter: Used in Multiplier
	 */
	public void setTuningParameter(double tuningParameter) {
		this._tuning = tuningParameter;
	}
	
	/*
	 * sigmaParameter:  Used in CC.
	 */
	public void setSigmaParameter (double sigmaParameter) {
		this._sigmaN = sigmaParameter;
	}
	
	
	/**
	 * Returns proposal ratio for the last perturbation.
	 */
	public double getProposalRatio() {
		return _proposalRatio;
	}
	
	/**
	 * Indicator of whether last perturbation changed topology.
	 */
	public boolean hasTopologyChanged(){
		return _topologyHasChanged;
	}

	
	/**
	 *  Tree Perturbation Code Follows. 
	 *  Implementation of stNNI, SPR, SS.
	 */	

	/*
	 * stNNI: stochastic NNI.
	 */
	public STITree<Clone> stNNI (STITree<Clone> otree) {

		//Create copy of the original tree.
		STITree<Clone> tree = new STITree<Clone> (otree);
		
		STINode<Clone> u = tree.selectRandomNode(false, false);
//		STINode<Clone> u = selectRandomNode(tree, false, false, _rng);
		STINode<Clone> v = u.getParent();
		
		// Get children of "u".
		STINode<Clone> a = getRandomChild(u);
		STINode<Clone> b = getOtherChild(u,a);
		
		// Get other child of "v".
		STINode<Clone> c = getOtherChild(v,u);
		
		double randDouble = _rng.nextDouble();
		//Swap "a" with "c".
		if (randDouble < 1.0/3) {
			v.adoptChild(a);
			u.adoptChild(c);	
			_proposalRatio = 1;
			_topologyHasChanged = true;
		}
		//Swap "b" with "c".
		else if (randDouble < 2.0/3){
			v.adoptChild(b);
			u.adoptChild(c);
			_proposalRatio = 1;
			_topologyHasChanged = true;
		}
		else {
			_topologyHasChanged = false;
//			ArrayList<STINode<Clone>> nodes = new ArrayList<STINode<Clone>>();
//			nodes.add(u);
//			nodes.add(a);
//			nodes.add(b);
//			nodes.add(c);
//			tree = multiplier (tree, nodes, _tuningNNI);
		}
//		System.out.println(tree.toNewick());
//		System.out.println(_topologyHasChanged);
		ArrayList<STINode<Clone>> nodes = new ArrayList<STINode<Clone>>();
		nodes.add(u);
		nodes.add(a);
		nodes.add(b);
		nodes.add(c);
		tree = multiplier (tree, nodes, _tuningNNI);

		return tree;
	}
	
	private STINode<Clone> getOtherChild(STINode<Clone> node, STINode<Clone> c) {
		assert node.getChildCount() <= 3;		
		ArrayList<STINode<Clone>> a = new ArrayList<STINode<Clone>>();
		
		for (STINode<Clone> child : node.getChildren()) {
			if (child.getID() != c.getID()) {
				a.add(child);
			}
		}			
		
		int randIndex = (int) (_rng.nextDouble() * a.size());
		return a.get(randIndex);
	}

	private STINode<Clone> getRandomChild(STINode<Clone> node) {
		assert node.getChildCount() >= 2;
		ArrayList<STINode<Clone>> a = new ArrayList<STINode<Clone>>();
		for (STINode<Clone> child : node.getChildren()) {
			a.add(child);
		}

		int randIndex = (int) (_rng.nextDouble() * a.size());
		return a.get(randIndex);
	}

	/*
	 * Changes a list of branch lengths.
	 */
	public STITree<Clone> multiplier (STITree<Clone> otree) {
		//Create copy of the original tree.
		STITree<Clone> tree = new STITree<Clone> (otree);
		
		STINode<Clone> u = tree.selectRandomNode(false, false);
//		STINode<Clone> u = selectRandomNode(tree, false, false, _rng);
		STINode<Clone> v = u.getParent();	
		
		// Get children of "u".
		STINode<Clone> a = getRandomChild(u);
		STINode<Clone> b = getOtherChild(u,a);

		// Get other child of "v".
		STINode<Clone> c = getOtherChild(v,u);		
		
		ArrayList<STINode<Clone>> nodes = new ArrayList<STINode<Clone>>();
		nodes.add(u);
		nodes.add(a);
		nodes.add(b);
		nodes.add(c);
		tree = multiplier (tree, nodes, _tuning);
		
		return tree;
	}
	
	public STITree<Clone> changeBranchLength(STITree<Clone> otree){
		//Create copy of the original tree.
		STITree<Clone> tree = new STITree<Clone> (otree);
		STINode<Clone> u = tree.selectRandomNode(true, false);
//		STINode<Clone> u = selectRandomNode(tree, true, false, _rng);
		double oldB = u.getParentDistance();
		double newB = 0.0, m = 0.0;
		
		do {
			m = Math.exp(_tuning * (_rng.nextDouble()-0.5));
			newB = oldB * m;
		} while (newB >= _maxV);
		u.setParentDistance(newB);
		return tree;
	}
	
	/**
	 * Get the log prior probability on the branch lengths of a tree
	 * \sum_e log P(e| e_prior)
	 * @param intree
	 * @return
	 * Created On: Sep 18, 2017
	 */
	public double getBranchLengthLogPriorProb(Tree intree){
		double logProb = 0.0;
		for (TNode n: intree.getNodes()){
			if (n.isRoot() == false){
//				System.out.println(n.getParentDistance() + " "+ blPriorDist.density(n.getParentDistance()));
				logProb += Math.log(blPriorDist.density(n.getParentDistance()));
			}
		}
		return logProb;
	}
	
	/**
	 * Compute the log-prior ratio for proposing a new tree
	 * @param nuTree
	 * @param oldTree
	 * @return
	 * Created On: Sep 18, 2017
	 */
	public double getBranchLengthLogPriorRatio(Tree nuTree, Tree oldTree){
		double nuTreeLogPrior = getBranchLengthLogPriorProb(nuTree);
		double oldTreeLogPrior = getBranchLengthLogPriorProb(oldTree);
		return nuTreeLogPrior - oldTreeLogPrior;
	}
	
	/**
	 * Compute the log prior probability of a tree
	 * from uniform distribution
	 * @param tree
	 * @return
	 * Created On: Sep 18, 2017
	 */
	public double getTreeLogPriorProb(Tree tree){
		double topologyProb = getTreeTopologyLogProb(tree.getLeafCount());
		double branchLengthProb = getBranchLengthLogPriorProb(tree);
		return topologyProb + branchLengthProb;
	}
	
	public double getTreeLogPriorProb(Tree tree, int nCell){
		double topologyProb = getTreeTopologyLogProb(tree.getLeafCount(), nCell);
		double branchLengthProb = getBranchLengthLogPriorProb(tree);
		return topologyProb + branchLengthProb;
	}
	
	/**
	 * Compute the uniform prior density of a tree topology
	 * @param nLeaf
	 * @return
	 * Created On: Sep 18, 2017
	 */
	public double getTreeTopologyLogProb(int nLeaf){
//		return (nLeaf - 2)*Math.log(2) + Gamma.logGamma(nLeaf - 2) - Gamma.logGamma(2*nLeaf - 3);
		return (nLeaf - 2) * Math.log(2) + org.apache.commons.math3.special.Gamma.logGamma(nLeaf - 2) - org.apache.commons.math3.special.Gamma.logGamma(2*nLeaf - 3);
	}
	
	public double getTreeTopologyLogProb(int nLeaf, int nCell){
		return 0 - Math.log(MultifurcatedPrior.computeNumberTrees(nCell, nLeaf));
	}
	
	/**
	 * Function that performs "multiplier"
	 */
	private STITree<Clone> multiplier (STITree<Clone> otree, ArrayList<STINode<Clone>> nodes, double tuningMultiplier) {
		//Create copy of the original tree.
		STITree<Clone> tree = new STITree<Clone> (otree);
		
		_proposalRatio = 1.0;
//		_topologyHasChanged = false;
		for (int i = 0; i < nodes.size(); i++) {
			STINode<Clone> curNode = tree.getNode(nodes.get(i).getID());
			double oldB = curNode.getParentDistance();
			double newB = 0.0, m = 0.0;
			do {
				m = Math.exp(tuningMultiplier * (_rng.nextDouble()-0.5));
				newB = oldB * m;
			} while (newB >= _maxV);
			
			curNode.setParentDistance(newB);
			_proposalRatio *= m;
		}

		return tree;
	}
	
	/*
	 * Continuous Change: Jow et al. 2002
	 */
	public STITree<Clone> CC (STITree<Clone> otree) {
		//Create copy of the original tree.
		STITree<Clone> tree = new STITree<Clone> (otree);
	
		STINode<Clone> u = tree.selectRandomNode(false, false);
//		STINode<Clone> u = selectRandomNode(tree, false, false, _rng);
		STINode<Clone> v = getRandomChild(u);
		
		//Normal Distribution: N(0, sigma)
		DoubleRandomEngine rngEngine = new DoubleMersenneTwister(_rng.nextInt());
		Normal rngN = new Normal(0.0,_sigmaN,rngEngine);
		
		double newV = v.getParentDistance() + rngN.nextDouble();
		
		//Branch length should not exceed maximum allowed value.
		if (Math.abs(newV) < _maxV) {
			v.setParentDistance(Math.abs(newV));
		}

		_proposalRatio = 1;
		_topologyHasChanged = false;
		
		//Change Topology
		if (newV < 0 && !v.isLeaf()) {
			_topologyHasChanged = true;
			STINode<Clone> a = getRandomChild(v);
			STINode<Clone> c = getOtherChild(u,v);
			
			v.adoptChild(c);
			u.adoptChild(a);
		}
			
		return tree;
	}
	
	/**
	 * Change topology using Random Subtree Pruning and Regrafting (rSPR) move.
	 * @param otree
	 * @return
	 */
	public STITree<Clone> rSPR (STITree<Clone> otree) {
		//Create copy of the original tree.
		STITree<Clone> tree = new STITree<Clone> (otree);
		
		STINode<Clone> u = tree.selectRandomNode(false, false);
//		STINode<Clone> u = selectRandomNode(tree, false, false, _rng);
		//Node "u" has two children: "v" and "w".
		STINode<Clone> v = getRandomChild(u);
		STINode<Clone> w = getOtherChild(u,v);
		
		//We want to prune subtree rooted at "v".
		//We want to re-graft the pruned subtree between node "x" and it's child "cx".
		STINode<Clone> x = null;
		STINode<Clone> cx = null;
		
		//"x" cannot be any node in the pruned subtree OR the same node from which it was pruned.
		do {
			x = tree.selectRandomNode(false, true);
//			x = selectRandomNode(tree, false, true, _rng);
		} while(v.isAncestor(x) || (u.getID() == x.getID()));

		u.getParent().adoptChild(w);
		x.adoptChild(u);

		//"cx" cannot be: "w" as it gives rise to the same topology again; "u" as it is part of re-grafting.
		do {
			cx = getRandomChild(x);
		} while (cx.getID() == u.getID() || cx.getID() == w.getID());
		u.adoptChild(cx);

		//Use Multiplier:
		ArrayList<STINode<Clone>> nodes = new ArrayList<STINode<Clone>>();
		nodes.add(v);
		nodes.add(cx);
		tree = multiplier (tree, nodes, _tuning);	
		
		_proposalRatio *= 1.0;
		_topologyHasChanged = true;
		return tree;
	}
	
	/**
	 * Change topology using Extending Subtree Pruning and Regrafting (eSPR) move.
	 * @param otree
	 * @return
	 */
	public STITree<Clone> eSPR (STITree<Clone> otree) {
		STITree<Clone> tree = new STITree<Clone> (otree);
		
		//New change: Prevent waste of computation when "w" is leaf.
		STINode<Clone> u, v, w;
		do {
			u = tree.selectRandomNode(false, false);
//			u = selectRandomNode(tree, false, false, _rng);
			//Node "u" has two children: "v" and "w".
			v = getRandomChild(u);
			w = getOtherChild(u,v);			
		} while (w.isLeaf());		
		
		boolean isConstrainedPrunedBranch = v.isLeaf();
		boolean isConstrainedRegraftingBranch = false;
		
		STINode<Clone> x = null;
		STINode<Clone> cx = null;
		
		_proposalRatio = 1;
		_topologyHasChanged = false;
		//When "w" is a leaf, no change possible.
		if (!w.isLeaf()) {
			x = w;
			cx = getRandomChild(w); 
			while (!cx.isLeaf()) {
				if (_rng.nextDouble() < (1 - _extension)) {
					break;
				}
				x = cx;
				cx = getRandomChild(cx);
			}
			
			u.getParent().adoptChild(w);
			x.adoptChild(u);
			u.adoptChild(cx);			
		
			isConstrainedRegraftingBranch = cx.isLeaf();	
			//Set proposal ratio.
			if (isConstrainedPrunedBranch && !isConstrainedRegraftingBranch) {
				_proposalRatio = 1 - _extension;
			}
			else if (!isConstrainedPrunedBranch && isConstrainedRegraftingBranch) {
				_proposalRatio = 1 / (1 - _extension);
			}
			
			
			//Use Multiplier:
			double topologyRatio = _proposalRatio;
			ArrayList<STINode<Clone>> nodes = new ArrayList<STINode<Clone>>();
			nodes.add(v);
			nodes.add(cx);
			tree = multiplier (tree, nodes, _tuning);	
			_proposalRatio *= topologyRatio;
			_topologyHasChanged = true;
			
		}
		
		return tree;
	}	
	
	/* 
	 * Change topology using Random Subtree Swapping (rSTS) move.
	 */
	public STITree<Clone> rSTS (STITree<Clone> otree) {
		//Create copy of the original tree.
		STITree<Clone> tree = new STITree<Clone> (otree);
		
		STINode<Clone> u = tree.selectRandomNode(false, false);
//		STINode<Clone> u = selectRandomNode(tree, false, false, _rng);
		//Node "u" has two children: "v" and "w". As "w" is not involved, it is not defined.
		STINode<Clone> v = getRandomChild(u);
		
		// One of the two subtrees taking part in swap is rooted at "v".
		STINode<Clone> x = null;
		STINode<Clone> cx = null;
		
		//"x" cannot be any node in the pruned subtree OR the same node from which it was pruned.
		do {
			x = tree.selectRandomNode(false, true);
//			x = selectRandomNode(tree, false, true, _rng);
		} while(v.isAncestor(x) || (u.getID() == x.getID()));

		do {
			cx = getRandomChild(x);
		} while (cx.getID() == u.getID() || cx.isAncestor(v));
		
		
		//We want to swap "v" (child of "u") with "cx" (child of "x")
		u.adoptChild(cx);
		x.adoptChild(v);

		
		//Use Multiplier:
		ArrayList<STINode<Clone>> nodes = new ArrayList<STINode<Clone>>();
		nodes.add(v);
		nodes.add(cx);
		tree = multiplier (tree, nodes, _tuning);	
		
		_proposalRatio *= 1.0;
		_topologyHasChanged = true;
		return tree;			
	}

	
	/* 
	 * Change topology using Extending Subtree Swapping (eSTS) move.
	 */
	public STITree<Clone> eSTS (STITree<Clone> otree) {	
		//Create copy of the original tree.
		STITree<Clone> tree = new STITree<Clone> (otree);		
		
		//New change: Prevent waste of computation when "w" is leaf.
		STINode<Clone> u, v, w;
		do {
			u = tree.selectRandomNode(false, false);
//			u = selectRandomNode(tree, false, false, _rng);
			//Node "u" has two children: "v" and "w".
			v = getRandomChild(u);
			w = getOtherChild(u,v);			
		} while (w.isLeaf());
		
		boolean isConstrainedBranchOne = v.isLeaf();
		boolean isConstrainedBranchTwo = false;
		
		STINode<Clone> x = null;
		STINode<Clone> cx = null;
		
		_proposalRatio = 1;
		//When "w" is a leaf, no change possible.
		if (!w.isLeaf()) {
			x = w;
			cx = getRandomChild(w); 
			while (!cx.isLeaf()) {
				if (_rng.nextDouble() < (1 - _extension)) {
					break;
				}
				x = cx;
				cx = getRandomChild(cx);
			}
			
			u.adoptChild(cx);
			x.adoptChild(v);
			
			isConstrainedBranchTwo = cx.isLeaf();	
			//Set proposal ratio.
			if (isConstrainedBranchOne && !isConstrainedBranchTwo) {
				_proposalRatio = 1 - _extension;
			}
			else if (!isConstrainedBranchOne && isConstrainedBranchTwo) {
				_proposalRatio = 1 / (1 - _extension);
			}

			//Use Multiplier:
			double topologyRatio = _proposalRatio;
			ArrayList<STINode<Clone>> nodes = new ArrayList<STINode<Clone>>();
			nodes.add(v);
			nodes.add(cx);
			tree = multiplier (tree, nodes, _tuning);	
			_proposalRatio *= topologyRatio;
			_topologyHasChanged = true;				
		}

		return tree;
	}	
	
	public STITree<Clone> addEdge (STITree<Clone> intree) {
		//Create copy of the original tree.
		STITree<Clone> tree = new STITree<Clone> (intree);
		ArrayList<STINode<Clone>> nonBinaryPars = new ArrayList<>();
		for (STINode<Clone> intNode : tree.getNodes()){
			if (!intNode.isLeaf()){
				if (intNode.getChildCount() > 2)
					nonBinaryPars.add(intNode);
			}
		}
		for (STINode<Clone> node: nonBinaryPars)
			System.out.println(node.getName());
		return null;
		
	}
	
	/**
	 * proposes Tree 
	 * @param curTree
	 * @param nCell
	 * @return
	 */
	public STITree<Clone> proposeTree(STITree<Clone> curTree, int nCell){
		double randDouble = _rng.nextDouble();
		STITree<Clone> newTree = curTree;
		
		//For smaller problems rely on small subset of perturbations:
		if (nCell <= 7) {			
			if (randDouble < 1.0/2) {
				newTree = this.stNNI(curTree);
			}
			else {
				newTree = this.multiplier(curTree);
			}			
		}
		
		//For larger problems rely on mixture of all perturbations:
		else{
			if (randDouble < 1.0/4) {
				newTree = this.multiplier(curTree);
			}
			else if (randDouble < 3.0/5) {
				newTree = this.stNNI(curTree);
			}
			else if (randDouble < 0.7) {
//				System.out.println("eSPR happens");
				newTree = this.eSPR(curTree);
			}
			else if (randDouble < 0.8) {
				newTree = this.rSPR(curTree);
			}
			else if (randDouble < 0.9) {
				newTree = this.eSTS(curTree);
			}
			else {
				newTree = this.rSTS(curTree);
			}			
		}
		return newTree;
	}
	
	public STITree<Clone> proposeTreeParsimony(STITree<Clone> curTree, int nCell){
		double randDouble = _rng.nextDouble();
		STITree<Clone> newTree = curTree;
		if (randDouble < 0.25){
//			System.out.println("stNNI happens");
			newTree = this.stNNI(curTree);
		}
		else if (randDouble < 0.5) {
//			System.out.println("eSPR happens");
			newTree = this.eSPR(curTree);
		}
		else if (randDouble < 0.75) {
			newTree = this.rSPR(curTree);
		}
		else {
			newTree = this.eSTS(curTree);
		}
		return newTree;
	}
	
	/**
	 * @param args
	 * Created On: Aug 10, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
//		String t = "((((C0:1.0,C1:3.0):4.0,C2:1.1):6.0,C3:2.5):2.9,C4:8.0):0.5;";

		String t = "((((C0:1.0,C1:3.0):4.0,C2:1.1):6.0,C3:2.5):2.9,(C6:8.383374399001655,C4:0.45139344027564565):7.548606559724354):0.5;";

//		String t = "((((C0:1.0,C1:3.0):4.0,C2:1.1):6.0,C3:2.5):2.9,C6:8.383374399001655):0.5;";
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
		STITree<Clone> polytree = SCF.getCloneTree(t);
		TreeProposal TP = new TreeProposal();
		Clone nuClone = new Clone(9);
		nuClone.cloneName = "C6";
//		STITree<Clone> nuTree = TP.addNode(polytree, nuClone);
		STITree<Clone> nuTree = TP.removeNode(polytree, nuClone);
		
		
		System.out.println(nuTree.toNewick());
		
		String tt = "((C0:1.0,C1:3.0):4.0,C2:1.1);"; 
		STITree<Clone> tree = SCF.getCloneTree(tt);
		for (STINode<Clone> c : tree.getNodes())
		{if (c.isRoot() == false){
			System.out.println("Node " + c.getID() + "," +c.getName() + " has branch length = " + c.getParentDistance() + " and parent = " + c.getParent().getID());
		}
		else
			System.out.println("Node " + c.getID() + "," +c.getName());
		}
		
		Clone c2 = new Clone(2);
		c2.setNameID(2);
//		System.out.println(c2.cloneName);
		Clone c1 = new Clone(1);
		c1.setNameID(1);
		STITree<Clone> c2Rt = TP.removeNode(tree, c2);
//		tree.rerootTreeAtNode(tree.getNode("C0"));
		System.out.println(c2Rt.toNewick());
		
		STITree<Clone> c1Rt = TP.removeNode(c2Rt, c1);
		System.out.println(c1Rt.toNewick());
		
		
		String ttt = "C0:1.0;";
		STITree<Clone> ttree = SCF.getCloneTree(ttt);
		STITree<Clone> nuRt = TP.addNode(c1Rt, nuClone);
		System.out.println(nuRt.toNewick());
		
//		for (STINode<Clone> c : tree.getNodes())
//		{if (c.isRoot() == false){
//			System.out.println("Node " + c.getID() + "," + c.getName() + " has branch length = " + c.getParentDistance() + " and parent = " + c.getParent().getID());
//		}
//		else
//			System.out.println("Node " + c.getID() + "," +c.getName());
//		}
//		for (STINode<Clone> n: nuTree.getNodes())
//			System.out.println(n.getName());
		STITree<Clone> secondtree = TP.addNode(nuTree, nuClone);
//		System.out.println(secondtree.toNewick());
		
		secondtree.getNode("C0").createChild("C0_1");
		secondtree.getNode("C0").createChild("C0_2");
		
//		System.out.println(TP.getBranchLengthLogPriorProb(nuTree));
		
//		System.out.println(TP._jacobian);
//		System.out.println(TP._branchLengthPriorRatio);
//		System.out.println(secondtree.toNewick());
		
		
//		System.out.println(TP.getHastingsRatioAddClone(1, 5, 4));
//		System.out.println(TP.getTreePriorRatio(5, 1));
//		
//		ExponentialDistribution Ep = new ExponentialDistribution(0.5);
//		System.out.println(Ep.density(0.1));
		
		
//		System.out.println(Gamma.logGamma(4.2));
		System.out.println(TP.getTreePriorRatio(8, 1));
		


	}

}
