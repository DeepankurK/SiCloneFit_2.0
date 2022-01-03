/**
 * Aug 4, 2017
 */
package utils;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.Set;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

import SiFit.utility.BasicUtilityFunctions;
import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.SingleCell;

/**
 * @author hz22
 * Aug 4, 2017
 */
public class SCFUtilityFunctions extends BasicUtilityFunctions {

	/**
	 * 
	 */
	public Random rng = new Random();
	public SCFUtilityFunctions() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param rng
	 */
	/**
	 * Constructs a list of SingleCell Objects, assigns ID (Zero-indexed) and observed GT vector to each cell
	 * Index such as sc1 gets index 0, sc2 gets 1,...., scm gets m-1.
	 * @param cellGTVectorMap, HashMap<cellName, cellGTvector>. {c1:c1_GTV, c2:c2_GTV,....,cn:cn_GTV}
	 * @param singleCellNames, list containing the names of the single cells
	 * @param nCell, number of cells in the dataset, m
	 * @return
	 * Aug 4, 2017
	 */
	public ArrayList<SingleCell> constructSingleCellList(HashMap<String, Integer[]> cellGTVectorMap, 
														 ArrayList<String> singleCellNames, int nCell){
		ArrayList<SingleCell> listSingleCell = new ArrayList<>();
		for (int i = 0; i < nCell; i++){
			String cellName = singleCellNames.get(i);
			SingleCell singleCellObj = new SingleCell(cellName);
			singleCellObj.assignGTVectorID(i, cellGTVectorMap.get(cellName));
			listSingleCell.add(singleCellObj);
		}
		return listSingleCell;
	}
	
	/**
	 * Remove the empty clone (having index in the list) from the list of clones and update names and IDs of the clones
	 * @param listClone
	 * @param index
	 * @return
	 * Created On: Aug 4, 2017
	 */
//	public ArrayList<Clone> removeOneClone(ArrayList<Clone> listClone, int index){
//		listClone.remove(index);
//		updateNameIDCloneList(listClone);
//		return listClone;
//	}
	
	public void removeOneClone(ArrayList<Clone> listClone, int index){
		listClone.remove(index);
	}
	
	/**
	 * Return a new Clone list after removing the empty clones
	 * @param listClone
	 * @return
	 * Created On: Aug 11, 2017
	 */
	public ArrayList<Clone> removeEmptyClones(ArrayList<Clone> listClone){
		ArrayList<Clone> nuList = new ArrayList<>();
		for (Clone C : listClone){
			if (C.memberCellList.size() > 0)
				nuList.add(C);
		}
		return nuList;
	}
	
	/**
	 * Update the name and ID of the Clones in the list
	 * @param listClone
	 * Created On: Aug 4, 2017
	 */
	public void updateNameIDCloneList(ArrayList<Clone> listClone){
		for (int i = 0; i < listClone.size(); i++){
			listClone.get(i).setId(i);
		}
	}
	
	/**
	 * Obtain the genotypeObservation for clone genotypes
	 * To use for calculating the tree likelihood of clone genotypes
	 * @param listClone
	 * @param cloneNames
	 * @param nMut
	 * @return
	 * Created On: Aug 11, 2017
	 */
	public ArrayList<GenotypeObservation> getCloneGenotypeObs(ArrayList<Clone> listClone, ArrayList<String> cloneNames, int nMut){
		ArrayList<GenotypeObservation> cloneGenotypeObsMatrix = new ArrayList<>();
		for (int i = 0; i < nMut; i++){
			Integer[] thisMutVector = new Integer[listClone.size()];
			for (int j = 0; j < listClone.size(); j++){
				thisMutVector[j] = listClone.get(j).cloneGTVector[i];
			}
			cloneGenotypeObsMatrix.add(new GenotypeObservation(cloneNames, thisMutVector));
		}
		return cloneGenotypeObsMatrix;
	}
	
	
	public STITree<Clone> generateRandomCloneTree(ArrayList<String> cloneNames){
		return randomCloneTreeHelper(cloneNames, cloneNames.size());
	}
	
	private STITree<Clone> randomCloneTreeHelper(ArrayList<String> cloneNames, int nClone) {
		// If there is one cell, create a tree with one node
		if (nClone == 1){
			STITree<Clone> cell = new STITree<>(cloneNames.get(0), true);
			cell.getRoot().setParentDistance(rng.nextDouble());
			return cell;
		}
		
		// Randomly split the list of leafs into two lists.
		ArrayList<String> leftCells = new ArrayList<>();
		ArrayList<String> rightCells = new ArrayList<>();
		for (int i = 0; i < nClone; i++){
			if (rng.nextDouble() <= 0.5){
				leftCells.add(cloneNames.get(i));
			}
			else
				rightCells.add(cloneNames.get(i));				
		}
		
		// left split is empty
		if (leftCells.size() == 0){
			return randomCloneTreeHelper(rightCells, rightCells.size());
		}
		// right split is empty
		if (rightCells.size() == 0){
			return randomCloneTreeHelper(leftCells, leftCells.size());
		}
		
		// Root the randomly generated subtrees together.
		STITree<Clone> leftTree = randomCloneTreeHelper(leftCells, leftCells.size());
		STITree<Clone> rightTree = randomCloneTreeHelper(rightCells, rightCells.size());
		STITree<Clone> newTree = new STITree<>(true);
		newTree.getRoot().adoptChild(leftTree.getRoot());
        newTree.getRoot().adoptChild(rightTree.getRoot());
        
        // Randomly assign branch lengths between the subtrees' roots and the new tree's root.
        leftTree.getRoot().setParentDistance(rng.nextDouble());
        rightTree.getRoot().setParentDistance(rng.nextDouble());
        
        return newTree;
	}

	/**
	 * Get tree from newick string
	 * @param treeString
	 * @return
	 * Created On: Aug 11, 2017
	 */
	public STITree<Clone> getCloneTree (String treeString) {
		NewickReader nr = new NewickReader(new StringReader(treeString));
        STITree<Clone> tree = new STITree<Clone>(true);
        try {
            nr.readTree(tree);
        }
        catch(Exception ex) {
            ex.printStackTrace();
        }
        
        return tree;
	}
	
	/**
	 * Return the Clone where the cell belongs to
	 * @param listClone
	 * @param s
	 * @return
	 * Created On: Aug 12, 2017
	 */
	public Clone getCellClone(ArrayList<Clone> listClone, SingleCell s){
		for (Clone C: listClone){
			if (C.cloneID == s.cloneID)
				return C;
		}
		return null;
	}
	
	/**
	 * Return a new list of clones by appending the nuClone at the end of the old Clone list
	 * @param oldCloneList
	 * @param nuClone
	 * @return
	 * Created On: Aug 14, 2017
	 */
	public ArrayList<Clone> getNuCloneList(ArrayList<Clone> oldCloneList, Clone nuClone, Clone oldClone, SingleCell cell, int nMut){
		ArrayList<Clone> nuCloneList = new ArrayList<>();
		for (Clone C: oldCloneList){
			if (C.cloneID != oldClone.cloneID){
				Clone C_copy = this.copyClone(C, nMut);
				nuCloneList.add(C_copy);
			}
				
			else{
				Clone C_copy = this.copyClone(C, nMut);
				C_copy.removeCell(cell.cellID);
				nuCloneList.add(C_copy);
			}
		}			
		nuCloneList.add(nuClone);		
		return nuCloneList;
	}
	
	/**
	 * Add a clone to clone list
	 * @param oldCloneList
	 * @param nuClone
	 * @param nMut
	 * @return
	 * Created On: Oct 3, 2017
	 */
	public ArrayList<Clone> getNuCloneList(ArrayList<Clone> oldCloneList, Clone nuClone, int nMut){
		ArrayList<Clone> nuCloneList = new ArrayList<>();
		for (Clone C: oldCloneList){
			Clone C_copy = this.copyClone(C, nMut);
			nuCloneList.add(C_copy);
		}
		nuCloneList.add(nuClone);		
		return nuCloneList;
	}
	
	/**
	 * Returns a list of clones other than the one with position index in the original list
	 * @param listClone
	 * @param index
	 * @param nMut
	 * @return
	 * Created On: Aug 15, 2017
	 */
	public ArrayList<Clone> getOtherCloneList(ArrayList<Clone> listClone, int index, int nMut){
		ArrayList<Clone> nuCloneList = new ArrayList<>();
		for (int k = 0; k < listClone.size(); k++){
			if (k != index){
				Clone C_copy = this.copyClone(listClone.get(k), nMut);
				nuCloneList.add(C_copy);
			}
		}
		return nuCloneList;
	}
	
	/**
	 * Create a copy of the Clone
	 * @param C
	 * @param nMut
	 * @return
	 * Created On: Aug 15, 2017
	 */
	public Clone copyClone(Clone C, int nMut){
		Clone C_copy = new Clone(C.cloneID);
		C_copy.setNameID(C.cloneID);
		C_copy.cloneGTVector = new Integer[nMut];
		for (int i = 0; i < nMut; i++)
			C_copy.cloneGTVector[i] = C.cloneGTVector[i];
		C_copy.memberCellList = new ArrayList<>();
		for (int cell : C.memberCellList)
			C_copy.memberCellList.add(cell);
		return C_copy;
		
	}
	
	/**
	 * Return a copy of the list of clones
	 * @param cloneL
	 * @param nMut
	 * @return
	 * Created On: Oct 4, 2017
	 */
	public ArrayList<Clone> copyCloneList(ArrayList<Clone> cloneL, int nMut){
		ArrayList<Clone> nuCloneList = new ArrayList<>();
		for (Clone C: cloneL){
			Clone C_copy = copyClone(C, nMut);
			nuCloneList.add(C_copy);
		}
		return nuCloneList;
	}
	
	/**
	 * Create a copy of the single cell
	 * @param s
	 * @param nMut
	 * @return
	 * Created On: Oct 4, 2017
	 */
	public SingleCell copySingleCell(SingleCell s, int nMut){
		SingleCell s_copy = new SingleCell(s.name);
		s_copy.setID(s.cellID);
		s_copy.assignClone(s.cloneID);
		s_copy.observedGTVector = s.observedGTVector;
		s_copy.cloneGTVector = new Integer[nMut];
//		for (int i = 0; i < nMut; i++)
//			s_copy.cloneGTVector[i] = s.cloneGTVector[i];
		s_copy.cloneName = s.cloneName;
		return s_copy;
	}
	
	/**
	 * Create a copy of list of Cells
	 * @param cellList
	 * @param nMut
	 * @return
	 * Created On: Oct 4, 2017
	 */
	public ArrayList<SingleCell> copyCellList(ArrayList<SingleCell> cellList, int nMut){
		ArrayList<SingleCell> nuCellList = new ArrayList<>();
		for (SingleCell s : cellList){
			SingleCell s_copy = copySingleCell(s, nMut);
			nuCellList.add(s_copy);
		}
		return nuCellList;			
	}
	
	/**
	 * Return the list of clone names
	 * @param listClone
	 * @return
	 * Created On: Oct 4, 2017
	 */
	public ArrayList<String> getCloneNames(ArrayList<Clone> listClone){
		ArrayList<String> cloneNames = new ArrayList<>();
		for (Clone C: listClone){
			cloneNames.add(C.cloneName);
		}
		return cloneNames;
	}
	
	/**
	 * Update the list of single cells after updating the clone info of the cell
	 * @param listCells
	 * @param nuClone
	 * @param cell
	 * Created On: Aug 15, 2017
	 */
	public void updateSingleCellList(ArrayList<SingleCell> listCells, Clone nuClone, SingleCell cell){
		for (SingleCell s : listCells){
			if (s.cellID == cell.cellID){
				s.cloneID = nuClone.cloneID;
				s.cloneGTVector = nuClone.cloneGTVector;
			}
		}
	}
	
	public void updateSingleCellList(ArrayList<SingleCell> listCells, ArrayList<Clone> nuListClone){
		for (Clone C : nuListClone){
			for (int cellID : C.memberCellList){
				listCells.get(cellID).cloneID = C.cloneID;
				listCells.get(cellID).cloneName = C.cloneName;
				listCells.get(cellID).cloneGTVector = C.cloneGTVector;
			}
		}
	}
	
	public void updateCloneList(ArrayList<Clone> nuListClone, Clone C_j_new, int j){
		for (Clone C: nuListClone){
			if (C.cloneID == C_j_new.cloneID){
				C.addCell(j);
			}
		}
	}
	/**
	 * Return a new list of clone names by appending the nuClone at the end of the old Clone list
	 * @param oldCloneList
	 * @param nuClone
	 * @return
	 * Created On: Aug 14, 2017
	 */
	public ArrayList<String> getNuCloneNameList(ArrayList<String> oldCloneList, String nuClone){
		ArrayList<String> nuCloneList = new ArrayList<>();
		for (String C: oldCloneList)
			nuCloneList.add(C);
		nuCloneList.add(nuClone);
		return nuCloneList;
	}
	
	/**
	 * Sample the concentration parameter during Gibbs Sampling
	 * Using the auxiliary variable method of Escobar and West
	 * @param currAlpha0, current value of alpha0
	 * @param nClones, k
	 * @param nCells, n
	 * @param GammaA, a
	 * @param GammaB, b
	 * @return
	 * Created On: Aug 8, 2017
	 */
	public double sampleConcentrationParam(double currAlpha0, int nClones, int nCells, double GammaA, double GammaB){
		// ita ~ Beta(alpha+1,n)
		double ita = new BetaDistribution(currAlpha0+1, nCells).sample();
		double categoricalPi = (GammaA + nClones - 1)/(GammaA + nClones - 1 + nCells * (GammaB - Math.log(ita)));
		double rr = new Random().nextDouble();
		if (rr <= categoricalPi){
			double alphaNu = new GammaDistribution(GammaA + nClones, GammaB - Math.log(ita)).sample();
			return alphaNu;
		}
		else {
			double alphaNu = new GammaDistribution(GammaA + nClones - 1, GammaB - Math.log(ita)).sample();
			return alphaNu;
		}
		
//		double t1 = new GammaDistribution(GammaA + nClones, GammaB - Math.log(ita)).sample();
//		double t2 = new GammaDistribution(GammaA + nClones - 1, GammaB - Math.log(ita)).sample();
////		double alphaNu = ((GammaA + nClones - 1) * t1 + nCells * (GammaB - Math.log(ita)) * t2)/(GammaA + nClones - 1 + nCells * (GammaB - Math.log(ita)));
//		return alphaNu;
	}
	
	/**
	 * Write genotype matrix to file
	 * @param filename
	 * @param cellList
	 * @param cellsGenotypeArrMap
	 * @param genotypeMat
	 * @param nMut
	 * @throws FileNotFoundException
	 * @throws UnsupportedEncodingException
	 * Created On: Sep 22, 2017
	 */
	public void writeGenotypeMatrix2File(String filename, ArrayList<Integer> cellList,
									     HashMap<Integer, Integer[]> cellsGenotypeArrMap,
									     Integer[][] genotypeMat, int nMut) throws FileNotFoundException, UnsupportedEncodingException{
		PrintWriter writer = new PrintWriter(filename, "UTF-8");
		for (int i = 0; i < nMut; i++){
			String this_mut_row = Integer.toString(i);
			genotypeMat[i][0] = i;
			for (int j = 0; j < cellList.size(); j++){
				genotypeMat[i][j+1] = cellsGenotypeArrMap.get(j)[i];
				this_mut_row = this_mut_row + " " + String.valueOf(genotypeMat[i][j+1]);
			}
			writer.printf("%s\n", this_mut_row);
		}
		writer.close();
	}
	
	/**
	 * Write genotype matrix to file
	 * @param filename
	 * @param nCell
	 * @param cellsGenotypeArrMap
	 * @param genotypeMat
	 * @param nMut
	 * @throws FileNotFoundException
	 * @throws UnsupportedEncodingException
	 * Created On: Oct 30, 2017
	 */
	public void writeGenotypeMatrix2File(String filename, int nCell,
		     HashMap<Integer, Integer[]> cellsGenotypeArrMap,
		     Integer[][] genotypeMat, int nMut) throws FileNotFoundException, UnsupportedEncodingException{
		PrintWriter writer = new PrintWriter(filename, "UTF-8");
		for (int i = 0; i < nMut; i++){
			String this_mut_row = Integer.toString(i);
			genotypeMat[i][0] = i;
			for (int j = 0; j < nCell; j++){
				genotypeMat[i][j+1] = cellsGenotypeArrMap.get(j)[i];
				this_mut_row = this_mut_row + " " + String.valueOf(genotypeMat[i][j+1]);
			}
			writer.printf("%s\n", this_mut_row);
		}
		writer.close();
	}
	
	public void writeSCGGenotypeMatrix2File(String filename, ArrayList<ArrayList<String>> scgMatrix) throws UnsupportedEncodingException, IOException{
		FileOutputStream output = new FileOutputStream(filename);
		try {
			Writer writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
			try {
				for (ArrayList<String> line : scgMatrix){
					String line2p = String.join("\t", line);
					line2p += "\n";
					writer.write(line2p);
				}
			} finally {
				writer.close();
			}
		} finally {
			output.close();
		}
	}
	
	/**
	 * Get the index of the maximum element in array
	 * @param array
	 * @return
	 * Created On: Aug 11, 2017
	 */
	public int getMaxIndexArray(int[] array){
		if (array.length == 0) {
	        return -1; // array contains no elements
	    }
	    int max = array[0];
	    int pos = 0;

	    for(int i=1; i<array.length; i++) {
	        if (max < array[i]) {
	            pos = i;
	            max = array[i];
	        }
	    }
	    return pos;
	}

	/**
	 * @param args
	 * Aug 4, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		System.out.print("hello");
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
		for (int i = 1; i < 10; i++){
		double a = SCF.sampleConcentrationParam(1, i, 40, 1, 1);
		System.out.println(a);
		}

	}

	public ArrayList<String> getNuCloneNameList(ArrayList<Clone> nuCloneList) {
		ArrayList<String> cloneNames =  new ArrayList<>();
		for (Clone C: nuCloneList)
			cloneNames.add(C.cloneName);
		return cloneNames;
	}

	/**
	 * Assign 2nd clone randomly to cell s_j
	 * @param s_j
	 * @param listClone
	 * @param cloneNames
	 * @return
	 * Created On: Oct 8, 2017
	 */
	public int assignRandomClone2(SingleCell s_j, ArrayList<Clone> listClone) {
		ArrayList<Integer> otherCloneList = new ArrayList<>();
		for (Clone C: listClone){
			if (C.cloneID != s_j.cloneID)
				otherCloneList.add(C.cloneID);
		}
		int nuCloneIndex = rng.nextInt(otherCloneList.size());
		return otherCloneList.get(nuCloneIndex);
	}
	
	/**
	 * Return the Clone that has same ID as cloneID
	 * @param listClone
	 * @param cloneID
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public Clone getCloneFrmList(ArrayList<Clone> listClone, int cloneID){
		for (Clone C : listClone){
			if (C.cloneID == cloneID)
				return C;
		}
		return null;
	}
	
	/**
	 * Create a copy of list of integers
	 * @param c2CellCloneList
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public ArrayList<Integer> getCopyC2CellCloneList(ArrayList<Integer> c2CellCloneList){
		ArrayList<Integer> copyC2CellCloneList = new ArrayList<>();
		for (Integer i : c2CellCloneList)
			copyC2CellCloneList.add(i);
		return copyC2CellCloneList;
	}
	
	/**
	 * Create a copy of integer array
	 * @param arr
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public int[] getCopyCellDoubletFlagArr(int[] arr){
		int[] copyArr = new int[arr.length];
		for (int i = 0; i < arr.length; i++)
			copyArr[i] = arr[i];
		return copyArr;
	}
	
	public Integer[] getMutatedGenotypeArrFSMNu(
			STINode<Double> child,
			Integer[] parentGenotypeArr,
			HashMap<Integer, Set<Integer>> MutationTypeMap,
			Integer[] mutFlagArr,
			double recurProb, double omega,
			double deletion, double branchLength,
			int nMut){
		Integer[] childGenotypeArr = new Integer[nMut];
		int nuMut2introduce;
		if (!child.isLeaf()){
//			int PoissonParam = (int) Math.round(branchLength*nMut);
//			if (PoissonParam == 0)
//				nuMut2introduce = (int) Math.round(branchLength*nMut);
//			else{
			int PoissonParam = 1+(int) Math.round(branchLength*nMut);
//				PoissonParam += 1;
			if (PoissonParam == 0)
				PoissonParam = 1;
			nuMut2introduce = new PoissonDistribution(PoissonParam).sample();
//			}
		}
		else
			nuMut2introduce = (int) Math.round(branchLength*nMut);
//		System.out.printf("param=%s\tmut=%d\n", child.getName(), nuMut2introduce);
		
		// LOH are introduced
		for (int i = 0; i< parentGenotypeArr.length; i++){
			if (parentGenotypeArr[i] == 1){
				double rrw = rng.nextDouble();
				if (rrw <= omega){
					System.out.println("LOH happening");
					double rrw1 = rng.nextDouble();
					if (rrw1 <= 0.5){
					childGenotypeArr[i] = 2;
					MutationTypeMap.get(2).add(i);
					}
					else{
						childGenotypeArr[i] = 0;
					}
				}
				else
					childGenotypeArr[i] = 1;
			}
			else
				childGenotypeArr[i] = parentGenotypeArr[i];
		}
		double rr = rng.nextDouble();
		// Deletion happening
		if (rr <= deletion){
			double rr1 = rng.nextDouble();
			// Deletion results in homozygous reference genotype
			if (rr1 <= 0.5){
				for (int i = 0; i < nuMut2introduce; i++){
					if (MutationTypeMap.get(1).size() == 0)
						continue;
					ArrayList<Integer> mutatedPosList = new ArrayList<>();
					mutatedPosList.addAll(MutationTypeMap.get(1));
					int indexInSet = rng.nextInt(mutatedPosList.size());
					int site_index = mutatedPosList.get(indexInSet);
					if (parentGenotypeArr[site_index] == 1){
						childGenotypeArr[site_index] = 0;
						System.out.println("deletion to 0");
					}
					else if (parentGenotypeArr[site_index] == 2)
						childGenotypeArr[site_index] = 0;
					else
						childGenotypeArr[site_index] = parentGenotypeArr[site_index];
				}
			}
			// Deletion results in homozygous variant genotype
			else{
				for (int i = 0; i < nuMut2introduce; i++){
					if (MutationTypeMap.get(1).size() == 0)
						continue;
					ArrayList<Integer> mutatedPosList = new ArrayList<>();
					mutatedPosList.addAll(MutationTypeMap.get(1));
					int indexInSet = rng.nextInt(mutatedPosList.size());
					int site_index = mutatedPosList.get(indexInSet);
					if (parentGenotypeArr[site_index] == 1){
						childGenotypeArr[site_index] = 2;
						System.out.println("deletion to 2");
					}
					else if (parentGenotypeArr[site_index] == 2)
						childGenotypeArr[site_index] = 0;
					else
						childGenotypeArr[site_index] = parentGenotypeArr[site_index];
				}
			}
		}
		else{
			for (int i = 0; i < nuMut2introduce; i++){
				double rr2 = rng.nextDouble();
				// Recurrent mutation happens
				if (rr2 <= recurProb){
					if (MutationTypeMap.get(1).size() == 0)
						continue;
					System.out.println("recurrent mutation happening");
					ArrayList<Integer> mutatedPosList = new ArrayList<>();
					mutatedPosList.addAll(MutationTypeMap.get(1));
					int indexInSet = rng.nextInt(mutatedPosList.size());
					int site_index = mutatedPosList.get(indexInSet);
					if (parentGenotypeArr[site_index] == 1)
						childGenotypeArr[site_index] = 0;         // back mutation
					else if (parentGenotypeArr[site_index] == 0)  // parallel mutation
						childGenotypeArr[site_index] = 1;
					else
						childGenotypeArr[site_index] = 1;
				}
				else{
					ArrayList<Integer> mutatedPosList = new ArrayList<>();
					// Nu mutation happens
					if (MutationTypeMap.get(0).size() > 0){
						mutatedPosList.addAll(MutationTypeMap.get(0));
//						System.out.println(mutatedPosList.size());
						int indexInSet = rng.nextInt(mutatedPosList.size());
						int site_index = mutatedPosList.get(indexInSet);
						childGenotypeArr[site_index] = 1;
						MutationTypeMap.get(0).remove(site_index);
						MutationTypeMap.get(1).add(site_index);
					}
					// No more new mutation, recurrent mutation happens
					else{
						mutatedPosList.addAll(MutationTypeMap.get(1));
						System.out.println("no nu site to mutate, recurrent mutation");
						int indexInSet = rng.nextInt(mutatedPosList.size());
						int site_index = mutatedPosList.get(indexInSet);
						if (parentGenotypeArr[site_index] == 1)
							childGenotypeArr[site_index] = 0;         // back mutation
						else if (parentGenotypeArr[site_index] == 0)  // parallel mutation
							childGenotypeArr[site_index] = 1;
						else
							childGenotypeArr[site_index] = 1;
					}
				}
			}
			
		}
//		for (int i = 0; i < nuMut2introduce; i++){
//			
//			
//			// Recurrent mutation
//			else if (rr <= 0.01 + recurProb){
//				if (MutationTypeMap.get(1).size() == 0)
//					continue;
//				System.out.println("recurrent mutation happening");
//				ArrayList<Integer> mutatedPosList = new ArrayList<>();
//				mutatedPosList.addAll(MutationTypeMap.get(1));
//				int indexInSet = _rng.nextInt(mutatedPosList.size());
//				int site_index = mutatedPosList.get(indexInSet);
//				if (parentGenotypeArr[site_index] == 1)
//					childGenotypeArr[site_index] = 0;         // back mutation
//				else if (parentGenotypeArr[site_index] == 0)  // parallel mutation
//					childGenotypeArr[site_index] = 1;
//				else
//					childGenotypeArr[site_index] = 1;
//			}
//			else{
//				ArrayList<Integer> mutatedPosList = new ArrayList<>();
//				if (MutationTypeMap.get(0).size() > 0){
//					mutatedPosList.addAll(MutationTypeMap.get(0));
////					System.out.println(mutatedPosList.size());
//					int indexInSet = _rng.nextInt(mutatedPosList.size());
//					int site_index = mutatedPosList.get(indexInSet);
//					childGenotypeArr[site_index] = 1;
//					MutationTypeMap.get(0).remove(site_index);
//					MutationTypeMap.get(1).add(site_index);
//				}
//				else{
//					mutatedPosList.addAll(MutationTypeMap.get(1));
////					System.out.println(mutatedPosList.size());
//					int indexInSet = _rng.nextInt(mutatedPosList.size());
//					int site_index = mutatedPosList.get(indexInSet);
//					if (parentGenotypeArr[site_index] == 1)
//						childGenotypeArr[site_index] = 0;         // back mutation
//					else if (parentGenotypeArr[site_index] == 0)  // parallel mutation
//						childGenotypeArr[site_index] = 1;
//					else
//						childGenotypeArr[site_index] = 1;
//				}
//			}
//		}
		return childGenotypeArr;
	}

	/**
	 * Add missing data to a genotype matrix
	 * @param ingenotypeMat
	 * @param outgenotypeMat
	 * @param missing
	 * Created On: Oct 25, 2017
	 */
	public void addMissingData2GenotypeMat(Integer[][] ingenotypeMat, Integer[][] outgenotypeMat, double missing){
		int nCell = ingenotypeMat[0].length;
		for (int i = 0; i < ingenotypeMat.length; i++){
			outgenotypeMat[i][0] = ingenotypeMat[i][0];
			for (int j = 1; j < nCell; j++){
				double rr = rng.nextDouble();
				if (rr <= missing)
					outgenotypeMat[i][j] = 3;
				else
					outgenotypeMat[i][j] = ingenotypeMat[i][j];
			}
		}
	}
	
	/**
	 * Print genotype matrix to a file
	 * @param filename
	 * @param genotypeMat
	 * @throws FileNotFoundException
	 * @throws UnsupportedEncodingException
	 * Created On: Oct 25, 2017
	 */
	public void writeGenotypeMat2File(String filename, Integer[][] genotypeMat) throws FileNotFoundException, UnsupportedEncodingException{
		PrintWriter writer = new PrintWriter(filename, "UTF-8");
		for (int i = 0; i < genotypeMat.length; i++){
			String this_mut_row = String.valueOf(genotypeMat[i][0]);
			for (int j = 0; j < genotypeMat[0].length-1; j++){
				this_mut_row = this_mut_row + " " + String.valueOf(genotypeMat[i][j+1]);
			}
			writer.printf("%s\n", this_mut_row);
		}
		writer.close();
	}
}
