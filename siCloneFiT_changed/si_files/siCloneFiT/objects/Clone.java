/**
 * Aug 4, 2017
 */
package objects;

import java.util.ArrayList;
import java.util.NoSuchElementException;

import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Aug 4, 2017
 */
public class Clone {

	/**
	 * @param args
	 * Aug 4, 2017
	 */
	public String cloneName;
	public int cloneID;
	public Integer[] cloneGTVector;
	public ArrayList<Integer> memberCellList = new ArrayList<>();
	
	public Clone(int ID){
		this.cloneID = ID;
	}
	
	/**
	 * Assign the ID of the Clone
	 * @param ID
	 * @return
	 * Aug 4, 2017
	 */
	public boolean setId(int ID){
		this.cloneID = ID;
		return true;
	}
	
	/**
	 * Assign Name of the Clone
	 * @param ID
	 * @return
	 * Aug 4, 2017
	 */
	public boolean setName(int ID){
		this.cloneName = "C" + Integer.toString(ID);
		return true;
	}
	
	/**
	 * Assign name and ID of a Clone
	 * @param ID
	 * @return
	 * Created On: Aug 4, 2017
	 */
	public boolean setNameID(int ID){
		this.setName(ID);
		this.setId(ID);
		return true;
	}
	
	/**
	 * Assign a random Genotype Vector to the clone
	 * Each locus has the genotype that gets maximum vote 
	 * among the cells of this clone
	 * @param nMut
	 * @param listSingleCell
	 * @param SCF
	 * @return
	 * Created On: Aug 11, 2017
	 */
	public boolean assignRandomGTVector(int nMut, ArrayList<SingleCell> listSingleCell, SCFUtilityFunctions SCF){
		this.cloneGTVector = new Integer[nMut];
		for (int i = 0; i < nMut; i++){
			int[] cellGTCount = new int[3]; // Array that keeps count of genotypes
			for (Integer cellID: this.memberCellList){
				int gt = listSingleCell.get(cellID).observedGTVector[i];
				if (gt != 3){ 
//					System.out.println(gt);
					cellGTCount[gt] += 1;
				}
				else
					cellGTCount[0] += 1;				
			}
			int cloneGt_i = SCF.getMaxIndexArray(cellGTCount);	
			cloneGTVector[i] = cloneGt_i;
		}
		return true;
	}
	
	/**
	 * Add a cell to the Clone
	 * @param cell, SingleCell Object
	 * Created On: Aug 4, 2017
	 */
	public void addCell(SingleCell cell){
		memberCellList.add(cell.cellID);
	}
	
	/**
	 * Add a cell to the Clone
	 * @param cellID, ID of the cell
	 * Created On: Aug 4, 2017
	 */
	public void addCell(int cellID){
		memberCellList.add(cellID);
	}
	
	/**
	 * Remove a cell from the Clone
	 * @param cellID
	 * Created On: Aug 4, 2017
	 */
	public void removeCell(int cellID){
		int index = memberCellList.indexOf(cellID);
		if (index == -1){
			throw new NoSuchElementException("Invalid membership, cell not in clone");
		}
		memberCellList.remove(index);
	}
	
	/**
	 * Check if the Clone has only one member
	 * @return
	 * Created On: Aug 12, 2017
	 */
	public boolean isSingleton(){
		if (memberCellList.size() == 1)
			return true;
		else
			return false;
	}
	
	
	
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
		
		ArrayList<Clone> cloneList = new ArrayList<>();
		for (int i = 0; i < 4; i++){
			Clone clone = new Clone(i);
			clone.setNameID(i);
			for (int j = i; j < i+i+1; j++){
				clone.addCell(j);
			}
//			System.out.println(clone.memberCellList);
			cloneList.add(clone);
		}
		cloneList.get(0).removeCell(0);
		
		
		
		ArrayList<Clone> nuclones = new ArrayList<>();
		for (Clone c : cloneList){
			if (c.memberCellList.size() > 0){
				nuclones.add(c);
			}
			System.out.printf("clone %s has ID = %d\n", c.cloneName, c.cloneID);
			System.out.println(c.memberCellList);
		}
		System.out.println(" ");
		
		SCF.removeOneClone(cloneList, 1);
		
		for (int i = 0; i < nuclones.size(); i++){
			nuclones.get(i).setNameID(i);
		}
		for (Clone c: cloneList){
			System.out.printf("clone %s has ID = %d\n", c.cloneName, c.cloneID);
			System.out.println(c.memberCellList);
		}

	}

}
