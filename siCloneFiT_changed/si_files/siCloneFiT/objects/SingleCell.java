/**
 * Aug 4, 2017
 */
package objects;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import siCloneFiT.io.GenotypeMatrixReader;

/**
 * @author hz22
 * Aug 4, 2017
 */
public class SingleCell {

	/**
	 * @param args
	 * Aug 4, 2017
	 */
	public String name;
	public int cellID;
	public Integer[] observedGTVector;
	public int cloneID;
	public String cloneName;
	public Integer[] cloneGTVector;
	
	public SingleCell(String cellName){
		this.name = cellName;
	}
	
	/**
	 * Assign the ID of the cell
	 * @param cID
	 * @return
	 * Aug 4, 2017
	 */
	public boolean setID(int cID){
		this.cellID = cID;
		return true;
	}
	
	/**
	 * Assign SingleCell to a clone with index cloneIndex
	 * @param cloneIndex
	 * @return
	 * Aug 4, 2017
	 */
	public boolean assignClone(int cloneIndex){
		this.cloneID = cloneIndex;
		return true;
	}
	
	/**
	 * Assign Observed GT Vector to the cell
	 * @param cellGTVector
	 * @return
	 * Aug 4, 2017
	 */
	public boolean assignObsGTVector(Integer[] cellGTVector){
		this.observedGTVector = cellGTVector;
		return true;
	}
	
	/**
	 * Set ID and assign observed GT Vector to the cell
	 * @param ID
	 * @param cellGTVector
	 * @return
	 * Aug 4, 2017
	 */
	public boolean assignGTVectorID(int ID, Integer[] cellGTVector){
		this.setID(ID);
		this.assignObsGTVector(cellGTVector);
		return true;
	}
	
	
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		// This is a test file
//		String cell = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/examples/10cells_cellnames.txt";
		String GtMatFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/examples/10cells.txt";

		int nMut = 10;
		int nCell = 10;
		GenotypeMatrixReader GMR = new GenotypeMatrixReader(GtMatFile);
		HashMap<String, Integer[]> gt = GMR.getCellGTVectorMap(GtMatFile, null, nMut, nCell);
		ArrayList<String> singleCellNames = GMR.scNames;
		
		ArrayList<SingleCell> listSingleCell = new ArrayList<>();
		for (int i = 0; i < nCell; i++){
			SingleCell sc = new SingleCell(singleCellNames.get(i));
			sc.assignGTVectorID(i, gt.get(singleCellNames.get(i)));
			listSingleCell.add(sc);
		}
		for (int i = 0; i < nCell; i++){
			System.out.printf("cell %s has ID = %d\n", listSingleCell.get(i).name, listSingleCell.get(i).cellID);
			System.out.println(Arrays.toString(listSingleCell.get(i).observedGTVector));
		}

//		ArrayList<Integer> test = new ArrayList<>();
//		for (int i = 0; i < nCell; i++){
//			test.add(i);
//		}
//		System.out.println(test);
	}

}
