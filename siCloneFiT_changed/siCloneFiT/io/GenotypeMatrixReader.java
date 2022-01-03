/**
 * 
 */
package siCloneFiT.io;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import SiFit.io.VariantMatrixReader;
import cern.colt.Arrays;-

/**
 * @author hz22
 *
 */
public class GenotypeMatrixReader extends VariantMatrixReader {

	/**
	 * @param filename
	 */
	public GenotypeMatrixReader(String filename) {
		super(filename);
		// TODO Auto-generated constructor stub
	}
	
	/**
	 * Construct a HashMap<cellName, cellGTvector>. {c1:c1_GTV, c2:c2_GTV,....,cn:cn_GTV}
	 * @param GTMatrixFile, File containing GenotypeMatrix (n by m+1 matrix)
	 * @param cellNamesFile, File containing the names of the cells in a single row
	 * @param nMut, number of mutations, n
	 * @return HashMap<cellName, cellGTvector>
	 * @throws IOException
	 * Aug 3, 2017
	 */
	public HashMap<String, Integer[]> getCellGTVectorMap(String GTMatrixFile, String cellNamesFile, int nMut, int nCell) throws IOException{
		HashMap<String, Integer[]> cellGTVectorMap = new HashMap<>();
		String[] a="";
		// Process the names of the cells
		// Cell names not given
		if (cellNamesFile == null){
			for (int i = 1; i <= nCell; i++){
            	this.scNames.add("sc" + Integer.toString(i-1));
            }
		}
		// Cell names given
		else{
			BufferedReader brCells = new BufferedReader(new FileReader(cellNamesFile));
	        while ((line = brCells.readLine()) != null) {
            	a=line.split(",");
                this.scNames.add(a[0]);
            }
			brCells.close();
		}
		
		// Initiate the hashmap
		for (String cell : this.scNames){
			Integer[] cellGtVector = new Integer[nMut];
			cellGTVectorMap.put(cell, cellGtVector);
		}
		// Populate cellGTVectorMap from the input genotype matrix file
		try (BufferedReader br = new BufferedReader(new FileReader(GTMatrixFile))) {
			String line=br.readLine();
            while ((line = br.readLine()) != null) { 
            	String[] thisMutVector = line.split("	");
            	for (int i = 1; i <= nMut; i++){
            		cellGTVectorMap.get(thisMutVector[0])[i] = Integer.parseInt(thisMutVector[i]);
            	}
            }
		}
		return cellGTVectorMap;
	}
	

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
//		String cell = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/colon_cancer_mutation_matrix/CO5_cellnames.txt";
//		String GtMatFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/colon_cancer_mutation_matrix/CO5.master_matrix.sifit";
		
		// This is a test file
		String cell= "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/examples/10cells_cellnames.txt";
		String DMatFile= "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/examples/10cells.txt";

		int nMut = 10;
		int nCell = 10;
		HashMap<String, Integer[]> d_mat = new GenotypeMatrixReader(GtMatFile).getCellGTVectorMap(DMatFile, null, nMut, nCell);
		for (String c: d_mat.keySet()){
			System.out.println(c);
			System.out.println(Arrays.toString(d_mat.get(c)));
		}
			
	}

}
