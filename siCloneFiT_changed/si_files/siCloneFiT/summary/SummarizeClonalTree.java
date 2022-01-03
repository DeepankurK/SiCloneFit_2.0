package summary;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.ArrayList;
//import java.util.Collections;
//import java.util.Comparator;
//import java.util.LinkedHashMap;
//import java.util.LinkedList;
//import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import SiFit.metric.CompareTrees;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import siCloneFiT.metric.ClonalTreeDistance;
import siCloneFiT.utils.SCFUtilityFunctions;

public class SummarizeClonalTree {
	
	public String treeFile;
	public String cellNamesFile;
	
	public SummarizeClonalTree(String tfile, String cfile){
		this.treeFile = tfile;
		this.cellNamesFile = cfile;
	}
	
	public String getSameTree(String candidateTree, HashMap<String, Integer> treeMap, SCFUtilityFunctions SCF, int nCell) throws IOException{
		if (treeMap.isEmpty()){
			treeMap.put(candidateTree, 1);
			return candidateTree;
		}
		else{
			STITree<Double> icandidateTree = ClonalTreeDistance.getIntegerBLTree(SCF.getTree(candidateTree));
			for (String tree: treeMap.keySet()){
				STITree<Double> thisTree = SCF.getTree(tree);
				STITree<Double> ithisTree =  ClonalTreeDistance.getIntegerBLTree(thisTree);
				double dist = ClonalTreeDistance.getPairwiseCellSPTreeDist(ithisTree, icandidateTree, this.cellNamesFile, nCell);
				if (dist == 0.0){
					treeMap.put(tree, treeMap.get(tree) + 1);
					return candidateTree;
				}
			}
		}
		treeMap.put(candidateTree, 1);
		return candidateTree;
	}
	
	/**
	 * Check if a candidateTree exists in the treeSet 
	 * @param candidateTree
	 * @param treeMap
	 * @param SCF
	 * @param nCell
	 * @return
	 * @throws IOException
	 * Created On: Mar 23, 2018
	 */
	public String getSameTreeSD(String candidateTree, HashMap<String, Integer> treeMap, SCFUtilityFunctions SCF, int nCell) throws IOException{
		if (treeMap.isEmpty()){
			treeMap.put(candidateTree, 1);
			return candidateTree;
		}
		else{
			STITree<Double> icandidateTree = SCF.getTree(candidateTree);
			double FPdist, FNdist;
			for (String tree: treeMap.keySet()){
				STITree<Double> thisTree = SCF.getTree(tree);
				CompareTrees CTObj = new CompareTrees(thisTree);
				FPdist = CTObj.calcFPDist(icandidateTree);
				FNdist = CTObj.calcFNDist(icandidateTree);
				if (FPdist == 0.0 && FNdist == 0.0){
					treeMap.put(tree, treeMap.get(tree) + 1);
					return candidateTree;
				}
			}			
		}
		treeMap.put(candidateTree, 1);
		return candidateTree;
	}
	
	public String getSameTreeSD(String candidateTree, HashMap<String, ArrayList<String>> treeListMap, HashMap<String, Integer> treeMap, SCFUtilityFunctions SCF) throws IOException{
		if (treeMap.isEmpty()){
			ArrayList<String> thisTreeList = new ArrayList<>();
			thisTreeList.add(candidateTree);
			treeListMap.put(candidateTree, thisTreeList);
			treeMap.put(candidateTree, 1);
			return candidateTree;
		}
		else{
			STITree<Double> icandidateTree = SCF.getTree(candidateTree);
			double FPdist, FNdist;
			for (String tree: treeMap.keySet()){
				STITree<Double> thisTree = SCF.getTree(tree);
				CompareTrees CTObj = new CompareTrees(thisTree);
				FPdist = CTObj.calcFPDist(icandidateTree);
				FNdist = CTObj.calcFNDist(icandidateTree);
				if (FPdist == 0.0 && FNdist == 0.0){
					treeMap.put(tree, treeMap.get(tree) + 1);
					treeListMap.get(tree).add(candidateTree);
					return candidateTree;
				}
			}			
		}
		treeMap.put(candidateTree, 1);
		ArrayList<String> thisTreeList = new ArrayList<>();
		thisTreeList.add(candidateTree);
		treeListMap.put(candidateTree, thisTreeList);
		return candidateTree;
	}
	/*
    private Map<String, Integer> sortByComparator(Map<String, Integer> unsortMap, final boolean order)
    {

        List<Entry<String, Integer>> list = new LinkedList<Entry<String, Integer>>(unsortMap.entrySet());

        // Sorting the list based on values
        Collections.sort(list, new Comparator<Entry<String, Integer>>()
        {
            public int compare(Entry<String, Integer> o1,
                    Entry<String, Integer> o2)
            {
                if (order)
                {
                    return o1.getValue().compareTo(o2.getValue());
                }
                else
                {
                    return o2.getValue().compareTo(o1.getValue());

                }
            }
        });

        // Maintaining insertion order with the help of LinkedList
        Map<String, Integer> sortedMap = new LinkedHashMap<String, Integer>();
        for (Entry<String, Integer> entry : list)
        {
            sortedMap.put(entry.getKey(), entry.getValue());
        }

        return sortedMap;
    }
    */
    public void printMap(Map<String, Integer> map)
    {
    	ArrayList<Integer> values = new ArrayList<>();
    	ArrayList<String> keys = new ArrayList<>();
        for (Entry<String, Integer> entry : map.entrySet())
        {
//        	if (values.size() > 10){
        		
//        		System.out.println(" Value : "+ entry.getValue());
        		System.out.println(entry.getKey());
        		values.add(entry.getValue());
        		keys.add(entry.getKey());
//        		break;
//        	}
        		for (int i = 0; i<10; i++){
        			System.out.println(values.get(i));
        			System.out.println(keys.get(i));
        		}
            
        }
    }
    
    public void writeValidTrees(String ipFile, String opFile, SCFUtilityFunctions SCF) throws FileNotFoundException, IOException{
    	try (BufferedReader br = new BufferedReader(new FileReader(ipFile))) {
			String line;
			FileWriter fileWriter = new FileWriter(opFile);
			PrintWriter printWriter = new PrintWriter(fileWriter);
			while ((line = br.readLine()) != null) { 
				String thisTreeS = line.trim();
				STITree<Double> thisTree = SCF.getTree(thisTreeS);
//				System.out.println(thisTree.getNodeCount());
				if (thisTree.getNodeCount() == 197){
					System.out.println(thisTree.getNodeCount());
					printWriter.printf("%s\n", thisTreeS);
				}
			}
			printWriter.close();
    	}
    }

	public static void main(String[] args) throws FileNotFoundException, IOException {
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
		// String t = "(((((MA_27:1.0E-15,MA_30:1.0E-15,MA_31:1.0E-15,MA_32:1.0E-15,MA_34:1.0E-15,MA_36:1.0E-15,MA_38:1.0E-15,MA_39:1.0E-15,MA_41:1.0E-15,MA_42:1.0E-15,MA_43:1.0E-15,MA_44:1.0E-15,MA_46:1.0E-15,MA_47:1.0E-15,MA_48:1.0E-15,MA_87:1.0E-15,MA_88:1.0E-15,MA_89:1.0E-15,MA_90:1.0E-15,MA_92:1.0E-15,MA_93:1.0E-15,MA_94:1.0E-15,MA_95:1.0E-15)C2:0.7377703242240996,(MD_12:1.0E-15,PA_26:1.0E-15,PA_27:1.0E-15,PA_28:1.0E-15,PA_29:1.0E-15,PA_32:1.0E-15,PA_38:1.0E-15,PA_39:1.0E-15,PA_40:1.0E-15,PA_41:1.0E-15,PA_42:1.0E-15,PA_43:1.0E-15,PA_46:1.0E-15,PA_47:1.0E-15,PA_48:1.0E-15,PA_49:1.0E-15,PA_51:1.0E-15,PA_52:1.0E-15,PA_53:1.0E-15,PA_55:1.0E-15,PA_57:1.0E-15,PA_58:1.0E-15,PA_59:1.0E-15,PA_61:1.0E-15,PA_62:1.0E-15,PA_63:1.0E-15,PA_64:1.0E-15,PA_65:1.0E-15,PA_68:1.0E-15,PA_70:1.0E-15,PA_71:1.0E-15,PA_73:1.0E-15,PA_74:1.0E-15,PA_77:1.0E-15)C3:0.09006529665215349):0.00837874433941705,(MA_28:1.0E-15,MA_29:1.0E-15,MA_33:1.0E-15,MA_35:1.0E-15,MA_37:1.0E-15,MA_40:1.0E-15,MA_45:1.0E-15,MA_85:1.0E-15,MA_86:1.0E-15,MA_91:1.0E-15)C1:0.8095991683221228):0.60350617160311,(PD_1:1.0E-15,PD_14:1.0E-15,PD_16:1.0E-15,PD_17:1.0E-15,PD_18:1.0E-15,PD_20:1.0E-15,PD_9:1.0E-15)C4:0.06595838766664627):0.07164617725088579,(MA_25:1.0E-15,MA_26:1.0E-15,MD_10:1.0E-15,MD_11:1.0E-15,MD_13:1.0E-15,MD_3:1.0E-15,MD_4:1.0E-15,MD_5:1.0E-15,MD_7:1.0E-15,MD_8:1.0E-15,MD_9:1.0E-15,PA_54:1.0E-15,PD_5:1.0E-15,PD_7:1.0E-15,PD_8:1.0E-15,PDD_1:1.0E-15,PDD_10:1.0E-15,PDD_11:1.0E-15,PDD_12:1.0E-15,PDD_13:1.0E-15,PDD_14:1.0E-15,PDD_15:1.0E-15,PDD_16:1.0E-15,PDD_17:1.0E-15,PDD_18:1.0E-15,PDD_19:1.0E-15,PDD_2:1.0E-15,PDD_20:1.0E-15,PDD_22:1.0E-15,PDD_23:1.0E-15,PDD_24:1.0E-15,PDD_25:1.0E-15,PDD_26:1.0E-15,PDD_27:1.0E-15,PDD_28:1.0E-15,PDD_3:1.0E-15,PDD_30:1.0E-15,PDD_31:1.0E-15,PDD_32:1.0E-15,PDD_33:1.0E-15,PDD_34:1.0E-15,PDD_35:1.0E-15,PDD_36:1.0E-15,PDD_37:1.0E-15,PDD_38:1.0E-15,PDD_39:1.0E-15,PDD_4:1.0E-15,PDD_40:1.0E-15,PDD_41:1.0E-15,PDD_42:1.0E-15,PDD_43:1.0E-15,PDD_44:1.0E-15,PDD_45:1.0E-15,PDD_46:1.0E-15,PDD_47:1.0E-15,PDD_48:1.0E-15,PDD_49:1.0E-15,PDD_5:1.0E-15,PDD_50:1.0E-15,PDD_51:1.0E-15,PDD_52:1.0E-15,PDD_53:1.0E-15,PDD_54:1.0E-15,PDD_55:1.0E-15,PDD_56:1.0E-15,PDD_57:1.0E-15,PDD_58:1.0E-15,PDD_59:1.0E-15,PDD_6:1.0E-15,PDD_60:1.0E-15,PDD_61:1.0E-15,PDD_62:1.0E-15,PDD_63:1.0E-15,PDD_64:1.0E-15,PDD_65:1.0E-15,PDD_67:1.0E-15,PDD_68:1.0E-15,PDD_69:1.0E-15,PDD_7:1.0E-15,PDD_70:1.0E-15,PDD_71:1.0E-15,PDD_72:1.0E-15,PDD_73:1.0E-15,PDD_74:1.0E-15,PDD_75:1.0E-15,PDD_76:1.0E-15,PDD_77:1.0E-15,PDD_78:1.0E-15,PDD_79:1.0E-15,PDD_8:1.0E-15,PDD_80:1.0E-15,PDD_81:1.0E-15,PDD_82:1.0E-15,PDD_83:1.0E-15,PDD_84:1.0E-15,PDD_85:1.0E-15,PDD_86:1.0E-15,PDD_87:1.0E-15,PDD_88:1.0E-15,PDD_89:1.0E-15,PDD_9:1.0E-15,PDD_90:1.0E-15,PDD_91:1.0E-15,PDD_92:1.0E-15,PDD_93:1.0E-15,PDD_94:1.0E-15,PDD_95:1.0E-15,PDD_96:1.0E-15)C0:5.857486293396192E-4);";
		// STITree<Double> mapTree = SCF.getTree(t);
		// STITree<Double> imapTree = ClonalTreeDistance.getIntegerBLTree(mapTree);
//		String filename = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO8/Complete/allSample_trees.txt";
//		String opfilename = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO8/Complete/MPEAR_trees_geq6.txt";
////		String dir = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/set7/";
//		String dir = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/500_cells/100_sites/10_clones/";
////		String filename = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/100_cells/100_sites/10_clones/dataset2/best/sample_trees.txt";
//		String cellNames = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO8/CO8_cellnames.txt";
		
		
//		String filename = "/Volumes/My Passport/SiCloneFiT/realData/CO5/samples/all_tree_samples.txt";
//		String opfilename = "/Volumes/My Passport/SiCloneFiT/realData/CO5/samples/MPEAR_trees.txt";
//		String cellNames = "/Volumes/My Passport/SiCloneFiT/realData/CO5/CO5_cellnames.txt";
		
		String filename = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO5/samples/all_tree_samples.txt";
		String opfilename = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO5/samples/MPEAR_trees.txt";
		String cellNames = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO5/CO5_cellnames.txt";
		
//		String filename = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO5/best/best_trees.txt";
//		String cellNames = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/realData/CO5/CO5_cellnames.txt";
		
		
//		for (int j = 1; j <= 10; j++){
//			String dataset = "dataset" + j;
//			String filename = dir + dataset + "/30p_missing_samples/best/sample_trees.txt";
//			String outfile = dir + dataset + "/30p_missing_samples/best/nu_MAP_tree.txt";
//		
//		
//		
//		HashMap<String, Integer> treeMap = new HashMap<>();
//		HashMap<String, ArrayList<String>> treeListMap = new HashMap<>();
		SummarizeClonalTree SCT = new SummarizeClonalTree(filename, cellNames);
		SCT.writeValidTrees(filename, opfilename, SCF);
//		
//		try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
//			String line;
//			int i = 0;
//			while ((line = br.readLine()) != null) { 
//				String thisTreeS = line.trim();
//				String tr = SCT.getSameTreeSD(thisTreeS, treeListMap, treeMap, SCF);
////				STITree<Double> thisTree = SCF.getTree(thisTreeS);
////				STITree<Double> ithisTree =  ClonalTreeDistance.getIntegerBLTree(thisTree);
////				double dist = ClonalTreeDistance.getPairwiseCellSPTreeDist(imapTree, ithisTree, cellNames, nCell);
////				System.out.println(dist);
////				i += 1;
////				System.out.println(i);
//			}
//		}
//		
//		Map<String, Integer> sortedMap = SCT.sortByComparator(treeMap, false);
////		SCT.printMap(treeMap);
//		
////		int best = 0;
////		String bestTree = "";
//		ArrayList<Integer> values = new ArrayList<>();
//		ArrayList<String> keys = new ArrayList<>();
//		
//		for (String tree: sortedMap.keySet()){
//			
////			if (treeMap.get(tree) > best){
////				best = treeMap.get(tree);
////				bestTree = tree;
////			}
////			System.out.println(tree);
//			values.add(sortedMap.get(tree));
//			keys.add(tree);
////			System.out.println(sortedMap.get(tree));
//		}
//		
//		SCF.writeNewickTree(outfile, keys.get(0));
//		
////		for (int i = 0; i < 2; i++){
////			String tt = keys.get(i);
//////			for (String s:treeListMap.get(tt)){
//////				System.out.println(s);
//////			}
////			System.out.println(values.get(i));
////			System.out.println(keys.get(i));
////		}
//		}
		
//		System.out.println(best);
//		System.out.println(bestTree);

//		String tr = SCT.getSameTree(t, treeMap);
//		System.out.println(treeMap);
	}

}
