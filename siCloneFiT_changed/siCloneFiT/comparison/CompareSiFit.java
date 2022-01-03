/**
 * Dec 12, 2018
 */
package siCloneFiT.comparison;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

import cern.colt.Arrays;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import jeigen.DenseMatrix;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Dec 12, 2018
 */
public class CompareSiFit {
	public STITree<Double> phyloTree;
	public String newickTree;
	public static TNode[] leafset;
	
	/**
	 * Get the distance matrix of leaves
	 * @param filename
	 * @return
	 * @throws IOException
	 * Created On: Dec 12, 2018
	 */
	public static double[][] getLeafDistances(String filename) throws IOException{
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
		String treeN = SCF.readNewickString(filename);
		STITree<Double> t = SCF.getTree(treeN);
		Hashtable<TNode,Integer> idLookUp = new Hashtable<TNode,Integer>();
		int nLeaf = t.getLeafCount();
		TNode[] node_lookup = new TNode[nLeaf];
		leafset = node_lookup;
		for (int i = 1; i <= nLeaf; i++){
			String cellName = "sc" + Integer.toString(i);
			TNode cellNode = t.getNode(cellName);
			node_lookup[i-1] = cellNode;
			idLookUp.put(cellNode, idLookUp.size());
		}
		
		double[][] leafDistMatrix = getLeafDistanceMatrix(t, idLookUp, node_lookup);
//		for (TNode node: node_lookup){
//			System.out.println(node.getName());
//		}
		return leafDistMatrix;
	}
	
	public static void writeDistMatrixToFile (String outfile, double[][] dist){
		try{
			File file = new File(outfile);
		    FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
			for (double[] arr : dist){
				String d = "";
				for (int j = 0; j < arr.length; j++){
					d = d + Double.toString(arr[j]) + " ";
				}
				d = d + "\n";
//				writer.println(d);
				bw.write(d);
			}
			bw.close();
		} catch (IOException e){
			System.out.println("file not found\n");
		}
	}
	
	/**
	 * Reads a newick string from a file, construct the leaf distance matrix for the tree and write it to a file
	 * writes in the order sc1, sc2, sc3, ..., scM (M is the number of leaves)
	 * @param newickFile
	 * @param outfile
	 * @throws IOException
	 * Created On: Dec 12, 2018
	 */
	public static void constructLeafDistMat_WriteToFile(String newickFile, String outfile) throws IOException{
		double[][] leafDistMatrix = getLeafDistances(newickFile);
		writeDistMatrixToFile(outfile, leafDistMatrix);
	}
	
	/**
	 * Helper function for getting the leaf distance matrix
	 * @param t
	 * @param id_lookup
	 * @param node_lookup
	 * @return
	 * Created On: Dec 12, 2018
	 */
	public static double[][] getLeafDistanceMatrix(Tree t, Hashtable<TNode,Integer> id_lookup, TNode[] node_lookup){
		// compute distances
		double[][] distances = new double[id_lookup.size()][id_lookup.size()];

		for(int i = 0; i < distances.length; i++) {
		   findDistances(i, distances[i], id_lookup, node_lookup);
		}

		return distances;
	}
	
	public static final void findDistances(int leaf_idx, double[] distances, Hashtable<TNode,Integer> id_lookup, TNode[] node_lookup) {
		distances[leaf_idx] = 0;
		TNode leaf = node_lookup[leaf_idx];

		recurComputeDist(leaf, leaf.getParent(), leaf.getParentDistance(), distances, id_lookup, node_lookup);
	}
	
	public static final void recurComputeDist(TNode incoming, TNode curr_node, double curr_dist, double[] distances, Hashtable<TNode,Integer> id_lookup, TNode[] node_lookup) {

		// if we're a leaf, compute the distance
		if(curr_node.isLeaf()) {
			int idx = id_lookup.get(curr_node);
			distances[idx] = curr_dist;

			return;
		}

		// recur to parent if the incoming node wasn't this node's parent
		if(curr_node.getParent() != incoming && curr_node.getParent() != null) {
			recurComputeDist(curr_node, curr_node.getParent(), curr_dist + curr_node.getParentDistance(), distances, id_lookup, node_lookup);
		}

		// recur to all children
		for(TNode child : curr_node.getChildren()) {

			if(child != incoming) {
				recurComputeDist(curr_node, child, curr_dist + child.getParentDistance(), distances, id_lookup, node_lookup);
			}
		}

		return;
	}

	/**
	 * @param args
	 * Created On: Dec 12, 2018
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
//		String newickFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/fsTest/recurTest/recur_0.2/dataset3/noisy_genotype_dataset3_mlTree.newick";
//		double[][] leafDistMatrix = CompareSiFit.getLeafDistances(newickFile);
//		System.out.println(new DenseMatrix(leafDistMatrix));

//		String outfile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/fsTest/recurTest/recur_0.2/dataset3/SiFit_dataset3_mlTree_leafDistMat.txt";
//		CompareSiFit.constructLeafDistMat_WriteToFile(newickFile, outfile);
//		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
		
//		String s = "((sc10:0.027467401004195827,((((sc76:1.249215952497375E-4,sc55:0.0010570373743423091):0.0058904876804110155,((sc84:0.0029301439451800866,sc8:0.0295968018909013):0.016096222067897053,(((((sc89:0.003630635115361815,sc6:4.284098018885956E-4):0.008116264173390337,sc15:0.0032339076109172025):4.393673288718997E-4,sc19:0.0025655870559964243):0.00978435444592332,((sc82:0.008533172742128694,((sc29:0.012767493839882104,sc11:0.019084015217594816):3.276599598470233E-6,sc12:0.0028783706179290333):0.008716222411580427):0.011285982239940979,((sc98:0.003999564138054065,(((sc69:8.13201308388012E-5,(sc96:5.026205088069045E-4,sc74:0.0025355401476660974):0.00770769341850678):2.9476384794400384E-5,sc14:3.8310534241144134E-4):1.1078491047891932E-5,sc49:0.0169736725081162):1.6938030281533116E-5):1.092294779948755E-4,sc42:0.0011781573527790833):1.0861423094174162E-4):5.272462787596569E-4):6.805049738542432E-6,sc62:1.310991191135016E-4):0.01326571410805215):0.00738292261138432):0.011767397158257932,sc46:3.4405033898708615E-5):0.07620228840764258,((sc20:0.0018966169361840741,sc71:0.0075044948977361):0.007017489316798752,((((sc61:0.0583815717731851,sc1:0.01364517667108458):1.5790053505395535E-5,((sc39:2.8565354039559617E-5,sc37:1.9263661890699186E-5):3.731618534388358E-6,(sc36:7.405753228719585E-5,sc23:7.392565123206612E-5):3.409860185216348E-5):2.9012419585184152E-5):0.015310224626610214,sc38:0.01980637614453357):1.1449811609840083E-4,sc40:0.0058565979142355725):0.019002407246505138):0.06266520222289979):0.05735675681813484):0.022974077333240746,(((sc93:0.04330336015536107,(((sc56:4.327982944939254E-6,((sc88:2.0020724553491106E-5,(((((((sc79:1.536237940572877E-4,(sc41:0.025111989253391125,sc64:0.012203461410561714):5.92954387708117E-4):2.3705556598320292E-4,sc81:1.8848332987505336E-4):0.014031953800256191,sc87:5.8921323993867366E-5):4.182003884680622E-6,(sc27:1.19621836805313E-4,sc21:0.005648651768354047):9.48709104019083E-6):1.5752446492314228E-6,((((sc100:1.7651573835141085E-4,sc60:0.0016947462131785627):0.015276290017595446,sc86:0.004625508838923255):5.758150572718546E-5,sc30:0.0012497773702911251):0.012616465395673848,sc16:4.815661998344153E-5):2.382886918729038E-5):2.0674810390081416E-7,sc80:1.5074651640442232E-5):3.875204548397266E-6,sc83:2.5752156244899607E-5):1.0013031219621073E-5):1.514077385893131E-5,(sc35:3.536152587076541E-5,sc34:0.004207997941223263):3.261487513013786E-5):8.515441900580117E-6):6.018214095096179E-6,sc65:1.3870090678670023E-5):0.010207415604496279,(sc31:0.0020799343640165763,sc72:0.011047124939702628):0.016536301838907145):0.010696704427708545):0.01880399236372766,(sc4:0.13640977244895,((((((sc7:0.005648918866039127,sc78:6.806031956422588E-5):3.3714105877939096E-5,sc97:6.914631438483516E-6):1.6130870799744483E-5,((((((((sc24:0.02465338009206649,sc90:0.0011144374045033891):0.005687649745308865,sc3:3.651127145659954E-5):0.011537738083870501,sc17:0.03540759592700107):0.012351886253821344,sc51:9.097185534179054E-4):3.545744812213938E-4,((sc57:0.002778975345111575,sc25:0.050003549351663956):0.0026834470852951336,sc91:0.0016649440871584183):0.012228302816823028):0.014335312984449997,sc66:6.29155489598326E-4):0.030739679331533733,sc94:0.004387421683221365):0.1421522516761966,sc5:0.005741532658796962):0.1535844431186566):6.011301167372379E-6,sc47:0.006901840696794918):9.652370523838258E-5,sc22:8.349523935021784E-5):1.02659458651869E-5,sc73:0.001158389248866784):0.10357057311493836):0.035033957188866274):9.516496185744906E-5,(((sc77:0.0019263965555983207,sc18:7.858622177658093E-5):3.4836569666756544E-4,(((sc44:0.002824855134034152,sc59:0.009415351161043051):8.441938947876892E-4,sc52:0.0030973073242941507):0.010543155258590009,(((((((sc26:0.0017075736231171148,((sc43:0.01904809431243147,sc63:0.0293757345574881):0.001913516609304971,sc95:0.02468335375961346):0.009897893382537782):0.01524168527454413,sc70:7.57798014847358E-4):5.991592506442458E-5,sc45:0.02486653078109287):0.012812296795776407,((sc53:0.006160583251245497,sc54:4.173281230400297E-4):0.008145094050926172,sc92:4.947113532160188E-5):5.368767466974046E-5):1.6062009563703263E-4,(((sc28:0.002555955154557731,sc32:0.0011526516717716621):0.0017259294204320331,sc33:0.004325676770496154):0.00956340977106434,sc13:0.055801771832687766):1.6489986669027937E-4):2.7600909500636624E-4,((sc68:2.1383188969048903E-4,((sc67:0.0017606899793517576,sc48:0.005686184351662954):5.666006140207195E-5,sc2:0.03341668094895808):0.006035104477817244):1.1868756283297274E-4,sc58:4.894763116729639E-4):0.009079280328230849):2.583804183329128E-4,sc85:4.4576151093705637E-4):0.02185101646863895):0.009260439389459972):0.18577695776146916,(sc75:5.689944825864688E-4,((sc50:0.0035370248647416877,sc99:0.01257742156293215):3.548831227946524E-4,sc9:0.0020363136014595108):0.031049433366252353):0.17596176184834417):0.021166198148471465):0.02629373696379806);";
//		STITree<Double> ts = SCF.getTree(s);
//		System.out.println(Arrays.toString(ts.getLeaves()));
		
		String dir = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/noDoublet/fnRateExpt/fn4/SiFit_n_SCITE/";
//		String dataset = "dataset" + 1;
//		String newickFile = dir + dataset + "/SiFit_dataset1_mutation_tree.newick";
//		String outFile = dir + dataset + "/SiFit_" + dataset + "_mutationTree_leafDistMat.txt";
//		CompareSiFit.constructLeafDistMat_WriteToFile(newickFile, outFile);
		
//		for (int i = 1; i <= 5; i++){
//			String dataset = "dataset" + i;
//			String newickFile = dir + dataset + "/noisy_genotype_" + dataset + "_mlTree.newick";
//			String outFile = dir + dataset + "/SiFit_" + dataset + "_mlTree_leafDistMat.txt";
//			CompareSiFit.constructLeafDistMat_WriteToFile(newickFile, outFile);
//		}
		
		// Leaf distance for the mutation tree
		for (int i = 1; i <= 10; i++){
			String dataset = "dataset" + i;
			String newickFile = dir + dataset + "/SiFit_" + dataset + "_mutation_tree.newick";
			String outFile = dir + dataset + "/SiFit_" + dataset + "_mutationTree_leafDistMat.txt";
			CompareSiFit.constructLeafDistMat_WriteToFile(newickFile, outFile);
		}
	}

}
