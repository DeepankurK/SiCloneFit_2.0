/**
 * Oct 30, 2017
 */
package siCloneFiT.main;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;

import SiFit.metric.CompareTrees;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import siCloneFiT.algorithm.GibbsSamplingDoublet;
import siCloneFiT.io.GenotypeMatrixReader;
import siCloneFiT.metric.ClonalTreeDistance;
import siCloneFiT.objects.DoubletPosteriorSample;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Oct 30, 2017
 */
public class RunSiCloneFiTDoublet {
	public static int nCell;							// Number of cells
	public static int nMut;								// Number of mutations
	public static double _fp = 0.0174;					// FP rate, alpha
	public static double _fnStart = 0.1375;				// FN rate, beta
	public static double _doublet = 0.1;				// doublet rate, delta
	public static double _delProb = 0.05;				// Deletion Rate
	public static double _LOHProb = 0.01;				// LOH Rate
	public static double _recurProb = 0.05;				// Recurrent mutation Rate
	public static double _missing = 0.0;
	public static int restart = 1;						// Number of restarts
	public static double GammaA = 1;					// Hyper-parameter for alpha_0
	public static double GammaB = 1;					// Hyper-parameter for alpha_0
	public static int iterT = 50000;					// Total number of iterations
	public static int iterP = 1000;						// Number of iterations after which likelihood will be printed
	public static int burnIn = 100;
	public static int sampleTreeIter = 200;
	public static int dataFlag = 0;						// Flag to indicate the data type, binary or ternary
	public static int modelFlag = 1;					// Flag to indicate which model to use
	public static String varMatFilename = null;			// Filename of the input matrix
	public static String trueTreeFilename = null;		// File containing the true tree in newick form
	public static String cellNamesFilename = null;		// File containing the cell names
	public static String workingDir = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/testData/100_cells/dataset8/";

	/**
	 * @param args
	 * Created On: Oct 30, 2017
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		readArguments(args);
//		nMut = 30;
//		nCell = 100;
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
//		varMatFilename = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Clonal_Phylogeny_SCS/testData/100_cells/dataset8/noisy_genotype_dataset8.txt";
		// Read the observed genotype matrix
		GenotypeMatrixReader GMR = new GenotypeMatrixReader(varMatFilename);
		HashMap<String, Integer[]> cellGTVectorMap = GMR.getCellGTVectorMap(varMatFilename, cellNamesFilename, nMut, nCell);
		ArrayList<String> singleCellNames = GMR.scNames;
		
		// Construct a list of single cell objects
		ArrayList<SingleCell> listSingleCells = SCF.constructSingleCellList(cellGTVectorMap, singleCellNames, nCell);
		HashMap<Integer, Integer[]> cellIDGTVectorMap = new HashMap<>(); // HashMap<cellID, observedVector>
		for (SingleCell S : listSingleCells){
			cellIDGTVectorMap.put(S.cellID, S.observedGTVector);
		}
		
		// Concentration parameter
		double alpha_0 = 1; 
		GammaDistribution alpha_0PriorDist = new GammaDistribution(GammaA, GammaB);

		// Work on the parameters
		// Distribution for Fn
        double fnPriorMean = _fnStart;
        double fnPriorSD = fnPriorMean*0.5;
        double betaPriora = ((1 - fnPriorMean)*fnPriorMean*fnPriorMean/(fnPriorSD*fnPriorSD)) - fnPriorMean;
        double betaPriorb = betaPriora * ((1/fnPriorMean) - 1);
        BetaDistribution fnPriorDist = new BetaDistribution(betaPriora, betaPriorb);

        // Distribution for Fp
        double fpPriorMean = _fp;
        double fpPriorSD = fpPriorMean * 0.5;
        double fpBetaPriora = ((1 - fpPriorMean)*fpPriorMean*fpPriorMean/(fpPriorSD*fpPriorSD)) - fpPriorMean;
        double fpBetaPriorb = fpBetaPriora * ((1/fpPriorMean) - 1);
        BetaDistribution fpPriorDist = new BetaDistribution(fpBetaPriora, fpBetaPriorb);
        
        // Distribution for delta (doublet rate)
        double _doubletPriora = SCF.getBetaPriora(_doublet, _doublet*0.5);
        double _doubletPriorb = _doubletPriora * ((1/_doublet) - 1);
        BetaDistribution doubletPriorDist = new BetaDistribution(_doubletPriora, _doubletPriorb);
        

        // Parameterize prior distribution for _delProb
        double _delProbPriora = SCF.getBetaPriora(_delProb, _delProb*0.1);
        double _delProbPriorb = _delProbPriora * ((1/_delProb) - 1);
        BetaDistribution delProbPriorDist = new BetaDistribution(_delProbPriora, _delProbPriorb); 
        
        // Parameterize prior distribution for _LOHProb
        double _LOHProbPriora = SCF.getBetaPriora(_LOHProb, _LOHProb*0.1);
        double _LOHProbPriorb = _LOHProbPriora * ((1/_LOHProb) - 1);
        BetaDistribution LOHProbPriorDist = new BetaDistribution(_LOHProbPriora, _LOHProbPriorb);
        
        // Parameterize prior distribution for recurProb
        double recurProbPriora = SCF.getBetaPriora(_recurProb, _recurProb*0.3); // More variance for recurProb
        double recurProbPriorb = recurProbPriora * ((1/_recurProb) - 1);
        BetaDistribution recurProbPriorDist = new BetaDistribution(recurProbPriora, recurProbPriorb);
        
        // List of posterior samples from each restart
        ArrayList<DoubletPosteriorSample> bestPosteriorSamples = new ArrayList<>();
        ArrayList<Double> listPosterior = new ArrayList<>();

        // Create the directory for storing samples
        String samplesDir;
        if (_missing > 0.0){
        	int percentMissing = (int) (_missing*100);
        	samplesDir = workingDir + Integer.toString(percentMissing) + "p_missing_samples/";
        }
        else
        	samplesDir = workingDir + "samples/";
		File samplesDirFP = new File(samplesDir);
		samplesDirFP.mkdir();
		
		
        // Looping over the number of restarts
        for (int i = 0; i < restart; i++){
        	System.out.printf("Run %d\n", i);
        	DoubletPosteriorSample currRunMAPSample = GibbsSamplingDoublet.samplePosterior(nCell, nMut, listSingleCells, singleCellNames, cellIDGTVectorMap, _fnStart, fnPriorDist, _fp, fpPriorDist, _doublet, doubletPriorDist, _delProb, delProbPriorDist, _LOHProb, LOHProbPriorDist, _recurProb, recurProbPriorDist, alpha_0, GammaA, GammaB, alpha_0PriorDist, burnIn, iterT, iterP, sampleTreeIter, dataFlag, i, samplesDir, SCF);
        	bestPosteriorSamples.add(currRunMAPSample);
        	listPosterior.add(currRunMAPSample.logPosteriorProb);
        }
        
        double bestPosterior = listPosterior.get(0);
		int bestConfigIndex = 0;
		for (int i = 1; i < restart; i++){
			if (listPosterior.get(i) > bestPosterior){
				bestConfigIndex = i;
				bestPosterior = listPosterior.get(i);
			}			
		}
		
		// Create the directory for storing best samples
		String bestSamplesDir = samplesDir + "best/";
		File bestSamplesDirFP = new File(bestSamplesDir);
		bestSamplesDirFP.mkdir();
		
		/*
		 * Write the best samples to the best directory
		 */
		copyBestSamples(samplesDir, bestSamplesDir, bestConfigIndex, bestPosteriorSamples.get(bestConfigIndex), SCF);
		
		// If true Tree is provided, print tree reconstruction error
		if (trueTreeFilename != null){
			String trueTreeNewick = SCF.readNewickString(trueTreeFilename);
			STITree<Double> trueTree = SCF.getTree(trueTreeNewick);
			STITree<Double> inferredTree = SCF.getTree(bestPosteriorSamples.get(bestConfigIndex).cloneTree.toNewick());
			CompareTrees CTObj = new CompareTrees(trueTree);
			double FPdist = CTObj.calcFPDist(inferredTree);
			double FNdist = CTObj.calcFNDist(inferredTree);
			System.out.printf("FP error = %f\n", FPdist);
			System.out.printf("FN error = %f\n", FNdist);
			double clonalDist = ClonalTreeDistance.getPairwiseCellSPTreeDist(trueTreeFilename, inferredTree, nCell, SCF);
			System.out.printf("shortest path distance = %f\n", clonalDist);
		}
		
		
	}
	
	/**
	 * Write the best samples to the best directory
	 * @param samplesDir
	 * @param bestSamplesDir
	 * @param bestConfigIndex
	 * @param bestMAPSample
	 * @param SCF
	 * @throws IOException
	 * Created On: Dec 21, 2017
	 */
	private static void copyBestSamples(String samplesDir, String bestSamplesDir, int bestConfigIndex, DoubletPosteriorSample bestMAPSample, SCFUtilityFunctions SCF) throws IOException{
		// Determine the source files
		String clonalClusterFile = samplesDir + "samples_" + bestConfigIndex + "_clusters_1.txt";
		String clonalCluster2File = samplesDir + "samples_" + bestConfigIndex + "_clusters_2.txt";
        String clonalGenotypesFile = samplesDir + "samples_" + bestConfigIndex + "_genotypes.txt";
        String clonalTreesFile = samplesDir + "samples_" + bestConfigIndex + "_trees.txt";
        String cellDoubleFlagFile = samplesDir + "samples_" + bestConfigIndex + "_doubletFlags.txt";
        String posteriorErrorRatesFile = samplesDir + "samples_" + bestConfigIndex + "_posteriors.txt";
        String predictedFile = samplesDir + "sample_" + bestConfigIndex + "_MAP_predicted_genotype.txt";
        
        // File pointers for source files 
        File clonalClusterSrc = new File(clonalClusterFile);
        File clonalCluster2Src = new File(clonalCluster2File);
        File clonalGenotypesSrc = new File(clonalGenotypesFile);
        File clonalTreesSrc = new File(clonalTreesFile);
        File cellDoubletFlagSrc = new File(cellDoubleFlagFile);
        File posteriorErrorRatesSrc = new File(posteriorErrorRatesFile);
        File predictedSrc = new File(predictedFile);
        
        // Destinate files path
        String bestclonalClusterFile = bestSamplesDir + "best_" + "clusters_1.txt";
        String bestclonalCluster2File = bestSamplesDir + "best_" + "clusters_2.txt";
        String bestclonalGenotypesFile = bestSamplesDir + "best_" + "genotypes.txt";
        String bestclonalTreesFile = bestSamplesDir + "best_" + "trees.txt";
        String bestcellDoubleFlagFile = bestSamplesDir + "best_" + "doubletFlags.txt";
        String bestposteriorErrorRatesFile = bestSamplesDir + "best_" + "posteriors.txt";
        String bestpredictedFile = bestSamplesDir + "best_" + "MAP_predicted_genotype.txt";
        
        // File pointers for destination files
        File clonalClusterDest = new File(bestclonalClusterFile);
        File clonalCluster2Dest = new File(bestclonalCluster2File);
        File clonalGenotypesDest = new File(bestclonalGenotypesFile);
        File clonalTreesDest = new File(bestclonalTreesFile);
        File cellDoubletFlagDest = new File(bestcellDoubleFlagFile);
        File posteriorErrorRatesDest = new File(bestposteriorErrorRatesFile);
        File predictedDest = new File(bestpredictedFile);
        
        // Copy from source to destination
        copyFileUsingApacheCommonsIO(clonalClusterSrc, clonalClusterDest);
        copyFileUsingApacheCommonsIO(clonalCluster2Src, clonalCluster2Dest);
        copyFileUsingApacheCommonsIO(clonalGenotypesSrc, clonalGenotypesDest);
        copyFileUsingApacheCommonsIO(clonalTreesSrc, clonalTreesDest);
        copyFileUsingApacheCommonsIO(cellDoubletFlagSrc, cellDoubletFlagDest);
        copyFileUsingApacheCommonsIO(posteriorErrorRatesSrc, posteriorErrorRatesDest);
        copyFileUsingApacheCommonsIO(predictedSrc, predictedDest);
        
        // Best MAP solutions
        String bestMAPTreeFile = bestSamplesDir + "best_MAP_tree.txt";
        String bestMAPcellCloneFile = bestSamplesDir + "best_MAP_clusters_1.txt";
        String bestMAPcellClone2File = bestSamplesDir + "best_MAP_clusters_2.txt";
        String bestMAPcellDoubletFlagFile = bestSamplesDir + "best_MAP_doubletFlags.txt";
        SCF.writeNewickTree(bestMAPTreeFile, bestMAPSample.cloneTree.toNewick());
        
        ArrayList<String> cellCloneID = new ArrayList<>();
        for (int i : bestMAPSample.cellCloneIDVector)
        	cellCloneID.add(Integer.toString(i));
        String bestMAPcellCloneStr = String.join(" ", cellCloneID); 
        SCF.writeNewickTree(bestMAPcellCloneFile, bestMAPcellCloneStr);
        
        String cellClone2Str = Integer.toString(bestMAPSample.cell2ndCloneIDVector.get(0));
        String cellDoubletFlagStr = Integer.toString(bestMAPSample.cellDoubletFlagArr[0]);
        for (int i = 1; i < bestMAPSample.cell2ndCloneIDVector.size(); i++){
        	cellClone2Str += " " + Integer.toString(bestMAPSample.cell2ndCloneIDVector.get(i));
        	cellDoubletFlagStr += " " + Integer.toString(bestMAPSample.cellDoubletFlagArr[i]);
        }
        SCF.writeNewickTree(bestMAPcellClone2File, cellClone2Str);
        SCF.writeNewickTree(bestMAPcellDoubletFlagFile, cellDoubletFlagStr);     

	}
	
	/**
	 * Read input arguments
	 * @param args
	 * Created On: Dec 21, 2017
	 */
	public static void readArguments(String[] args){
		int nPar = args.length;
		for (int i = 0; i < nPar; i = i+2){
			if (args[i].equals("-m") == true){
				nCell = Integer.parseInt(args[i+1]); 
			}
			else if (args[i].equals("-n") == true){
				nMut = Integer.parseInt(args[i+1]); 
			}
			else if (args[i].equals("-fp") == true){
				_fp = Double.parseDouble(args[i+1]);
				if (_fp < 0 || _fp > 1)
					throw new IllegalArgumentException("Invalid input for fp rate, should be between 0 and 1");
			}
			else if (args[i].equals("-fn") == true){
				_fnStart = Double.parseDouble(args[i+1]);
				if (_fnStart < 0 || _fnStart > 1)
					throw new IllegalArgumentException("Invalid input for fn rate, should be between 0 and 1");
			}
			else if (args[i].equals("-doublet") == true){
				_doublet = Double.parseDouble(args[i+1]);
				if (_doublet < 0 || _doublet > 1)
					throw new IllegalArgumentException("Invalid input for doublet rate, should be between 0 and 1");
			}
			else if (args[i].equals("-burnin") == true){
				burnIn = Integer.parseInt(args[i+1]);
			}
			else if (args[i].equals("-iter") == true){
				iterT = Integer.parseInt(args[i+1]);
			}
			else if (args[i].equals("-printIter") == true){
				iterP = Integer.parseInt(args[i+1]);
			}
			else if (args[i].equals("-treeIter") == true){
				sampleTreeIter = Integer.parseInt(args[i+1]);
			}
			else if (args[i].equals("-r") == true){
				restart = Integer.parseInt(args[i+1]);
			}
			else if (args[i].equals("-df") == true){
				dataFlag = Integer.parseInt(args[i+1]); 
				if (dataFlag > 1)
					throw new IllegalArgumentException("Invalid input for dataFlag parameter, should be 0 or 1");
			}
			else if (args[i].equals("-ipMat") == true){
				varMatFilename = args[i+1];
			}
			else if (args[i].equals("-missing") == true){
				_missing = Double.parseDouble(args[i+1]);
			}
			else if (args[i].equals("-trueTree") == true){
				trueTreeFilename = args[i+1];
			}
			else if (args[i].equals("-cellNames") == true){
				cellNamesFilename = args[i+1];
			}
			else if (args[i].equals("-outDir") == true){
				workingDir = args[i+1];
				if (!workingDir.endsWith("/"))
					workingDir += "/";
			}
		}
	}
	
	/**
	 * Copy file
	 * @param source
	 * @param dest
	 * @throws IOException
	 * Created On: Dec 21, 2017
	 */
	private static void copyFileUsingApacheCommonsIO(File source, File dest) throws IOException {
	    FileUtils.copyFile(source, dest);
	}

}
