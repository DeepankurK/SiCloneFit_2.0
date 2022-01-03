/**
 * Oct 26, 2017
 */
package siCloneFiT.main;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter; 
import java.io.UnsupportedEncodingException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;

import SiFit.metric.CompareTrees;
import cern.colt.Arrays;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import siCloneFiT.algorithm.GibbsSamplingSinglet;
import siCloneFiT.io.GenotypeMatrixReader;
import siCloneFiT.metric.ClonalTreeDistance;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.objects.SingletPosteriorSample;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Oct 26, 2017
 */
public class RunSiCloneFiTSinglet {

	public static int nCell;							// Number of cells
	public static int nMut;								// Number of mutations
	public static double _fp = 0.05;					// FP rate, alpha
	public static double _fnStart = 0.2;				// FN rate, beta
	public static double _delProb = 0.1;				// Deletion Rate
	public static double _LOHProb = 0.1;				// LOH Rate
	public static double _recurProb = 0.05;				// Recurrent mutation Rate
	public static double _missing = 0.0;
	public static double GammaA = 5;					// Hyper-parameter for alpha_0
	public static double GammaB = 5;					// Hyper-parameter for alpha_0
	public static int restart = 3;						// Number of restarts
	public static int burnIn = 100;
	public static int iterT = 500;					    // Total number of iterations
	public static int iterP = 1000;						// Number of iterations after which likelihood will be printed
	public static int sampleTreeIter = 200;
	public static int dataFlag = 0;						// Flag to indicate the data type, binary or ternary
	public static int modelFlag = 1;					// Flag to indicate which model to use
	public static String allele_MatFilename = null;		// Filename of the Allele input Matrix
	public static String ref_MatFilename = null;		// Filename of the Total input Matrix
	public static String trueTreeFilename = null;		// File containing the true tree in newick form
	public static String cellNamesFilename = null;		// File containing the cell names
	public static String workingDir = "/home/deepank/Downloads/Prof_Hamim/Data/";

	/**
	 * @param args
	 * Created On: Oct 26, 2017
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		readArguments(args);
//		nMut = 50;
//		nCell = 100;
		SCFUtilityFunctions SCF = new SCFUtilityFunctions();
//		varMatFilename = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/data_SiCloneFiT/test/dataset6/noisy_genotype_dataset6.txt";
		
		// Read the observed genotype matrix
		GenotypeMatrixReader GMR = new GenotypeMatrixReader(allele_MatFilename);
		HashMap<String, Integer[]> cellGTVectorMap_Mut = GMR.getCellGTVectorMap(allele_MatFilename, cellNamesFilename, nMut, nCell);
		ArrayList<String> singleCellNames = GMR.scNames;
		GenotypeMatrixReader GMR2 = new GenotypeMatrixReader(ref_MatFilename);
		HashMap<String, Integer[]> cellGTVectorMap_Ref = GMR2.getCellGTVectorMap(ref_MatFilename, cellNamesFilename, nMut, nCell);
		
		// Construct a list of single cell objects
		ArrayList<SingleCell> listSingleCells = SCF.constructSingleCellList(cellGTVectorMap_Mut, singleCellNames, nCell);
		HashMap<Integer, Integer[]> cellIDGTVectorMap = new HashMap<>(); // HashMap<cellID, observedVector>
		ArrayList<Integer> cellList = new ArrayList<>();
		for (SingleCell S : listSingleCells){
			cellIDGTVectorMap.put(S.cellID, S.observedGTVector);
			cellList.add(S.cellID);
		}
		ArrayList<SingleCell> listSingleCells_tot = SCF.constructSingleCellList(cellGTVectorMap_Ref, singleCellNames, nCell);
		HashMap<Integer, Integer[]> cellIDRefVectorMap = new HashMap<>(); // HashMap<cellID, observedVector>
		for (SingleCell S : listSingleCells_tot){
			cellIDRefVectorMap.put(S.cellID, S.observedGTVector);
		}



		
		// Concentration parameter
		double alpha_0 = 10; 
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
        ArrayList<SingletPosteriorSample> bestPosteriorSamples = new ArrayList<>();
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
		
//		ArrayList<Integer[]> test = new ArrayList<>();
//		test.add(new Integer[]{1,1,2,2});
//		test.add(new Integer[]{3,3,9,1});
//		
//		String newFilename = samplesDir + "test.txt";
//		PrintWriter writer = new PrintWriter(newFilename, "UTF-8");
//		writer.println(test);
////		writer.println("hello again");
//		writer.close();
		
//		SingletPosteriorSample[] MAPSamples = new SingletPosteriorSample[restart];
//		ArrayList<ArrayList<SingleCell>> copyListSingleCellsList = new ArrayList<>();
//		ArrayList<Double> alpha_0_list = new ArrayList<>();
//		ArrayList<Double> delProb_list = new ArrayList<>();
//		ArrayList<Double> LOHProb_list = new ArrayList<>();
//		ArrayList<Double> recurProb_list = new ArrayList<>();
//		
//		for (int i = 0; i < restart; i++){
//			ArrayList<SingleCell> copyListSingleCells = SCF.copyCellList(listSingleCells, nMut);
//			copyListSingleCellsList.add(copyListSingleCells);
//			alpha_0_list.add(alpha_0);
//			delProb_list.add(_delProb);
//			LOHProb_list.add(_LOHProb);
//			recurProb_list.add(_recurProb);
//		}
//		
//		// Launching HabaneroApp for parallel computation of likelihood
//		launchHabaneroApp(() -> {
//        
//        // Looping over the number of restarts
//			forallChunked(0, restart, (i) -> {
//				try {
//					MAPSamples[i] = GibbsSamplingSinglet.samplePosterior(nCell, nMut, copyListSingleCellsList.get(i), singleCellNames, cellIDGTVectorMap, fnPriorMean, fnPriorDist, fpPriorMean, fpPriorDist, delProb_list.get(i), delProbPriorDist, LOHProb_list.get(i), LOHProbPriorDist, recurProb_list.get(i), recurProbPriorDist, alpha_0_list.get(i), GammaA, GammaB, alpha_0PriorDist, burnIn, iterT, iterP, dataFlag, i, samplesDir, SCF);
//				} catch (FileNotFoundException e) {
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				} catch (UnsupportedEncodingException e) {
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				}
//			});
        for (int i = 0; i < restart; i++){
        	System.out.printf("Run %d\n", i);
        	
        	SingletPosteriorSample currRunMAPSample = GibbsSamplingSinglet.samplePosterior(nCell, nMut, listSingleCells,listSingleCells_tot, singleCellNames, cellIDGTVectorMap,cellIDRefVectorMap, fnPriorMean, fnPriorDist, fpPriorMean, fpPriorDist, _delProb, delProbPriorDist, _LOHProb, LOHProbPriorDist, _recurProb, recurProbPriorDist, alpha_0, GammaA, GammaB, alpha_0PriorDist, burnIn, iterT, iterP, sampleTreeIter, dataFlag, i, samplesDir, SCF);
        	bestPosteriorSamples.add(currRunMAPSample);
        	listPosterior.add(currRunMAPSample.logPosteriorProb);
        }
        
//		});
//		
//		for (SingletPosteriorSample currRunMAPSample: MAPSamples){
//			bestPosteriorSamples.add(currRunMAPSample);
//			System.out.println(currRunMAPSample.sampleFn);
//			listPosterior.add(currRunMAPSample.logPosteriorProb);
//		}
        
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
		System.out.printf("MAP delProb = %f\n", bestPosteriorSamples.get(bestConfigIndex).sampleDelProb);
		System.out.printf("MAP LOHProb = %f\n", bestPosteriorSamples.get(bestConfigIndex).sampleLOHProb);
		
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
	 * Created On: Oct 30, 2017
	 */
	private static void copyBestSamples(String samplesDir, String bestSamplesDir, int bestConfigIndex, SingletPosteriorSample bestMAPSample, SCFUtilityFunctions SCF) throws IOException{
		String clonalClusterFile = samplesDir + "samples_" + bestConfigIndex + "_clusters.txt";
        String clonalGenotypesFile = samplesDir + "samples_" + bestConfigIndex + "_genotypes.txt";
        String clonalTreesFile = samplesDir + "samples_" + bestConfigIndex + "_trees.txt";
        String posteriorErrorRatesFile = samplesDir + "samples_" + bestConfigIndex + "_posteriors.txt";
        String predictedFile = samplesDir + "sample_" + bestConfigIndex + "_MAP_predicted_genotype.txt";
        
        File clonalClusterSrc = new File(clonalClusterFile);
        File clonalGenotypesSrc = new File(clonalGenotypesFile);
        File clonalTreesSrc = new File(clonalTreesFile);
        File posteriorErrorRatesSrc = new File(posteriorErrorRatesFile);
        File predictedSrc = new File(predictedFile);
        
        String bestclonalClusterFile = bestSamplesDir + "best_" + "clusters.txt";
        String bestclonalGenotypesFile = bestSamplesDir + "best_" + "genotypes.txt";
        String bestclonalTreesFile = bestSamplesDir + "best_" + "trees.txt";
        String bestposteriorErrorRatesFile = bestSamplesDir + "best_" + "posteriors.txt";
        String bestpredictedFile = bestSamplesDir + "best_" + "MAP_predicted_genotype.txt";
        
        File clonalClusterDest = new File(bestclonalClusterFile);
        File clonalGenotypesDest = new File(bestclonalGenotypesFile);
        File clonalTreesDest = new File(bestclonalTreesFile);
        File posteriorErrorRatesDest = new File(bestposteriorErrorRatesFile);
        File predictedDest = new File(bestpredictedFile);
        
        copyFileUsingApacheCommonsIO(clonalClusterSrc, clonalClusterDest);
        copyFileUsingApacheCommonsIO(clonalGenotypesSrc, clonalGenotypesDest);
        copyFileUsingApacheCommonsIO(clonalTreesSrc, clonalTreesDest);
        copyFileUsingApacheCommonsIO(posteriorErrorRatesSrc, posteriorErrorRatesDest);
        copyFileUsingApacheCommonsIO(predictedSrc, predictedDest);
        
        String bestMAPTreeFile = bestSamplesDir + "best_MAP_tree.txt";
        String bestMAPcellCloneFile = bestSamplesDir + "best_MAP_clusters.txt";
        SCF.writeNewickTree(bestMAPTreeFile, bestMAPSample.cloneTree.toNewick());
        
        ArrayList<String> cellCloneID = new ArrayList<>();
        for (int i : bestMAPSample.cellCloneIDVector)
        	cellCloneID.add(Integer.toString(i));
        String bestMAPcellCloneStr = String.join(" ", cellCloneID); 
        SCF.writeNewickTree(bestMAPcellCloneFile, bestMAPcellCloneStr);
	}
	
	/**
	 * Read the input arguments and populate the 
	 * public variables to override the default values
	 * @param args
	 * Created On: Nov 1, 2017
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
			else if (args[i].equals("-allMat") == true){
				allele_MatFilename = args[i+1];
			}
			else if (args[i].equals("-refMat") == true){
				ref_MatFilename = args[i+1];
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
	 * Created On: Oct 30, 2017
	 */
	private static void copyFileUsingApacheCommonsIO(File source, File dest) throws IOException {
	    FileUtils.copyFile(source, dest);
	}

}
