/**
 * Aug 21, 2017
 */
package algorithm;

import java.util.ArrayList;
//import java.util.Random;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

//import siCloneFiT.likelihood.DoubletIndicatorLikelihood;
import siCloneFiT.likelihood.SCSAmplificationErrorLikelihood;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.SingleCell;
//import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * @author hz22
 * Aug 21, 2017
 */
public class SampleErrorRates {
	
	//private static Random _rng = new Random();
	
	/**********************************************************************
	 * Functions to sample new value of Fn, considering the singlet model *
	 **********************************************************************/
	
	/*
	 * Functions to sample Fn for the heuristic search
	 */
	/**
	 * Sample new Fn value based on the current Fn value (normal distribution)
	 * @param currFn
	 * @param fnSD
	 * @return
	 * Created On: Oct 15, 2017
	 */
	public static double sampleNewErrorRateHS(double currFn, double fnSD){
		NormalDistribution FnNormalDist = new NormalDistribution(currFn, fnSD);
		double newFn = FnNormalDist.sample();
		if (newFn < 0)
			newFn = Math.abs(newFn);
		if (newFn > 1)
			newFn = newFn - 2*(newFn - 1);
		return newFn;
	}
	
	public static ArrayList<Double> sampleNewFnHS(ArrayList<SingleCell> listCells,ArrayList<SingleCell> listCells_tot, ArrayList<Clone> listClone,
			 SCSAmplificationErrorLikelihood ampErrLikelihoodObj,
			 double currFp, BetaDistribution fnPriorDist, int nMut){

		// Obtain the envelope (Y_max) of the posterior distribution of Fn
		// Uniform distributions for doing rejection sampling
		UniformRealDistribution fn_fX = new UniformRealDistribution(0.01, 0.99);         // Uniform distribution x ~ [x_min, x_max]
		UniformRealDistribution fn_gY = new UniformRealDistribution(0, 1);  // Uniform distribution ~ [0, Y_max] 
		double x;  // Sample x ~ U(x_min, x_max)
		// calculate f(x), here f is conditional posterior
		double x_posterior;
		double y ;// Sample y ~ U(0, Y_max)
		ArrayList<Double> newFn=new ArrayList<>();
//		int nSamples = 100;
//		ArrayList<Double> samples = new ArrayList<>();
//		ArrayList<Double> probabilities = new ArrayList<>();
//		int accepted = 0;
		for(int i=0;i<nMut;i++){
			double FnPosteriorMax = getFnPosteriorMax(listCells, listCells_tot, listClone, ampErrLikelihoodObj, fnPriorDist, currFp, nMut, i);
			while(true){
				x = fn_fX.sample();
				y = Math.log(fn_gY.sample()) + FnPosteriorMax;
				x_posterior = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixIndividualMutation(listCells,listCells_tot, listClone, currFp, x, nMut,i) + Math.log(fnPriorDist.density(x));
				// Return the sample x if y <= f(x)
				if (y <= x_posterior){
					newFn.add(x);
					break;
				}
	//			//// Repeat the sampling
	//			else{
	//				x = fn_fX.sample();
	//				x_posterior = Math.exp(ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listCells, listClone, currFp, x, nMut)) * fnPriorDist.density(x);
	//				y = fn_gY.sample();
	//			}
		}
		}
		return newFn;
//		double bestSample = samples.get(0);
//		double bestProb = probabilities.get(0);
//		for (int i = 1; i < nSamples; i++){
//			if (probabilities.get(i) > bestProb){
//				bestSample = samples.get(i);
//				bestProb = probabilities.get(i);
//			}
//		}
		
	}
	
	/**
	 * Sample a new value of Fn from the conditional posterior distribution
	 * Implements a Rejection Sampler
	 * @param listCells
	 * @param listClone
	 * @param ampErrLikelihoodObj
	 * @param currFp
	 * @param fnPriorDist
	 * @param nMut
	 * @return
	 * Created On: Aug 21, 2017
	 */
	/*
	public static double sampleNewFn(ArrayList<SingleCell> listCells,ArrayList<SingleCell> listCells_tot, ArrayList<Clone> listClone,
									 SCSAmplificationErrorLikelihood ampErrLikelihoodObj,
									 double currFp, BetaDistribution fnPriorDist, int nMut){
		
		// Obtain the envelope (Y_max) of the posterior distribution of Fn
		double FnPosteriorMax = getFnPosteriorMax(listCells,listCells_tot, listClone, ampErrLikelihoodObj, fnPriorDist, currFp, nMut);
		System.out.println("FnMax = "+ FnPosteriorMax);
		// Uniform distributions for doing rejection sampling
		UniformRealDistribution fn_fX = new UniformRealDistribution(0.01, 0.99);         // Uniform distribution x ~ [x_min, x_max]
		UniformRealDistribution fn_gY = new UniformRealDistribution(0, Math.exp(FnPosteriorMax));  // Uniform distribution ~ [0, Y_max] 
		double x = fn_fX.sample();  // Sample x ~ U(x_min, x_max)
		// calculate f(x), here f is conditional posterior
		double x_posterior = Math.exp(ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listCells, listCells_tot, listClone, currFp, x, nMut)) * fnPriorDist.density(x);
		double y = fn_gY.sample(); // Sample y ~ U(0, Y_max)
		double newFn;
		
		int nSamples = 100;
		ArrayList<Double> samples = new ArrayList<>();
		ArrayList<Double> probabilities = new ArrayList<>();
		int accepted = 0;
		while(accepted < nSamples){
			x = fn_fX.sample();
			y = fn_gY.sample();
			double likelihood = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listCells, listClone, currFp, x, nMut);
			double prior = Math.log(fnPriorDist.density(x));
			
			
			x_posterior = likelihood + prior; 
			// Return the sample x if y <= f(x)
			if (y <= Math.exp(x_posterior)){
//				System.out.println("likelihood = "+ likelihood);
//				System.out.println("prior = " + prior);
				samples.add(x);
				probabilities.add(x_posterior);
				accepted += 1;
			}
//			// Repeat the sampling
//			else{
//				x = fn_fX.sample();
//				x_posterior = Math.exp(ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listCells, listClone, currFp, x, nMut)) * fnPriorDist.density(x);
//				y = fn_gY.sample();
//			}
		}
		System.out.println(samples);
		System.out.println(probabilities);
		double bestSample = samples.get(0);
		double bestProb = probabilities.get(0);
		for (int i = 1; i < nSamples; i++){
			if (probabilities.get(i) > bestProb){
				bestSample = samples.get(i);
				bestProb = probabilities.get(i);
			}
		}
		return bestSample;
	}
	*/
	
	/**
	 * Get the upper bound of the conditional posterior distribution of Fn
	 * Calculates the value of envelope function for doing rejection sampling
	 * @param listCells
	 * @param listClone
	 * @param ampErrLikelihoodObj
	 * @param fnPriorDist
	 * @param currFp
	 * @param nMut
	 * @return
	 * Created On: Aug 21, 2017
	 */
	public static double getFnPosteriorMax(ArrayList<SingleCell> listCells, ArrayList<SingleCell> listCells_tot ,ArrayList<Clone> listClone, 
										   SCSAmplificationErrorLikelihood ampErrLikelihoodObj,
										   BetaDistribution fnPriorDist, double currFp, int nMut,int i){
		double maxFnPosterior = Double.NEGATIVE_INFINITY;
		
		double beta = 0.01;
		for (int f = 0; f<99; f++){
			double likelihood = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixIndividualMutation(listCells, listCells_tot, listClone, currFp, beta, nMut,i);
			double prior = Math.log(fnPriorDist.density(beta));
			double posterior = likelihood + prior;
//			System.out.println(posterior);
			if (posterior > maxFnPosterior)
        		maxFnPosterior = posterior;
			beta += 0.01;
		}
		return maxFnPosterior;
	}
	
	/**
	 * Sample a new Fn value from the conditional posterior using metropolis-hastings
	 * @param listCells
	 * @param listClone
	 * @param ampErrLikelihoodObj
	 * @param currFp
	 * @param currFn
	 * @param fnPriorDist
	 * @param nMut
	 * @param sampleIter
	 * @return
	 * Created On: Oct 15, 2017
	 */
	/*
	public static double sampleNewFnPosteriorMH(ArrayList<SingleCell> listCells, ArrayList<Clone> listClone,
									 				   SCSAmplificationErrorLikelihood ampErrLikelihoodObj,
									                   double currFp, double currFn, BetaDistribution fnPriorDist, 
									                   int nMut, int sampleIter){
		ArrayList<Double> sampledFnList = new ArrayList<>();
		ArrayList<Double> sampledLikelihoodList = new ArrayList<>();
		sampledFnList.add(currFn);
		double likelihood = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixFpFn(listCells, listClone, currFp, currFn, nMut);
		sampledLikelihoodList.add(likelihood);
		for (int st = 1; st < sampleIter; st++){
			int nSamples = sampledFnList.size();
			double nowFn = sampledFnList.get(nSamples-1);
			double nowLLH = sampledLikelihoodList.get(nSamples-1);
			double fnSD = nowFn*0.5;
			double newFn = sampleNewErrorRateHS(nowFn, fnSD);
			double newLogLikelihood = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixFpFn(listCells, listClone, currFp, newFn, nMut);
			double logLikelihoodRatio = newLogLikelihood - nowLLH;
			double logPriorRatio = getLogPriorRatio(newFn, nowFn, fnPriorDist);
			double proposalRatio = getProposalRatio(nowFn, newFn, fnSD);
			double acRatio = Math.min(1, Math.exp(logLikelihoodRatio + logPriorRatio + proposalRatio));
			double rr = _rng.nextDouble();
			if (rr <= acRatio){
				sampledFnList.add(newFn);
				sampledLikelihoodList.add(newLogLikelihood);
			}
		}
		
//		System.out.println(sampledFnList);
//		System.out.println(sampledLikelihoodList);
		double bestSample = sampledFnList.get(0);
		double bestProb = sampledLikelihoodList.get(0);
		for (int i = 1; i < sampledFnList.size(); i++){
			if (sampledLikelihoodList.get(i) > bestProb){
				bestSample = sampledFnList.get(i);
				bestProb = sampledLikelihoodList.get(i);
			}
		}
		return bestSample;
//		return 0;
	}
	*/
	/**
	 * Compute the logPrior Ratio of two error values
	 * @param newVal
	 * @param currVal
	 * @param priorDist
	 * @return
	 * Created On: Oct 15, 2017
	 */
	public static double getLogPriorRatio(double newVal, double currVal, BetaDistribution priorDist){
		return Math.log(priorDist.density(newVal)) - Math.log(priorDist.density(currVal));
	}
	
	/**
	 * Compute the logProposal ratio when proposing a new error value
	 * @param currFn
	 * @param newFn
	 * @param fnSD
	 * @return
	 * Created On: Oct 15, 2017
	 */
	public static double getProposalRatio(double currFn, double newFn, double fnSD){
		NormalDistribution forwardNormalDist = new NormalDistribution(currFn, fnSD);
		NormalDistribution backwardNormalDist = new NormalDistribution(newFn, fnSD);
		return Math.log(backwardNormalDist.density(currFn)) - Math.log(forwardNormalDist.density(newFn));
	}
	
	/**********************************************************************
	 * Functions to sample new value of Fn, considering the doublet model *
	 **********************************************************************/
	
	/**
	 * Sample a new value of Fn from the conditional posterior distribution
	 * Implements a Rejection Sampler, considers the presence of doublets
	 * @param listCells
	 * @param listClone
	 * @param c2CellCloneList
	 * @param cellDoubletFlagArr
	 * @param ampErrLikelihoodObj
	 * @param currFp
	 * @param fnPriorDist
	 * @param nMut
	 * @param SCF
	 * @param dataFlag
	 * @return
	 * Created On: Oct 9, 2017
	 */
	/*
	public static double sampleNewFnDoublet(ArrayList<SingleCell> listCells, ArrayList<Clone> listClone,
											ArrayList<Integer> c2CellCloneList, int[] cellDoubletFlagArr,
											 SCSAmplificationErrorLikelihood ampErrLikelihoodObj,
											 double currFp, BetaDistribution fnPriorDist, int nMut,
											 SCFUtilityFunctions SCF, int dataFlag){
		// Obtain the envelope (Y_max) of the posterior distribution of Fn
		double FnPosteriorMax = getFnPosteriorMaxDoublet(listCells, listClone, c2CellCloneList, cellDoubletFlagArr, ampErrLikelihoodObj, fnPriorDist, currFp, nMut, SCF, dataFlag);
		// Uniform distributions for doing rejection sampling
		UniformRealDistribution fn_fX = new UniformRealDistribution(0.01, 0.99);         // Uniform distribution x ~ [x_min, x_max]
		UniformRealDistribution fn_gY = new UniformRealDistribution(0, 1);  // Uniform distribution ~ [0, Y_max] 
		double x = fn_fX.sample();  // Sample x ~ U(x_min, x_max)
		// calculate f(x), here f is conditional posterior
		double x_posterior = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listCells, listClone, c2CellCloneList, cellDoubletFlagArr, currFp, x, nMut, SCF, dataFlag) + Math.log(fnPriorDist.density(x));
		double y = fn_gY.sample(); // Sample y ~ U(0, Y_max)
		double newFn;
		while(true){
			x = fn_fX.sample();
			y = Math.log(fn_gY.sample()) + FnPosteriorMax;
			x_posterior = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listCells, listClone, c2CellCloneList, cellDoubletFlagArr, currFp, x, nMut, SCF, dataFlag) + Math.log(fnPriorDist.density(x));
			// Return the sample x if y <= f(x)
			if (y <= x_posterior){
				newFn = x;
				break;
			}
			// Repeat the sampling
//			else{
//				x = fn_fX.sample();
//				x_posterior = Math.exp(ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listCells, listClone, c2CellCloneList, cellDoubletFlagArr, currFp, x, nMut, SCF, dataFlag)) * fnPriorDist.density(x);
//				y = fn_gY.sample();
//			}
		}
		return newFn;
	}
	*/
	
	/**
	 * Get the upper bound of the conditional posterior distribution of Fn
	 * Calculates the value of envelope function for doing rejection sampling
	 * Considers the presence of doublets
	 * @param listCells
	 * @param listClone
	 * @param c2CellCloneList
	 * @param cellDoubletFlagArr
	 * @param ampErrLikelihoodObj
	 * @param fnPriorDist
	 * @param currFp
	 * @param nMut
	 * @param SCF
	 * @param dataFlag
	 * @return
	 * Created On: Oct 9, 2017
	 */
	/*
	public static double getFnPosteriorMaxDoublet(ArrayList<SingleCell> listCells, ArrayList<Clone> listClone, 
												 ArrayList<Integer> c2CellCloneList, int[] cellDoubletFlagArr,
										   		 SCSAmplificationErrorLikelihood ampErrLikelihoodObj,
										   		 BetaDistribution fnPriorDist, double currFp, int nMut,
										   		 SCFUtilityFunctions SCF, int dataFlag){
		double maxFnPosterior = Double.NEGATIVE_INFINITY;
		double beta = 0.01;
		for (int f = 0; f<99; f++){
			double likelihood = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listCells, listClone, c2CellCloneList, cellDoubletFlagArr, currFp, beta, nMut, SCF, dataFlag);
			double prior = Math.log(fnPriorDist.density(beta));
			double posterior = likelihood + prior;
			if (posterior > maxFnPosterior)
        		maxFnPosterior = posterior;
			beta += 0.01;
		}
		return maxFnPosterior;
	}
	*/
	
	
	/**********************************************************************
	 * Functions to sample new value of Fp, considering the singlet model *
	 **********************************************************************/
	
	/**
	 * Sample a new Fp value from conditional posterior distribution
	 * @param listCells
	 * @param listClone
	 * @param ampErrLikelihoodObj
	 * @param currFp
	 * @param currFn
	 * @param fpPriorDist
	 * @param nMut
	 * @param sampleIter
	 * @return
	 * Created On: Oct 15, 2017
	 */
	/*
	public static double sampleNewFpPosteriorMH(ArrayList<SingleCell> listCells, ArrayList<Clone> listClone,
			   SCSAmplificationErrorLikelihood ampErrLikelihoodObj,
            double currFp, double currFn, BetaDistribution fpPriorDist, 
            int nMut, int sampleIter){
		ArrayList<Double> sampledFpList = new ArrayList<>();
		ArrayList<Double> sampledLikelihoodList = new ArrayList<>();
		sampledFpList.add(currFp);
		double likelihood = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixFpFn(listCells, listClone, currFp, currFn, nMut);
		sampledLikelihoodList.add(likelihood);
		for (int st = 1; st < sampleIter; st++){
			int nSamples = sampledFpList.size();
			double nowFp = sampledFpList.get(nSamples-1);
			double nowLLH = sampledLikelihoodList.get(nSamples-1);
			double fpSD = nowFp*0.5;
			double newFp = sampleNewErrorRateHS(nowFp, fpSD);
			double newLogLikelihood = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixFpFn(listCells, listClone, newFp, currFn, nMut);
			double logLikelihoodRatio = newLogLikelihood - nowLLH;
			double logPriorRatio = getLogPriorRatio(newFp, nowFp, fpPriorDist);
			double proposalRatio = getProposalRatio(nowFp, newFp, fpSD);
			double acRatio = Math.min(1, Math.exp(logLikelihoodRatio + logPriorRatio + proposalRatio));
			double rr = _rng.nextDouble();
			if (rr <= acRatio){
				sampledFpList.add(newFp);
				sampledLikelihoodList.add(newLogLikelihood);
			}
		}
		//System.out.println(sampledFnList);
		//System.out.println(sampledLikelihoodList);
		double bestSample = sampledFpList.get(0);
		double bestProb = sampledLikelihoodList.get(0);
		for (int i = 1; i < sampledFpList.size(); i++){
			if (sampledLikelihoodList.get(i) > bestProb){
				bestSample = sampledFpList.get(i);
				bestProb = sampledLikelihoodList.get(i);
			}
		}
		return bestSample;
		//return 0;
	}
	*/
	public static double sampleNewFpHS(ArrayList<SingleCell> listCells,ArrayList<SingleCell> listCells_tot, ArrayList<Clone> listClone, 
			   SCSAmplificationErrorLikelihood ampErrLikelihoodObj,
			   BetaDistribution fpPriorDist, ArrayList<Double> currFn, int nMut){
		double FpPosteriorMax = getFpPosteriorMax(listCells, listCells_tot,listClone, ampErrLikelihoodObj, fpPriorDist, currFn, nMut);
		// Uniform distributions for doing rejection sampling
		UniformRealDistribution fn_fX = new UniformRealDistribution(0.01, 0.99);         // Uniform distribution x ~ [x_min, x_max]
		UniformRealDistribution fn_gY = new UniformRealDistribution(0, 1);  // Uniform distribution ~ [0, Y_max] 
		double x = fn_fX.sample();  // Sample x ~ U(x_min, x_max)
		// calculate f(x), here f is conditional posterior
		double x_posterior;
		double y = fn_gY.sample(); // Sample y ~ U(0, Y_max)
		double newFp;
//		int nSamples = 100;
//		ArrayList<Double> samples = new ArrayList<>();
//		ArrayList<Double> probabilities = new ArrayList<>();
//		int accepted = 0;
		while(true){
			//while(true){
			x = fn_fX.sample();
			x_posterior = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listCells, listCells_tot,listClone, x, currFn, nMut) + Math.log(fpPriorDist.density(x));
			y = Math.log(fn_gY.sample()) + FpPosteriorMax;
			// Return the sample x if y <= f(x)
			if (y <= x_posterior){
				newFp = x;
				break;
			}
			//// Repeat the sampling
			//else{
			//x = fn_fX.sample();
			//x_posterior = Math.exp(ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listCells, listClone, x, currFn, nMut)) * fpPriorDist.density(x);
			//y = fn_gY.sample();
			//}
		}
//		double bestSample = samples.get(0);
//		double bestProb = probabilities.get(0);
//		for (int i = 1; i < nSamples; i++){
//			if (probabilities.get(i) > bestProb){
//				bestSample = samples.get(i);
//				bestProb = probabilities.get(i);
//			}
//		}
		return newFp;
	}
	
	/**
	 * Sample a new value of Fp from the conditional posterior distribution
	 * Implements a Rejection Sampler
	 * @param listCells
	 * @param listClone
	 * @param ampErrLikelihoodObj
	 * @param fpPriorDist
	 * @param currFn
	 * @param nMut
	 * @return
	 * Created On: Aug 21, 2017
	 */
	/*
	public static double sampleNewFp(ArrayList<SingleCell> listCells, ArrayList<Clone> listClone, 
										   SCSAmplificationErrorLikelihood ampErrLikelihoodObj,
										   BetaDistribution fpPriorDist, double currFn, int nMut){
		double FpPosteriorMax = getFpPosteriorMax(listCells, listClone, ampErrLikelihoodObj, fpPriorDist, currFn, nMut);
		// Uniform distributions for doing rejection sampling
		UniformRealDistribution fn_fX = new UniformRealDistribution(0.01, 0.99);         // Uniform distribution x ~ [x_min, x_max]
		UniformRealDistribution fn_gY = new UniformRealDistribution(0, FpPosteriorMax);  // Uniform distribution ~ [0, Y_max] 
		double x = fn_fX.sample();  // Sample x ~ U(x_min, x_max)
		// calculate f(x), here f is conditional posterior
		double x_posterior = Math.exp(ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listCells, listClone, x, currFn, nMut)) * fpPriorDist.density(x);
		double y = fn_gY.sample(); // Sample y ~ U(0, Y_max)
		
		int nSamples = 100;
		ArrayList<Double> samples = new ArrayList<>();
		ArrayList<Double> probabilities = new ArrayList<>();
		int accepted = 0;
		while(accepted < nSamples){
//		while(true){
			x = fn_fX.sample();
			x_posterior = Math.exp(ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listCells, listClone, x, currFn, nMut)) * fpPriorDist.density(x);
			y = fn_gY.sample();
			// Return the sample x if y <= f(x)
			if (y <= x_posterior){
				samples.add(x);
				probabilities.add(x_posterior);
				accepted += 1;
			}
//			// Repeat the sampling
//			else{
//				x = fn_fX.sample();
//				x_posterior = Math.exp(ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listCells, listClone, x, currFn, nMut)) * fpPriorDist.density(x);
//				y = fn_gY.sample();
//			}
		}
		System.out.println(samples);
		System.out.println(probabilities);
		double bestSample = samples.get(0);
		double bestProb = probabilities.get(0);
		for (int i = 1; i < nSamples; i++){
			if (probabilities.get(i) > bestProb){
				bestSample = samples.get(i);
				bestProb = probabilities.get(i);
			}
		}
		return bestSample;
	}
	*/
	/**
	 * Get the upper bound of the conditional posterior distribution of Fp
	 * Calculates the value of envelope function for doing rejection sampling
	 * @param listCells
	 * @param listClone
	 * @param ampErrLikelihoodObj
	 * @param fpPriorDist
	 * @param currFn
	 * @param nMut
	 * @return
	 * Created On: Aug 21, 2017
	 */
	public static double getFpPosteriorMax(ArrayList<SingleCell> listCells,ArrayList<SingleCell> listCells_tot, ArrayList<Clone> listClone, 
			   SCSAmplificationErrorLikelihood ampErrLikelihoodObj,
			   BetaDistribution fpPriorDist, ArrayList<Double> currFn, int nMut){
		double maxFpPosterior = Double.NEGATIVE_INFINITY;
		double alpha = 0.01;
		for (int f = 0; f<99; f++){
			double likelihood = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrix(listCells, listCells_tot, listClone, alpha, currFn, nMut);
			double prior = Math.log(fpPriorDist.density(alpha));
			double posterior = likelihood + prior;
			if (posterior > maxFpPosterior)
				maxFpPosterior = posterior;
			alpha += 0.01;
		}
		return maxFpPosterior;
	}
	
	/**********************************************************************
	 * Functions to sample new value of Fp, considering the doublet model *
	 **********************************************************************/
			                               
	/**
	 * Sample a new value of Fp from the conditional posterior distribution
	 * Implements a Rejection Sampler, considers the presence of doublets
	 * @param listCells
	 * @param listClone
	 * @param c2CellCloneList
	 * @param cellDoubletFlagArr
	 * @param ampErrLikelihoodObj
	 * @param fpPriorDist
	 * @param currFn
	 * @param nMut
	 * @param SCF
	 * @param dataFlag
	 * @return
	 * Created On: Oct 9, 2017
	 */
	/*
	public static double sampleNewFpDoublet(ArrayList<SingleCell> listCells, ArrayList<Clone> listClone, 
											ArrayList<Integer> c2CellCloneList, int[] cellDoubletFlagArr,
											   SCSAmplificationErrorLikelihood ampErrLikelihoodObj,
											   BetaDistribution fpPriorDist, double currFn, int nMut,
											   SCFUtilityFunctions SCF, int dataFlag){
		double FpPosteriorMax = getFpPosteriorMaxDoublet(listCells, listClone, c2CellCloneList, cellDoubletFlagArr, ampErrLikelihoodObj, fpPriorDist, currFn, nMut, SCF, dataFlag);
		// Uniform distributions for doing rejection sampling
		UniformRealDistribution fn_fX = new UniformRealDistribution(0.01, 0.99);         // Uniform distribution x ~ [x_min, x_max]
		UniformRealDistribution fn_gY = new UniformRealDistribution(0, 1);  // Uniform distribution ~ [0, Y_max] 
		double x = fn_fX.sample();  // Sample x ~ U(x_min, x_max)
		// calculate f(x), here f is conditional posterior
		double x_posterior = Math.exp(ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listCells, listClone, c2CellCloneList, cellDoubletFlagArr, x, currFn, nMut, SCF, dataFlag)) * fpPriorDist.density(x);
		double y = fn_gY.sample(); // Sample y ~ U(0, Y_max)
		double newFp;
		while(true){
			x = fn_fX.sample();
			x_posterior = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listCells, listClone, c2CellCloneList, cellDoubletFlagArr, x, currFn, nMut, SCF, dataFlag) + Math.log(fpPriorDist.density(x));
			y = Math.log(fn_gY.sample()) + FpPosteriorMax;
			// Return the sample x if y <= f(x)
			if (y <= x_posterior){
				newFp = x;
				break;
			}
			// Repeat the sampling
//			else{
//				x = fn_fX.sample();
//				x_posterior = Math.exp(ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listCells, listClone, c2CellCloneList, cellDoubletFlagArr, x, currFn, nMut, SCF, dataFlag)) * fpPriorDist.density(x);
//				y = fn_gY.sample();
//			}
		}
		return newFp;
	}
	*/
	/**
	 * Get the upper bound of the conditional posterior distribution of Fp
	 * Calculates the value of envelope function for doing rejection sampling
	 * Considers the presence of doublets
	 * @param listCells
	 * @param listClone
	 * @param c2CellCloneList
	 * @param cellDoubletFlagArr
	 * @param ampErrLikelihoodObj
	 * @param fpPriorDist
	 * @param currFn
	 * @param nMut
	 * @param SCF
	 * @param dataFlag
	 * @return
	 * Created On: Oct 9, 2017
	 */
	/*
	public static double getFpPosteriorMaxDoublet(ArrayList<SingleCell> listCells, ArrayList<Clone> listClone, 
												 ArrayList<Integer> c2CellCloneList, int[] cellDoubletFlagArr,
										   		 SCSAmplificationErrorLikelihood ampErrLikelihoodObj,
										   		 BetaDistribution fpPriorDist, double currFn, int nMut,
										   		 SCFUtilityFunctions SCF, int dataFlag){
		double maxFpPosterior = Double.NEGATIVE_INFINITY;
		double alpha = 0.01;
		for (int f = 0; f<99; f++){
			double likelihood = ampErrLikelihoodObj.computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(listCells, listClone, c2CellCloneList, cellDoubletFlagArr, alpha, currFn, nMut, SCF, dataFlag);
			double prior = Math.log(fpPriorDist.density(alpha));
			double posterior = likelihood + prior;
			if (posterior > maxFpPosterior)
				maxFpPosterior = posterior;
			alpha += 0.01;
		}
		return maxFpPosterior;
	}
	*/
	/************************************************************************
	 * Function to sample new value of delta (doublet rate), doublet model *
	 ************************************************************************/
	
	/**
	 * Sample new value of delta from the conditional posterior
	 * \delta | Y, a_del, b_del  ~  P(Y | \delta) * Beta(\delta | a_del, b_del) ~ Beta (a_del + \sum_j Y_j, b_del + nCell - \sum_j Y_j)
	 * @param cellDoubletFlagArr
	 * @param doubletPriorDist
	 * @param nCell
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static double sampleNewDeltaDoublet(int[] cellDoubletFlagArr, BetaDistribution doubletPriorDist, int nCell){
		int nDoublet = 0;
		for (int i : cellDoubletFlagArr)
			nDoublet += i;
		BetaDistribution deltaPosterior = new BetaDistribution(doubletPriorDist.getAlpha() + nDoublet, doubletPriorDist.getBeta() + nCell - nDoublet);
		double newDelta = deltaPosterior.sample();
		return newDelta;
	}
	

	/**
	 * @param args
	 * Created On: Aug 21, 2017
	 */
	public static void main(String[] args) {
		System.out.println(Math.exp(-906.5846291760181));
	}

}
