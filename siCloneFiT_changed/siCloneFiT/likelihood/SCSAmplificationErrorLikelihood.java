/**
 * Aug 3, 2017
 */
package siCloneFiT.likelihood;

import java.util.ArrayList;
import java.util.HashMap;

import jeigen.DenseMatrix;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.DoubletGenotype;
import siCloneFiT.objects.SingleCell;
import siCloneFiT.utils.SCFUtilityFunctions;

/**
 * Class that contains the calculations regarding the SCS Amplification Error Model
 * @author hz22
 * Aug 3, 2017
 */
public class SCSAmplificationErrorLikelihood {

	/**
	 * 
	 */
	private static double _fp;
    private static ArrayList<Double> _fn;
    //private static double[][] errorMatrix;
    private static int dataFlag;
    
	public SCSAmplificationErrorLikelihood(double fp, ArrayList<Double> fn, int df) {
		_fp = fp;
		_fn = fn;
		dataFlag = df;
	}
	
	/**
	 * Sets a new false positive error and a new false negative error rate
	 * @param newFp, the new false positive error rate
	 * @param newFn, the new false negative error rate
	 * @param df, the dataFlag (0 if binary, 1 if ternary)
	 * @return true if error rates are valid and set, false otherwise
	 * Aug 3, 2017
	 */
	public boolean setFpFn(double newFp, ArrayList<Double> newFn, int df){
		if (newFp > 1.0 || newFp < 0.0) {
            return false;
        }
        for(Double beta : newFn){
        	if (beta > 1.0 || beta < 0.0) {
            	return false;
            }
        }

		_fp = newFp;
		_fn = newFn;
		return true;
	}
	
	/**
	 * Sets a new false positive error rate
     * @param newFp, the new false positive error rate
     * @param df, the dataFlag (0 if binary, 1 if ternary)
     * @return true if error rate is valid and set, false otherwise
	 * Aug 3, 2017
	 */
	public boolean setFp(double newFp, int df) {
        if (newFp > 1.0 || newFp < 0.0) {
            return false;
        }
        _fp = newFp;
        if (df == 1){
			errorMatrix[0][0] = (1 - _fp*_fn/2.0 - _fp);
			errorMatrix[0][1] = _fp;
			errorMatrix[0][2] = _fp*_fn/2.0;
        }
        else{
			errorMatrix[0][0] = 1 - _fp;
			errorMatrix[0][1] = _fp;
        }
        return true;
	}
	
	/**
     * Sets a new false negative error rate
     * @param newFn, the new false negative error rate
     * @param df, the dataFlag (0 if binary, 1 if ternary)
     * @return true if error rate is valid and set, false otherwise
     */
    public boolean setFn(ArrayList<Double> newFn, int df) {
        for(Double beta: newFn){
        	if (beta > 1.0 || beta < 0.0) {
            	return false;
        	}
        }
        _fn = newFn;
		return true;
    }
    
    /**
     * Return the error matrix
     * @return
     * Created On: Oct 3, 2017
     */
    public DenseMatrix getErrorMat(){
    	return new DenseMatrix(errorMatrix);
    }
    
    /**
     * Compute the Log-likelihood of a cell belonging to a clone, P(X_j | G_k, alpha, beta) 
     * @param cellCloneCandidateGTVector, the GT vector for the clone k, G_k
     * @param cellObsGTVector, X_j, the observed GT vector for cell j
     * @param fp, alpha
     * @param fn, beta
     * @param nMut, n, the number of mutations
     * @return
     * Aug 3, 2017
     */
    public double computeAmpErrorLogLikelihoodCellGTVector(Integer[] cellCloneCandidateGTVector, 
    													   Integer[] cellObsGTVector,Integer[] cellObsTotVector, double fp,ArrayList<Double> fn,
    													   int nMut){
    	double logProbVector = 0.0;
    	setFpFn(fp, fn, dataFlag);
    	for (int i = 0; i < nMut; i++){
//    		System.out.println(computeErrorLikelihoodGT(cellCloneCandidateGTVector[i], cellObsGTVector[i]));
    		logProbVector += computeErrorLikelihoodGT(cellCloneCandidateGTVector[i], cellObsGTVector[i],cellObsTotVector[i],_fn[i],_fp);    		
    	}
    	return logProbVector;
    }
    public double computeAmpErrorLogLikelihoodCellGTVectorInvidualMutation(Integer[] cellCloneCandidateGTVector, 
                                                           Integer[] cellObsGTVector,Integer[] cellObsTotVector, double fp,double fn,
                                                           int nMut,int j){
        double logProbVector = 0.0;
        if (fn > 1.0 || fn < 0.0 || fp > 1.0 || fp < 0.0) {
            logProbVector = computeErrorLikelihoodGT(cellCloneCandidateGTVector[j], cellObsGTVector[j],cellObsTotVector[j],_fn[j],_fp);
        else
            logProbVector = computeErrorLikelihoodGT(cellCloneCandidateGTVector[j], cellObsGTVector[j],cellObsTotVector[j],fn,fp);
        return logProbVector;
    }
    
    /**
     * Compute log-likelihood of the whole matrix
     * @param listClone
     * @param cloneGenotypeMat
     * @param cellIDGTVectorMap
     * @param nMut
     * @return
     * Created On: Oct 14, 2017
     */
    public double computeAmpErrorLoglikelihoodCellMatrix(ArrayList<Clone> listClone, ArrayList<Integer[]> cloneGenotypeMat, HashMap<Integer, Integer[]> cellIDGTVectorMap, int nMut){
    	double loglikelihood = 0.0;
    	for (int i = 0; i < listClone.size(); i++){
    		Integer[] cloneGtVector = cloneGenotypeMat.get(i);
    		for (Integer cellID: listClone.get(i).memberCellList){
    			loglikelihood += computeAmpErrorLoglikelihoodCellvector(cloneGtVector, cellIDGTVectorMap.get(cellID), nMut);
    		}
    	}
    	return loglikelihood; 
    }
    
    /**
     * Compute log likelihood of a clone
     * @param cloneGtVector
     * @param cellGtVector
     * @param nMut
     * @return
     * Created On: Oct 14, 2017
     */
    public double computeAmpErrorLoglikelihoodCellvector(Integer[] cloneGtVector, Integer[] cellGtVector, int nMut){
    	double logProbVector = 0.0;
    	for (int i = 0; i < nMut; i++){
    		logProbVector += computeErrorLikelihoodGT(cloneGtVector[i], cellGtVector[i]);   
    	}
    	return logProbVector;
    }
    
    /**
     * Compute the total Log-likelihood of all cells 
     * @param listCells
     * @param listClone
     * @param fp
     * @param fn
     * @param nMut
     * @return
     * Created On: Aug 12, 2017
     */
    public double computeAmpErrorLogLikelihoodCellMatrix(ArrayList<SingleCell> listCells, ArrayList<SingleCell> listCells_tot, ArrayList<Clone> listClone,
    													double fp, ArrayList<Double> fn, int nMut){
    	double log_likelihood = 0.0;
    	setFn(fn, dataFlag);
    	for (Clone C : listClone){
    		for (int cellID: C.memberCellList){
    			log_likelihood += computeAmpErrorLogLikelihoodCellGTVector(C.cloneGTVector, listCells.get(cellID).observedGTVector,listCells_tot.get(cellID).observedGTVector, fp, fn, nMut);
    		}
    	}
//    	for (SingleCell s: listCells){
//    		Clone cellClone = listClone.get(s.cloneID);
//    		log_likelihood += computeAmpErrorLogLikelihoodCellGTVector(cellClone.cloneGTVector, s.observedGTVector, fp, fn, nMut);
//    	}
    	return log_likelihood;
    }
    /**
     * Compute the total Log-likelihood of a specific mutation in the data.
     * @param listCells
     * @param listClone
     * @param fp
     * @param fn
     * @param nMut
     * @return
     * Created On: Aug 12, 2017
     */
    public double computeAmpErrorLogLikelihoodCellMatrixIndividualMutation(ArrayList<SingleCell> listCells, ArrayList<SingleCell> listCells_tot, ArrayList<Clone> listClone,
                                                        double fp, double fn, int nMut,int j){
        double log_likelihood = 0.0;
        for (Clone C : listClone){
            for (int cellID: C.memberCellList){
                log_likelihood += computeAmpErrorLogLikelihoodCellGTVectorInvidualMutation(C.cloneGTVector, listCells.get(cellID).observedGTVector,listCells_tot.get(cellID).observedGTVector, fp, fn, nMut,j);
            }
        }
        return log_likelihood;
    }
    
    public double computeAmpErrorLogLikelihoodCellMatrixFpFn(ArrayList<SingleCell> listCells,ArrayList<SingleCell> listCells_tot, ArrayList<Clone> listClone,
			double fp, ArrayList<Double> fn, int nMut){
    	double log_likelihood = 0.0;
    	setFpFn(fp, fn, dataFlag);
    	for (Clone C : listClone){
    		for (int cellID: C.memberCellList){
    			log_likelihood += computeAmpErrorLogLikelihoodCellGTVector(C.cloneGTVector, listCells.get(cellID).observedGTVector,listCells_tot.get(cellID).observedGTVector,fp, fn, nMut);
    		}
    	}
    	//for (SingleCell s: listCells){
    	//Clone cellClone = listClone.get(s.cloneID);
    	//log_likelihood += computeAmpErrorLogLikelihoodCellGTVector(cellClone.cloneGTVector, s.observedGTVector, fp, fn, nMut);
    	//}
    	return log_likelihood;
    }
    
    
    /**
     * Compute the error log-likelihood of an observed genotype matrix, 
     * given true clonal genotype matrix in the presence of doublets
     * \sum_{j=1}^m P(X_j| c1_j, c2_j, Y_j, G, \alpha, \beta)
     * @param listCells, list of cells
     * @param listClone, list of clones
     * @param c2CellCloneList, list of cloneIDs for 2nd clone of each cell (c2)
     * @param cellDoubletFlagArr, Y, array of flags for doublets
     * @param fp
     * @param fn
     * @param nMut
     * @param SCF
     * @param dataFlag
     * @return
     * Created On: Oct 9, 2017
     */
    public double computeAmpErrorLogLikelihoodCellMatrixDoubletFpFn(ArrayList<SingleCell> listCells, ArrayList<Clone> listClone,
    																ArrayList<Integer> c2CellCloneList, int[] cellDoubletFlagArr,
    																double fp, double fn, int nMut, SCFUtilityFunctions SCF, int dataFlag){
    	double log_likelihood = 0.0;
    	setFpFn(fp, fn, dataFlag);
    	for (Clone C : listClone){
    		for (int cellID: C.memberCellList){
    			// This cell is a doublet
    			// Compute the merged genotype 
    			if (cellDoubletFlagArr[cellID] == 1){
    				Clone doubletClone2 = SCF.getCloneFrmList(listClone, c2CellCloneList.get(cellID));
    				Integer[] cellDoubletGtVector = DoubletGenotype.getDoubletGtVector(C.cloneGTVector, doubletClone2.cloneGTVector, nMut, dataFlag);
    				log_likelihood += computeAmpErrorLogLikelihoodCellGTVector(cellDoubletGtVector, listCells.get(cellID).observedGTVector, fp, fn, nMut);
    			}
    			// cell is a singlet
    			else{
    				log_likelihood += computeAmpErrorLogLikelihoodCellGTVector(C.cloneGTVector, listCells.get(cellID).observedGTVector, fp, fn, nMut);
    			}
    		}
    	}
    	return log_likelihood;
    }
    
    
	
	/**
	 * compute the likelihood of a single observed genotype given the true genotype
	 * @param trueGT
	 * @param obsGT
	 * @param fp
	 * @param fn
	 * @return
	 * Aug 3, 2017
	 */

	public double binomial(int N, int k, double p) {
        double[][] b = new double[N+1][k+1];

        // base cases
        for (int i = 0; i <= N; i++)
            b[i][0] = Math.pow(1.0 - p, i);
        b[0][0] = 1.0;

        // recursive formula
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= k; j++) {
                b[i][j] = p * b[i-1][j-1] + (1.0 - p) *b[i-1][j];
            }
        }
        return b[N][k];
    }

	public double computeErrorLikelihoodGT(Integer trueGT, Integer obsGT,Integer obsTot,double fn,double fp){
		if trueGT==1{
			return Math.log(binomial(obsTot+obsGT,obsGT,fn));
		}
		else{
			return Math.log(binomial(obsTot+obsGT,obsGT,fp));
		}
	}
	/**
	 * Compute the Error Model Log-likelihood for a genotype gt of a clone at a locus i
	 * \sum_{j \in Clone_k} log E(X_ji | G_ki = gt) 
	 * @param gt, the candidate genotype, {0,1,2}
	 * @param mutIndex, index of the locus, 0,...,n-1
	 * @param clone, Clone 
	 * @param cellGTVectorMap, The HashMap of observed genotype matrix
	 * @return the log-likelihood value
	 * Created On: Aug 8, 2017
	 */
	public double computeCloneErrorLogLikelihood(Integer gt, int mutIndex, Clone clone, 
											  HashMap<Integer, Integer[]> cellGTVectorMap){
		double logLikelihood = 0.0;
		for (Integer cellID: clone.memberCellList){
			Integer obsCellGT = cellGTVectorMap.get(cellID)[mutIndex];
			double thisMutLikelihood = this.computeErrorLikelihoodGT(gt, obsCellGT);
			if (Double.isNaN(thisMutLikelihood)){
				System.out.println("cell gt = " + obsCellGT);
				System.out.println(gt);
				System.out.println(errorMatrix[gt][obsCellGT]);
			}
			logLikelihood += thisMutLikelihood;
		}
		return logLikelihood;		
	}
	
	/**
	 * Compute the log-likelihood ratio of two clonal genotypes (also vectors if only one position is changed)
	 * @param currGt
	 * @param newGt
	 * @param mutIndex
	 * @param clone
	 * @param cellGTVectorMap
	 * @return
	 * Created On: Oct 13, 2017
	 */
	public double computeCloneNewGenotypeErrorLogLikelihoodRatio(Integer currGt, Integer newGt, int mutIndex, Clone clone, HashMap<Integer, Integer[]> cellGTVectorMap){
		double currGtLogLikelihood = computeCloneErrorLogLikelihood(currGt, mutIndex, clone, cellGTVectorMap);
		double newGtLogLikelihood = computeCloneErrorLogLikelihood(newGt, mutIndex, clone, cellGTVectorMap);
		return newGtLogLikelihood - currGtLogLikelihood;
	}
	
	/**
	 * Compute the Error Model Log-likelihood for a genotype gt of a clone at a locus i
	 * considers presence of doublets
	 * \sum_{j \in Clone_k} log E(X_ji | G_ki = gt, c2_j_i, \alpha, \beta) 
	 * @param gt
	 * @param mutIndex
	 * @param clone
	 * @param listClone
	 * @param c2CellCloneList
	 * @param cellDoubletFlagArr
	 * @param df
	 * @param cellGTVectorMap
	 * @param SCF
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public double computeCloneErrorLogLikelihoodDoublet(Integer gt, int mutIndex, Clone clone, ArrayList<Clone> listClone, 
			   											ArrayList<Integer> c2CellCloneList, int[] cellDoubletFlagArr, int df,
			  											HashMap<Integer, Integer[]> cellGTVectorMap, SCFUtilityFunctions SCF){
		double logLikelihood = 0.0;
		for (Integer cellID: clone.memberCellList){
			Integer obsCellGT = cellGTVectorMap.get(cellID)[mutIndex];
			double thisMutLikelihood;
			// cell is singlet
			if (cellDoubletFlagArr[cellID] == 0)
				thisMutLikelihood = this.computeErrorLikelihoodGT(gt, obsCellGT);
			// cell is doublet
			else{
				Clone cellDoubletMateClone = SCF.getCloneFrmList(listClone, c2CellCloneList.get(cellID));
				int doubletGt = DoubletGenotype.getDoubletGt(gt, cellDoubletMateClone.cloneGTVector[mutIndex], df);
				thisMutLikelihood = this.computeErrorLikelihoodGT(doubletGt, obsCellGT);
			}
			logLikelihood += thisMutLikelihood;	
				
		}
		return logLikelihood;	
	}
	
	/**
	 * Compute the error likelihood ratio for a new proposed clone for the cell
	 * @param nuCloneGTVector
	 * @param oldCloneGTVector
	 * @param cellObsGTVector
	 * @param fp
	 * @param fn
	 * @param nMut
	 * @return
	 * Created On: Aug 14, 2017
	 */
	public double computeNewCloneErrorLikelihoodRatio(Integer[] nuCloneGTVector, Integer[] oldCloneGTVector,
    												  Integer[] cellObsGTVector, Integer[] cellObsTotVector,double fp, ArrayList<Double> fn,
    												  int nMut){
		double nuCloneLoglikelihood = this.computeAmpErrorLogLikelihoodCellGTVector(nuCloneGTVector, cellObsGTVector, cellObsTotVector,fp, fn, nMut);
//		System.out.println(" ");
		double oldCloneLoglikelihood = this.computeAmpErrorLogLikelihoodCellGTVector(oldCloneGTVector, cellObsGTVector, cellObsTotVector,fp, fn, nMut);
//		System.out.printf("nuClone = %f\toldClone = %f\n", nuCloneLoglikelihood, oldCloneLoglikelihood);
		return Math.exp(nuCloneLoglikelihood - oldCloneLoglikelihood);
	}
	
	/**
	 * Compute the Likelihood ratio for proposing new value of Fn
	 * @param listCells
	 * @param listClone
	 * @param fp
	 * @param currfn
	 * @param newFn
	 * @param nMut
	 * @return
	 * Created On: Aug 21, 2017
	 */
	public double computeNewFnErrorLikelihoodRatio(ArrayList<SingleCell> listCells, ArrayList<Clone> listClone,
			double fp, double currfn, double newFn, int nMut){
		double currFnLoglikelihood = this.computeAmpErrorLogLikelihoodCellMatrix(listCells, listClone, fp, currfn, nMut);
		setFn(newFn, dataFlag);
		double newFnLoglikelihood = this.computeAmpErrorLogLikelihoodCellMatrix(listCells, listClone, fp, newFn, nMut);
		setFn(currfn, dataFlag);
		return Math.exp(newFnLoglikelihood - currFnLoglikelihood);
	}

	/**
	 * @param args
	 * Aug 3, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
