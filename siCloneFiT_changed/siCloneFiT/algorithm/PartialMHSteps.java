/**
 * Aug 14, 2017
 */
package siCloneFiT.algorithm;

import siCloneFiT.likelihood.SCSAmplificationErrorLikelihood;
import siCloneFiT.objects.Clone;
import siCloneFiT.objects.DoubletGenotype;
import siCloneFiT.objects.SingleCell;

/**
 * @author hz22
 * Aug 14, 2017
 */
public class PartialMHSteps {

	/**
	 * Compute Acceptance Ratio for the partial MH steps of the sampler
	 * @param alpha_0
	 * @param nCell
	 * @param errorLikelihoodRatio
	 * @param treeLikelihoodRatio
	 * @param treePriorRatio
	 * @param hastingsRatio
	 * @param jacobian
	 * @param moveFlag, if a new clone is added then 1, else 0
	 * @return
	 * Created On: Aug 14, 2017
	 */
	public static double getAcceptanceRatio(double alpha_0, int nCell, double errorLikelihoodRatio, double treeLikelihoodRatio, 
									 double treePriorRatio, double hastingsRatio, double jacobian, int moveFlag){
		double CRPprior;
		if (moveFlag == 1)
			CRPprior = alpha_0/(nCell -1);
		else
			CRPprior = (nCell -1)/alpha_0;
//		System.out.println("CRPprior = "+CRPprior);
		double MHG_Ratio = CRPprior * errorLikelihoodRatio * treeLikelihoodRatio * treePriorRatio * hastingsRatio * jacobian;
//		System.out.println("MHG = "+ MHG_Ratio);
		double acceptanceRatio = Math.min(1, MHG_Ratio);
		return acceptanceRatio;
	}
	
	/**
	 * Compute the CRP prior times the error likelihood 
	 * n_j,c/n-1 * F(X_j | G_c)
	 * @param C
	 * @param cell
	 * @param nMut
	 * @param nCell
	 * @param fp
	 * @param fn
	 * @param errLObj
	 * @return
	 * Created On: Aug 15, 2017
	 */
	public static double computeCRPErrorLikelihood(Clone C, SingleCell cell,SingleCell cell_tot, int nMut, int nCell, double fp, ArrayList<Double> fn, SCSAmplificationErrorLikelihood errLObj){
		double CRPprior;
		if (C.memberCellList.contains(cell.cellID))
			CRPprior = (double) (C.memberCellList.size() - 1)/(nCell - 1);
		else
			CRPprior = (double) C.memberCellList.size()/(nCell - 1);
		double ampErrLogLikelihood = errLObj.computeAmpErrorLogLikelihoodCellGTVector(C.cloneGTVector, cell.observedGTVector, cell_tot.observedGTVectorm, fp, fn, nMut);
		return CRPprior * Math.exp(ampErrLogLikelihood);
	}
	
	/**
	 * Compute the CRP prior times the error likelihood considering doublet
	 * n_j,c/n-1 * F(X_j | G_c) if Y_j = 0
	 * n_j,c/n-1 * F(X_j | G_c1 + G_c2) if Y_j = 1
	 * @param C
	 * @param cell
	 * @param doubletMateClone
	 * @param cellDoubletFlag
	 * @param nMut
	 * @param nCell
	 * @param fp
	 * @param fn
	 * @param df
	 * @param errLObj
	 * @return
	 * Created On: Oct 9, 2017
	 */
	public static double computeCRPErrorLikelihoodDoublet(Clone C, SingleCell cell, Clone doubletMateClone, int cellDoubletFlag, int nMut, int nCell, double fp, double fn, int df, SCSAmplificationErrorLikelihood errLObj){
		double CRPprior;
		if (C.memberCellList.contains(cell.cellID))
			CRPprior = (double) (C.memberCellList.size() - 1)/(nCell - 1);
		else
			CRPprior = (double) C.memberCellList.size()/(nCell - 1);
		double ampErrLogLikelihood;
		if (cellDoubletFlag == 0)
			ampErrLogLikelihood = errLObj.computeAmpErrorLogLikelihoodCellGTVector(C.cloneGTVector, cell.observedGTVector, fp, fn, nMut);
		else {
			Integer[] cellDoubletGtVector = DoubletGenotype.getDoubletGtVector(C.cloneGTVector, doubletMateClone.cloneGTVector, nMut, df);
			ampErrLogLikelihood = errLObj.computeAmpErrorLogLikelihoodCellGTVector(cellDoubletGtVector, cell.observedGTVector, fp, fn, nMut);
		}
		return CRPprior * Math.exp(ampErrLogLikelihood);	
	}
	
	/**
	 * Compute the CRP prior times the error likelihood 
	 * n_j,c/(alpha_0 + n-1) * F(X_j | G_c)
	 * @param C
	 * @param cell
	 * @param nMut
	 * @param nCell
	 * @param fp
	 * @param fn
	 * @param errLObj
	 * @param alpha_0
	 * @return
	 * Created On: Oct 3, 2017
	 */
	public static double computeCRPErrorLikelihood(Clone C, SingleCell cell, int nMut, int nCell, double fp, double fn, SCSAmplificationErrorLikelihood errLObj, double alpha_0){
		double CRPprior;
		if (C.memberCellList.contains(cell.cellID))
			CRPprior = (double) (C.memberCellList.size() - 1)/(alpha_0 + nCell - 1);
		else
			CRPprior = (double) C.memberCellList.size()/(alpha_0 + nCell - 1);
		double ampErrLogLikelihood = errLObj.computeAmpErrorLogLikelihoodCellGTVector(C.cloneGTVector, cell.observedGTVector, fp, fn, nMut);
		return CRPprior * Math.exp(ampErrLogLikelihood);
	}
	
	/**
	 * Compute the CRP prior times the error likelihood 
	 * for the new clone in tree CRP
	 * (alpha_0/(nClone + 1) *(alpha_0 + nCell - 1)) *  F(X_j | G_c)
	 * @param cloneGTVector
	 * @param cell
	 * @param nMut
	 * @param nCell
	 * @param fp
	 * @param fn
	 * @param errLObj
	 * @param alpha_0
	 * @param currNClone
	 * @return
	 * Created On: Oct 3, 2017
	 */
	public static double computeNewCloneTreeCRPErrorLikelihood(Integer[] cloneGTVector, SingleCell cell, 
															   int nMut, int nCell, double fp, double fn, 
															   SCSAmplificationErrorLikelihood errLObj, 
															   double alpha_0, int currNClone){
		double CRPprior = alpha_0/(currNClone + 1)*(alpha_0 + nCell - 1);
		double ampErrLogLikelihood = errLObj.computeAmpErrorLogLikelihoodCellGTVector(cloneGTVector, cell.observedGTVector, fp, fn, nMut);
		return CRPprior * Math.exp(ampErrLogLikelihood);
	}
	/**
	 * @param args
	 * Created On: Aug 14, 2017
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
