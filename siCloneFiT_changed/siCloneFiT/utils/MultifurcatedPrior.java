package siCloneFiT.utils;

import cc.mallet.types.ConstantMatrix;
import cc.mallet.types.DenseMatrix;
import cern.colt.Arrays;

public class MultifurcatedPrior {
	
	/**
	 * Computes the possible number of trees with n taxa and m internal nodes
	 * Uses Felsenstein's recursive relation, T_n,m = (n+m -2)*T_n-1,m-1 + m*T_n-1,m
	 * @param n
	 * @param m
	 * @return
	 * Created On: Feb 7, 2018
	 */
	public static double computeNumberTrees(int n, int m){
		
		double[][] D = new double[m][n-1];
		for (int i = 0; i < n-1; i++){
			D[0][i] = 1;
		}
		for (int i = 1; i < m; i++){
			for (int j = i; j < n-1; j++){
				D[i][j] = (i+1+j+2 - 2)*D[i-1][j-1] + (i+1)*D[i][j-1];
			}
		}
		return D[m-1][n-2];
		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double a = MultifurcatedPrior.computeNumberTrees(100, 5);
		double b = MultifurcatedPrior.computeNumberTrees(100, 9);
		System.out.println(a);
		System.out.println(b);
		System.out.println(Math.log(a/b));
	}

}
