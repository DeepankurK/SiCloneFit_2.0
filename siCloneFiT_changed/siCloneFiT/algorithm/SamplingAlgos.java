/**
 * Aug 7, 2017
 */
package siCloneFiT.algorithm;

import java.util.ArrayList;
import java.util.Random;

import siCloneFiT.objects.Clone;

/**
 * @author hz22
 * Aug 7, 2017
 */
public class SamplingAlgos {

	/**
	 * @param args
	 * Created On: Aug 7, 2017
	 */
	public static Random _rng =  new Random();
	
	/**
	 * Sample from a discrete distribution using inverse CDF
	 * @param distribution
	 * @return
	 * Created On: Aug 7, 2017
	 */
	public static int sampleInverseCDFDiscrete(ArrayList<Double> distribution){
		int val = -2;
		ArrayList<Double> cdf = new ArrayList<>();
		cdf.add(distribution.get(0));
		for (int i = 1; i < distribution.size(); i++)
			cdf.add(cdf.get(i-1) + distribution.get(i));
//		System.out.println(cdf);
		double rr = _rng.nextDouble();
//		System.out.println("rr="+rr);
		if (rr <= cdf.get(0))
			val = 0;
		else{
			for (int i = 1; i < distribution.size(); i++){
				if ((rr > cdf.get(i-1)) & (rr <= cdf.get(i)))
					val = i;
			}
		}
		// simple hack now, may or may not work, need to figure out why distribution has NaN
		if (val == -2){
			val = 0;
//			System.out.println("dist " + distribution);
//			System.out.println("cdf " + cdf);
		}
		return val;
	}
	
	/**
	 * Sample a clone from the distribution
	 * @param distribution
	 * @param listClone
	 * @return
	 * Created On: Aug 15, 2017
	 */
	public static Clone sampleInverseCDFDiscrete(ArrayList<Double> distribution, ArrayList<Clone> listClone){
		int val = -1;
		ArrayList<Double> cdf = new ArrayList<>();
		cdf.add(distribution.get(0));
		for (int i = 1; i < distribution.size(); i++)
			cdf.add(cdf.get(i-1) + distribution.get(i));

		double rr = _rng.nextDouble();

		if (rr <= cdf.get(0))
			val = 0;
		else{
			for (int i = 1; i < distribution.size(); i++){
				if ((rr > cdf.get(i-1)) & (rr <= cdf.get(i)))
					val = i;
			}
		}
		return listClone.get(val);
	}
	
	/** 
	 * Normalize a distribution 
	 * @param dist
	 * @return
	 * Created On: Aug 12, 2017
	 */
	public static ArrayList<Double> normalizeDist(ArrayList<Double> dist){
		double sumDist = 0.0;
		for (double d : dist)
			sumDist += d;
		ArrayList<Double> normalizedDist = new ArrayList<>();
		for (double d : dist)
			normalizedDist.add(d/sumDist);
		return normalizedDist;
				
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ArrayList<Double> d = new ArrayList<>();
		d.add(0.12);
		d.add(0.2);
		d.add(0.5);
		ArrayList<Double> nd = SamplingAlgos.normalizeDist(d);
		System.out.println(nd);
		
		ArrayList<Double> distribution = new ArrayList<>();
		ArrayList<Clone> listClone = new ArrayList<>();
		for (int i = 0; i < 310; i++){
			distribution.add(SamplingAlgos._rng.nextDouble());
			Clone C = new Clone(i);
			listClone.add(C);
		}
		ArrayList<Double> normDist = SamplingAlgos.normalizeDist(distribution);
		Clone sample = SamplingAlgos.sampleInverseCDFDiscrete(normDist, listClone);
		System.out.println(sample.cloneID);

	}

}
