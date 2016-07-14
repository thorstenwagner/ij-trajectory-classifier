package de.biomedical_imaging.ij.trajectory_classifier.test;

import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.Test;

import de.biomedical_imaging.ij.trajectory_classifier.WeightedWindowedClassificationProcess;

public class WeightedWindowedClassificationProcessTest extends WeightedWindowedClassificationProcess {

	@Test
	public void applyWeighteningTest_allequal() {
		int tracklength = 1000;
		int windowsize = 30;
		String[] res = new String[tracklength-2*windowsize];
		Arrays.fill(res,"A");
		double[] confidence = new double[tracklength-2*windowsize];
		Arrays.fill(confidence,0.5);
		String[] result = this.applyWeightening(res, confidence, windowsize,tracklength);
		
		String[] exp = new String[tracklength];
		Arrays.fill(exp,"A");
		assertEquals(exp, result);
	}
	
	/*
	 * Because low confidence entry have a low weight, the result should
	 * only contain the high confidence values. Only the first and the last value
	 * can differ, because they only get one vote.
	 */
	@Test
	public void applyWeighteningTest_low_high_confidence() {
		int tracklength = 1000;
		int windowsize = 30;
		String[] res = new String[tracklength-2*windowsize];
		
		double[] confidence = new double[tracklength-2*windowsize];
		for(int i = 0; i < res.length; i++){
			if(i%2 == 0){
				res[i] = "A";
				confidence[i] = 0.9;
			}else{
				res[i] = "B";
				confidence[i] = 0.1;
			}
		}
		String[] result = this.applyWeightening(res, confidence, windowsize,tracklength);
		String[] exp = new String[tracklength];
		Arrays.fill(exp,"A");
		
		result = Arrays.copyOfRange(result, 1, 998);
		exp = Arrays.copyOfRange(exp, 1, 998);
		assertEquals(exp, result);
	}

}
