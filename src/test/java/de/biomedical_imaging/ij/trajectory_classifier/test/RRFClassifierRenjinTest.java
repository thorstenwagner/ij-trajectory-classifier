package de.biomedical_imaging.ij.trajectory_classifier.test;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;

import de.biomedical_imaging.ij.trajectory_classifier.RRFClassifierRenjin;
import de.biomedical_imaging.ij.trajectory_classifier.TraJClassifier_;
import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.simulation.AbstractSimulator;
import de.biomedical_imaging.traJ.simulation.CentralRandomNumberGenerator;
import de.biomedical_imaging.traJ.simulation.FreeDiffusionSimulator;

public class RRFClassifierRenjinTest {
	

	@Test
	public void Classify_FreeDiffusionSingleTrack_Test() {
		/*
		 * Export model!
		 */
		TraJClassifier_ tclass = new TraJClassifier_();
		String modelpath="";
		try {
			modelpath=tclass.ExportResource("/randomForestModel.RData");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		CentralRandomNumberGenerator.getInstance().setSeed(8);
		double diffusioncoefficient = 9.02*Math.pow(10,-2); //[µm^2/s];
		double timelag = 1.0/30;
		int simtracklength = 500;

		AbstractSimulator sim = new FreeDiffusionSimulator(diffusioncoefficient, timelag, 2, simtracklength);
		
		RRFClassifierRenjin c = new RRFClassifierRenjin(modelpath);
		c.start();
		String res = c.classify(sim.generateTrajectory());
		System.out.println("RES: " + res);
		c.stop();
		assertEquals("NORM. DIFFUSION", res);
	}
	
	@Test
	public void Classify_FreeDiffusionMultipleTracks_Test() {
		/*
		 * Export model!
		 */
		TraJClassifier_ tclass = new TraJClassifier_();
		String modelpath="";
		try {
			modelpath=tclass.ExportResource("/randomForestModel.RData");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		CentralRandomNumberGenerator.getInstance().setSeed(8);
		double diffusioncoefficient = 9.02*Math.pow(10,-2); //[µm^2/s];
		double timelag = 1.0/30;
		int simtracklength = 500;

		AbstractSimulator sim = new FreeDiffusionSimulator(diffusioncoefficient, timelag, 2, simtracklength);
		
		RRFClassifierRenjin c = new RRFClassifierRenjin(modelpath);
		c.start();
		ArrayList<Trajectory> tracks = new ArrayList<Trajectory>();
		tracks.add(sim.generateTrajectory());
		tracks.add(sim.generateTrajectory());
		String[] res = c.classify(tracks);
		System.out.println("RES: " + res[0]);
		c.stop();
		assertEquals("NORM. DIFFUSION", res[0]);
		assertEquals("NORM. DIFFUSION", res[1]);
	}

}
