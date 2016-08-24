package de.biomedical_imaging.ij.trajectory_classifier.test;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import de.biomedical_imaging.ij.trajectory_classifier.RRFClassifierRenjin;
import de.biomedical_imaging.ij.trajectory_classifier.TraJClassifier_;
import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.TrajectoryUtil;
import de.biomedical_imaging.traJ.simulation.AbstractSimulator;
import de.biomedical_imaging.traJ.simulation.ActiveTransportSimulator;
import de.biomedical_imaging.traJ.simulation.AnomalousDiffusionWMSimulation;
import de.biomedical_imaging.traJ.simulation.CentralRandomNumberGenerator;
import de.biomedical_imaging.traJ.simulation.ConfinedDiffusionSimulator;
import de.biomedical_imaging.traJ.simulation.FreeDiffusionSimulator;

public class RRFClassifierRenjinTest {
	static RRFClassifierRenjin c;
	private static double timelag= 1.0/30;
	
	@BeforeClass
	public static void setup(){
		
		TraJClassifier_ tclass = new TraJClassifier_();
		String modelpath="";
		try {
			modelpath=tclass.ExportResource("/randomForestModel.RData");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		c = new RRFClassifierRenjin(modelpath, timelag);
		c.start();
	
	}
	
	@AfterClass
	public static void shutdown(){
		c.stop();
	}

	@Test
	public void Classify_FreeDiffusionSingleTrack_Test() {
		CentralRandomNumberGenerator.getInstance().setSeed(8);
		double diffusioncoefficient = 9.02; //[µm^2/s];
		int simtracklength = 500;

		AbstractSimulator sim = new FreeDiffusionSimulator(diffusioncoefficient, timelag, 2, simtracklength);
		String res = c.classify(sim.generateTrajectory());
		

		assertEquals("NORM. DIFFUSION", res);
	}
	
	@Test
	public void Classify_FreeDiffusionMultipleTracks_Test() {

		CentralRandomNumberGenerator.getInstance().setSeed(9);
		double diffusioncoefficient = 9.02; //[µm^2/s];
		int simtracklength = 500;

		AbstractSimulator sim = new FreeDiffusionSimulator(diffusioncoefficient, timelag, 2, simtracklength);
		

		ArrayList<Trajectory> tracks = new ArrayList<Trajectory>();
		tracks.add(sim.generateTrajectory());
		tracks.add(sim.generateTrajectory());
		String[] res = c.classify(tracks);
		assertEquals("NORM. DIFFUSION", res[0]);
		assertEquals("NORM. DIFFUSION", res[1]);
	}
	
	@Test
	public void Classify_AnomalousDiffusionSingleTrack_Test() {
		
		
		CentralRandomNumberGenerator.getInstance().setSeed(8);
		double diffusioncoefficient = 9.02; //[µm^2/s];
		int simtracklength = 200;

		AbstractSimulator sim = new AnomalousDiffusionWMSimulation(diffusioncoefficient, timelag, 2, simtracklength, 0.5);
		Trajectory t = sim.generateTrajectory();
		String res = c.classify(t);

		assertEquals("SUBDIFFUSION", res);
	}
	
	@Test
	public void Classify_ConfinedDiffusionSingleTrack_Test() {
		
		
		CentralRandomNumberGenerator.getInstance().setSeed(8);
		double diffusioncoefficient = 9.02; //[µm^2/s];
		int simtracklength = 200;
		double radius_confined = Math.sqrt(-1*Math.log(0.9)*(4*diffusioncoefficient*simtracklength*timelag));
		AbstractSimulator sim = new ConfinedDiffusionSimulator(diffusioncoefficient, timelag, radius_confined, 2, simtracklength);
		
		String res = c.classify(sim.generateTrajectory());

		assertEquals("CONFINED", res);
	}
	
	@Test
	public void Classify_PureDirectedMotionSingleTrack_Test() {
		
		
		CentralRandomNumberGenerator.getInstance().setSeed(8);
		int simtracklength = 200;
		double velocity = 1.5; // µm/s
		AbstractSimulator sim = new ActiveTransportSimulator(velocity, Math.PI/4.0, timelag, 2, simtracklength);
		
		String res = c.classify(sim.generateTrajectory());

		assertEquals("DIRECTED/ACTIVE", res);
	}
	
	@Test
	public void Classify_DiffusionWithDriftSingleTrack_Test() {
		
		
		CentralRandomNumberGenerator.getInstance().setSeed(8);
		double diffusioncoefficient = 9.02*Math.pow(10,-2); //[µm^2/s];
		int simtracklength = 200;
		double aToDRatio = 10;
		double tracklength = simtracklength*timelag;
		double velocity = Math.sqrt(aToDRatio*4*diffusioncoefficient/tracklength);
		AbstractSimulator sim = new ActiveTransportSimulator(velocity, Math.PI/4.0, timelag, 2, simtracklength);
		Trajectory active = sim.generateTrajectory();
		sim = new FreeDiffusionSimulator(diffusioncoefficient, timelag, 2, simtracklength);
		Trajectory free = sim.generateTrajectory();
		
		Trajectory combined = TrajectoryUtil.combineTrajectory(free, active);
		String res = c.classify(combined);

		assertEquals("DIRECTED/ACTIVE", res);
	}
	

}
