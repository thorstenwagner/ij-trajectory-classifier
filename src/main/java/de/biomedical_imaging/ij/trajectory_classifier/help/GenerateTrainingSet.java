package de.biomedical_imaging.ij.trajectory_classifier.help;

import java.util.ArrayList;

import de.biomedical_imaging.traJ.ExportTools;
import de.biomedical_imaging.traJ.Trajectory;

import de.biomedical_imaging.traJ.features.BoundednessFeature;
import de.biomedical_imaging.traJ.simulation.AbstractSimulator;
import de.biomedical_imaging.traJ.simulation.ActiveTransportSimulator;
import de.biomedical_imaging.traJ.simulation.AnomalousDiffusionWMSimulation;
import de.biomedical_imaging.traJ.simulation.CentralRandomNumberGenerator;
import de.biomedical_imaging.traJ.simulation.CombinedSimulator;
import de.biomedical_imaging.traJ.simulation.ConfinedDiffusionSimulator;
import de.biomedical_imaging.traJ.simulation.FreeDiffusionSimulator;

import de.biomedical_imaging.traJ.simulation.SimulationUtil;

public class GenerateTrainingSet {
	/*
	 *  This script generates training trajecotires
	 *  
	 *  There are 4 different diffusion modes:
	 *   - Free diffusion
	 *   - Subdiffusion
	 *   - Confined diffusion (Radius is set to 5µm)
	 *   - Active transport + diffusion
	 *  
	 *  For each trajectory:
	 *  The signal to noise ratio is randomly choosen between 1 and 20
	 *  For confined diffusion trajectories, the boundedness is choosen randomly between 1 and 16
	 *  For subdiffusion trajectories, alpha is choosen randomly between 0.9 and 0.1
	 *  For active transport + diffusion, the ratio between active transport and diffusion is randomly choosen between 1 and 30.
	 *  
	 */
	enum SIM_TYPE {
		FREE,
	    ANOMALOUS,
	    CONFINED,
	    ACTIVE
	}
	
	private static CentralRandomNumberGenerator r;
	
	public static void main(String[] args) {
		
		final int MODE_TRAINING = 1;
		final int MODE_TEST = 2;
		final int MODE_VALIDATION = 3;
		
		
		//General Parameters
		int numberOfTracks = 0;					
		int seed = 0;
		int MODE = MODE_TRAINING; // SELECT WHICH TYPE OF DATA YOU WANT TO GENERATE
		String prefix = "";
		switch (MODE) {
		case MODE_TRAINING:
			numberOfTracks = 5000;
			seed = 22;
			prefix = "training";
			break;
		case MODE_TEST:
			numberOfTracks = 1000;
			prefix = "test";
			seed = 23;
			break;
		case MODE_VALIDATION:
			numberOfTracks = 250;
			prefix = "validation";
			seed = 24;
			break;

		default:
			break;
		}
		String path = "/home/thorsten/tracks_"+prefix+".RData";
		r  = CentralRandomNumberGenerator.getInstance();
		r.setSeed(seed); 

		double diffusioncoefficient = 9.02; //[µm^2/s];
		double timelag = 1.0/30; //s
		int dimension = 2;
		
		//Active transport / drift
		double angleVelocity = Math.PI/4.0; //rad/s
		

		SIM_TYPE[] types = SIM_TYPE.values();
		ArrayList<Trajectory> trajectorys= new ArrayList<Trajectory>();

		System.out.println("Generation of free, confined , subdiffsuion and active trajectories");
		int tCounter = 0;
				for (SIM_TYPE type : types) {
					
					for(int i = 0 ; i < numberOfTracks; i++){
							double tracklength = (1+r.nextDouble()*20);
							int numberOfSteps = (int)(tracklength * 1/timelag);
							double boundedness = 1 + r.nextDouble()*20;
							double alpha = 0.3+r.nextDouble()*0.4; // 0.1 - 0.9
							AbstractSimulator sim = null;
							String typestring = "";
							typestring += type.toString();
							Trajectory t = null;
							double diffusionToNoiseRatio = 1 + r.nextDouble()*8;
							double sigmaPosNoise = 1;
							switch (type) {
							case FREE:
								
								sim = new FreeDiffusionSimulator(diffusioncoefficient, timelag, dimension, numberOfSteps);
								sigmaPosNoise = Math.sqrt(diffusioncoefficient*timelag)/diffusionToNoiseRatio; 
								break;
							case CONFINED:
								
								double radius_confined = Math.sqrt(BoundednessFeature.a(numberOfSteps)*diffusioncoefficient*timelag/(4*boundedness));
								sim = new ConfinedDiffusionSimulator(diffusioncoefficient, timelag, radius_confined, dimension, numberOfSteps);
								sigmaPosNoise = Math.sqrt(diffusioncoefficient*timelag)/diffusionToNoiseRatio; 
								break;
							case ACTIVE:
								double aToDRatio = 2 + r.nextDouble()*16; // 1- 18
								double drift = Math.sqrt(aToDRatio*4*diffusioncoefficient/tracklength);
								AbstractSimulator sim1 = new ActiveTransportSimulator(drift, angleVelocity, timelag, dimension, numberOfSteps);
								AbstractSimulator sim2 = new FreeDiffusionSimulator(diffusioncoefficient, timelag, dimension, numberOfSteps);
								sim = new CombinedSimulator(sim1, sim2);
								sigmaPosNoise = Math.sqrt(diffusioncoefficient*timelag + drift*drift*timelag*timelag)/diffusionToNoiseRatio; 
								break;
							case ANOMALOUS:
								sim = new AnomalousDiffusionWMSimulation(diffusioncoefficient, timelag, dimension, 2000, alpha);
								sigmaPosNoise = Math.sqrt(diffusioncoefficient*timelag)/diffusionToNoiseRatio;
								break;
			
							default:
								break;
							}
							t = sim.generateTrajectory();
							if(type==SIM_TYPE.ANOMALOUS){
								t = t.subList(0, numberOfSteps+1);
							}
							
							
							
							t = SimulationUtil.addPositionNoise(t, sigmaPosNoise);
					
							t.setType(typestring);
							trajectorys.add(t);
							tCounter++;
							if(tCounter%10 == 0){
								System.out.println("T: "+ tCounter + " Type " + t.getType());
							}
					}
				}

	
		System.out.println("Tracks generated");
		ExportTools.exportTrajectoriesAsRData(trajectorys, path, timelag);
		trajectorys=null;
		System.out.println("Export done");
		
		
		System.out.println("Done!");

	}

}
