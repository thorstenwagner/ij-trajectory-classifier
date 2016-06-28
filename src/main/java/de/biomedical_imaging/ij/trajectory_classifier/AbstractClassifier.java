package de.biomedical_imaging.ij.trajectory_classifier;

import java.util.ArrayList;

import de.biomedical_imaging.traJ.Trajectory;

public abstract class AbstractClassifier {
	
	public abstract String classify(Trajectory t);
	
	public abstract void start();
	
	public abstract void stop();

	public abstract String[] classify(ArrayList<Trajectory> t);
	
	public abstract double[] getConfindence();
}
