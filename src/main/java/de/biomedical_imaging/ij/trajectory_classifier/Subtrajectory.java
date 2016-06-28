package de.biomedical_imaging.ij.trajectory_classifier;

import de.biomedical_imaging.traJ.Trajectory;

public class Subtrajectory extends Trajectory {
	
	private Trajectory parent;
	
	public Subtrajectory(Trajectory parent, int dimension) {
		super(dimension);
		this.parent = parent;
	}
	
	public Subtrajectory(Trajectory parent, int dimension, int relativeStartPoint) {
		super(dimension,relativeStartPoint);
		this.parent = parent;
	}
	
	public Trajectory getParent(){
		return parent;
	}

}
