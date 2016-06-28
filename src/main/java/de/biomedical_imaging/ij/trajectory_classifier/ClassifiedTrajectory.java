package de.biomedical_imaging.ij.trajectory_classifier;

import de.biomedical_imaging.traJ.Trajectory;

public class ClassifiedTrajectory extends Trajectory {
	private static final long serialVersionUID = 3845396333189368623L;
	private long parentID;
	
	public ClassifiedTrajectory(int dimension) {
		super(dimension);
		parentID = 0;
	}
	
	public ClassifiedTrajectory(int dimension, long parentID) {
		super(dimension);
		this.parentID = parentID;
	}
	
	public long getParentID(){
		return parentID;
	}
	
	public void setParentID(long id){
		parentID = id;
	}

}
