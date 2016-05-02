package de.biomedical_imaging.ij.trajectory_classifier;

import java.util.ArrayList;

import de.biomedical_imaging.traJ.Trajectory;

public class Track extends Trajectory {

	private static final long serialVersionUID = -3430994540564186417L;
	private ArrayList<Integer> frame;
	
	public Track(int dimension) {
		super(dimension);
	}

	
	
	public boolean add(double x, double y, Integer t){
		frame.add(t);
		return this.add(x, y, 0);
	}
	
	/**
	 * 
	 * @return List of the frame numbers
	 */
	public ArrayList<Integer> getListFrame(){
		return frame;
	}
}
