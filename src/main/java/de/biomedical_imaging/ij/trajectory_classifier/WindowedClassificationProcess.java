package de.biomedical_imaging.ij.trajectory_classifier;

import java.util.ArrayList;

import de.biomedical_imaging.traJ.Trajectory;

public class WindowedClassificationProcess {
	
	public String[] windowedClassification(Trajectory t, AbstractClassifier c, int n){
		int windowsize = 2*n+1;
		String[] types = new String[t.size()];
		for(int i = 0; i < n; i++){
			types[i] = "NONE";
		}
		for(int i = types.length- n; i < types.length; i++){
			types[i] = "NONE";
		}
		ArrayList<Trajectory> tracks = new ArrayList<Trajectory>();
		for(int i = 0; i < (t.size()-windowsize+1);i++){
			Trajectory sub = t.subList(i, i+windowsize-1);
			tracks.add(sub);
		}
		String[] res = c.classify(tracks);
		for(int i = 0; i < res.length; i++){
			types[i + n] = res[i];
		}
		return types;
	}

}
