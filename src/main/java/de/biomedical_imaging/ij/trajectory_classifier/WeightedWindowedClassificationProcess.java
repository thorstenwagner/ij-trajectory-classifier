package de.biomedical_imaging.ij.trajectory_classifier;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import org.apache.commons.math3.stat.Frequency;
import de.biomedical_imaging.traJ.Trajectory;

public class WeightedWindowedClassificationProcess {
	
	public String[] windowedClassification(Trajectory t, AbstractClassifier c, int n){
		int windowsize = 2*n+1;
		String[] types = new String[t.size()];
		for(int i = 0; i < types.length; i++){
			types[i] = "NONE";
		}
		int increment = 1;
		ArrayList<Trajectory> tracks = new ArrayList<Trajectory>();
		for(int i = 0; i < (t.size()-windowsize+increment);i=i+increment){
			Trajectory sub = t.subList(i, i+windowsize-1);
			tracks.add(sub);
		}
		System.out.println("Size: " + tracks.size());
		
		String[] res = c.classify(tracks);
		double[] confidence = c.getConfindence();
		
		//Build mapping
		ArrayList<String> restypes = new ArrayList<String>();
		for(int i =0; i < res.length; i++){
			restypes.add(res[i]);
		}
		restypes.add("NONE");
		HashSet<String> uniqueTypes = new HashSet<String>();
		uniqueTypes.addAll(restypes);
		HashMap<String, Integer> mapTypeToInt = new HashMap<String, Integer>();
		HashMap<Integer, String> mapIntToType = new HashMap<Integer, String>();
		Iterator<String> it = uniqueTypes.iterator();
		int key = 0;
		while(it.hasNext()){
			String type = it.next();
			mapTypeToInt.put(type, key);
			mapIntToType.put(key, type);
			key++;
		}
		
		ArrayList<ArrayList<Integer>> weightes = new  ArrayList<ArrayList<Integer>>();
		for(int i = 0; i < t.size(); i++){
			weightes.add(new ArrayList<Integer>());
		}
		
		for(int i = 0; i < res.length; i++){
			for(int j = i; j < i+2*n+1;j++){
				int typ = mapTypeToInt.get(res[i]);
				int weight = (int)(confidence[i]*100);
					for(int k = 0; k < weight; k++){
						weightes.get(j).add(typ);
						
					}
			}
		}
		
		for(int i = 0; i < types.length; i++){
			if(weightes.get(i).size()>0){
				double[] help = arrListToArray(weightes.get(i));
				double[] modes = getMode(help,0,help.length);
				int mode1 = (int)modes[0];
				String mode = mapIntToType.get(mode1);
				types[i] = mode;
		
			}
			else{
				types[i] = "NONE";
			}
		}
		
		return types;
	}
	
	private double[] getMode(double[] values, final int begin, final int length) {
        // Add the values to the frequency table
        Frequency freq = new Frequency();
        for (int i = begin; i < begin + length; i++) {
            final double value = values[i];
            if (!Double.isNaN(value)) {
                freq.addValue(Double.valueOf(value));
            }
        }
        List<Comparable<?>> list = freq.getMode();
        // Convert the list to an array of primitive double
        double[] modes = new double[2];
        modes[0] = ((Double) list.get(0)).doubleValue();
        modes[1] = freq.getCount(Double.valueOf(modes[0]))/(1.0*length);
        return modes;
    }
	
	public double[] arrListToArray(ArrayList<Integer> l){
		double[] a = new double[l.size()];
		for(int i = 0; i < l.size(); i++){
			a[i] = l.get(i).intValue();
		}
		return a;
	}

}
