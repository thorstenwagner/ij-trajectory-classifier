package de.biomedical_imaging.ij.trajectory_classifier;

import java.awt.Color;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.ArrayList;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.knowm.xchart.Chart;
import org.knowm.xchart.QuickChart;
import org.knowm.xchart.Series;
import org.knowm.xchart.SeriesMarker;
import org.knowm.xchart.SwingWrapper;

import sun.tools.tree.PreIncExpression;
import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.TrajectoryUtil;
import de.biomedical_imaging.traJ.simulation.ActiveTransportSimulator;
import de.biomedical_imaging.traJ.simulation.AnomalousDiffusionWMSimulation;
import de.biomedical_imaging.traJ.simulation.CentralRandomNumberGenerator;
import de.biomedical_imaging.traJ.simulation.ConfinedDiffusionSimulator;
import de.biomedical_imaging.traJ.simulation.FreeDiffusionSimulator;

public class TraJClassifier_Debug {

	public static void main(String[] args) {
		CentralRandomNumberGenerator.getInstance().setSeed(6);
		double diffusioncoefficient = 9.02*Math.pow(10,-2); //[µm^2/s];
		double[] driftspeed = {0, 0.27,0.8,2.4}; // µm/s
		double angleVelocity = Math.PI/4.0; //rad/s
		
		FreeDiffusionSimulator freesim = new FreeDiffusionSimulator(diffusioncoefficient, 1.0/30, 2, 500);
		
		ActiveTransportSimulator actsim = new ActiveTransportSimulator(driftspeed[2], angleVelocity, 1.0/30, 2, 500);
		
		AnomalousDiffusionWMSimulation anomsim = new AnomalousDiffusionWMSimulation(diffusioncoefficient, 1.0/30, 2, 500, 0.5);
		
		ArrayList<Trajectory> tracks = new ArrayList<Trajectory>();
		for(int i = 0; i < 100; i++){
			tracks.add(anomsim.generateTrajectory());
		}
		RRFClassifier pred = new RRFClassifier(1.0/30);
		String[] res = null;
		
		//pred.classify(tracks);
		//for(int i = 0; i < res.length; i++){
		//	System.out.println("Res: " + res[i]);
		//}
		
		int windowLength = 100;
		ArrayList<String> types = new ArrayList<String>();
	
			Trajectory t1 = TrajectoryUtil.concactTrajectorie(anomsim.generateTrajectory(), actsim.generateTrajectory());
			Trajectory t2 = TrajectoryUtil.concactTrajectorie(t1, anomsim.generateTrajectory());
			Trajectory t3 = TrajectoryUtil.concactTrajectorie(t2, actsim.generateTrajectory());
			Trajectory t4 = TrajectoryUtil.concactTrajectorie(t3, freesim.generateTrajectory());
			actsim = new ActiveTransportSimulator(driftspeed[2], angleVelocity, 1.0/30, 2, 500);
			freesim = new FreeDiffusionSimulator(diffusioncoefficient, 1.0/30, 2, 500);
			
			Trajectory t5 = TrajectoryUtil.combineTrajectory(freesim.generateTrajectory(), actsim.generateTrajectory());
	
	
			Trajectory t6 = TrajectoryUtil.concactTrajectorie(t4, t5);
			t6.showTrajectory();
			// SUB, ACTIVE, SUB, ACTIVE, FREE, FLOW
			pred = new RRFClassifier(1.0/30);
			
			
			System.out.println("Size Track: " + t6.size());
			res = windowedClassification(t6, pred, 80);
			res = movingMode(res, 40);
			System.out.println("Res size " +  res.length);
			printTypes(res);
			showTrack(t6,res);
			
			//res = movingMode(resList, uniqueTypes, 10);
			
			System.out.println("Complete");
	}
	
	public static String[] windowedClassification(Trajectory t, AbstractClassifier c, int n){
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
	
	public static void showTrack(Trajectory t, String[] classes){
		if(t.size() != classes.length){
			throw new IllegalArgumentException("Tracklength and the number of classes have to be equal");
		}
		ArrayList<String> cls = new ArrayList<String>();
		ArrayList<Trajectory> subtracks = new ArrayList<Trajectory>();
		String prevCls = classes[0];
		Trajectory sub = new Trajectory(t.getDimension());
		sub.add(t.get(0));
		for(int i = 1; i < t.size(); i++){
			if(prevCls == classes[i]){
				sub.add(t.get(i));
			}else{
				subtracks.add(sub);
				cls.add(prevCls);
				prevCls = classes[i];
				sub = new Trajectory(t.getDimension());
				sub.add(t.get(i));
			}
		}
		cls.add(prevCls);
		subtracks.add(sub);
		showTracks(subtracks, cls);
	}
	
	public static void showTracks(ArrayList<Trajectory> tracks, String[] classes){
		ArrayList<String> cls = new ArrayList<String>();
		for(int i = 0; i < classes.length; i++){
			cls.add(classes[i]);
		}
		showTracks(tracks, cls);
	}
	
	
	public static void showTracks(ArrayList<Trajectory> tracks, ArrayList<String> classes){
		Color[] colors = new Color[]
				{
				    Color.RED, Color.BLUE, Color.YELLOW, Color.GREEN,
				    Color.BLACK, Color.PINK, Color.ORANGE
				};
		HashSet<String> uniqueClasses = new HashSet<String>();
		uniqueClasses.addAll(classes);
		HashMap<String, Integer> mapTypeToInt = new HashMap<String, Integer>();
		HashMap<Integer, String> mapIntToType = new HashMap<Integer, String>();
		Iterator<String> it = uniqueClasses.iterator();
		int key = 0;
		while(it.hasNext()){
			String type = it.next();
			mapTypeToInt.put(type, key);
			mapIntToType.put(key, type);
			key++;
		}
		
		Chart chart = new Chart(700, 400);
		for(int i = 0; i < tracks.size(); i++){
			Trajectory t = tracks.get(i);
			String type = classes.get(i);
			if(t.getDimension()==2){
			 	double[] xData = new double[t.size()];
			    double[] yData = new double[t.size()];
			    for(int j = 0; j < t.size(); j++){
			    	xData[j] = t.get(j).x;
			    	yData[j] = t.get(j).y;
			    	
			    }
			    // Create Chart
			    Series s = chart.addSeries(i+". "+type+"("+t.size()+")", xData, yData);
			    s.setLineColor(colors[mapTypeToInt.get(type)]);
			    s.setMarker(SeriesMarker.NONE);
		
			    
			   
			} 	
		}
		
		//Show it
		 new SwingWrapper(chart).displayChart();
	}
	
	public static void printTypes(String[] res){
		String prev = res[0];
		int c = 1;
		for(int i = 1; i < res.length; i++){
			if(prev==res[i]){
				c++;
			}else{
				System.out.println("Type: " + prev + "("+c+")");
				c= 1;
				prev = res[i];
			}
			
		}
		System.out.println("Type: " + prev + "("+c+")");
	}
	
	public static void addToTypes(ArrayList<String> types, String type, int n){
		
		for(int i = 0; i < n; i++){
			types.add(type);
		}
	}
	public static String[] movingMode(String[] types, int n){
		ArrayList<String> ltypes = new ArrayList<String>();
		for(int i =0; i < types.length; i++){
			ltypes.add(types[i]);
		}
		return movingMode(ltypes, n);
	}
	
	public static String[] movingMode(ArrayList<String> types, int n){
		int windowsize = 2*n+1;
		HashSet<String> uniqueTypes = new HashSet<String>();
		uniqueTypes.addAll(types);
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

		String[] medTypes = new String[types.size()];
		
		for(int i = 0; i < n; i++){
			medTypes[i] = "NONE";
		}
		for(int i = types.size()- n; i < types.size(); i++){
			medTypes[i] = "NONE";
		}
		
		for(int i = 0; i < (types.size()-windowsize+1);i++){
			List<String> sub = types.subList(i, i+windowsize-1);
			double[] values = new double[sub.size()];
			for(int j = 0; j < sub.size(); j++){
				values[j] = mapTypeToInt.get(sub.get(j));
			}
			
			medTypes[i + n] = mapIntToType.get(((int)StatUtils.mode(values)[0]));
		}
		return medTypes;	
	}
	
	public static String[] movingMedian(ArrayList<String> types, HashSet<String> uniqueTypes, int windowsize){
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
		Median med = new Median();
		
		String[] medTypes = new String[(types.size()-windowsize)];
		for(int i = 0; i < (types.size()-windowsize);i++){
			List<String> sub = types.subList(i, i+windowsize-1);
			double[] values = new double[sub.size()];
			for(int j = 0; j < sub.size(); j++){
				values[j] = mapTypeToInt.get(sub.get(j));
			}
			
			medTypes[i] = mapIntToType.get(((int)med.evaluate(values)));
		}
		return medTypes;	
	}

}
