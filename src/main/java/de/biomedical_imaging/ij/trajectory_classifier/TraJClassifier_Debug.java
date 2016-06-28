package de.biomedical_imaging.ij.trajectory_classifier;

import ij.measure.CurveFitter;

import java.awt.Color;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math.stat.descriptive.moment.Mean;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math.stat.descriptive.rank.Max;
import org.apache.commons.math.stat.descriptive.rank.Min;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.knowm.xchart.Chart;
import org.knowm.xchart.Series;
import org.knowm.xchart.SeriesMarker;
import org.knowm.xchart.SwingWrapper;

import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.TrajectoryUtil;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.CovarianceDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.RegressionDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.features.ConfinedDiffusionParametersFeature;
import de.biomedical_imaging.traJ.features.MaxDistanceBetweenTwoPositionsFeature;
import de.biomedical_imaging.traJ.features.MeanSquaredDisplacmentFeature;
import de.biomedical_imaging.traJ.features.PowerLawFeature;
import de.biomedical_imaging.traJ.simulation.AbstractSimulator;
import de.biomedical_imaging.traJ.simulation.ActiveTransportSimulator;
import de.biomedical_imaging.traJ.simulation.AnomalousDiffusionWMSimulation;
import de.biomedical_imaging.traJ.simulation.CentralRandomNumberGenerator;
import de.biomedical_imaging.traJ.simulation.ConfinedDiffusionSimulator;
import de.biomedical_imaging.traJ.simulation.FreeDiffusionSimulator;
import de.biomedical_imaging.traJ.simulation.SimulationUtil;
import de.biomedical_imaging.traj.math.ConfinedDiffusionMSDCurveFit.FitMethod;

public class TraJClassifier_Debug {

	public static void main(String[] args) {
		CentralRandomNumberGenerator.getInstance().setSeed(8);
		double diffusioncoefficient = 9.02*Math.pow(10,-2); //[µm^2/s];
		double timelag = 1.0/30;
		double[] driftspeed = {0, 0.27,1,2.4}; // µm/s
		double angleVelocity = Math.PI/4; //rad/s
		int simtracklength = 500;
		int diffusionToNoiseRatio = 2;
		double sigmaPosNoise = Math.sqrt(2*diffusioncoefficient*timelag)/diffusionToNoiseRatio; 
		
		/*
		FreeDiffusionSimulator freesim = new FreeDiffusionSimulator(diffusioncoefficient, timelag, 2, simtracklength);
		
		ActiveTransportSimulator actsim = new ActiveTransportSimulator(driftspeed[2], angleVelocity, timelag, 2, simtracklength);
		
		AnomalousDiffusionWMSimulation anomsim = new AnomalousDiffusionWMSimulation(diffusioncoefficient, timelag, 2, simtracklength, 0.5);
		
		ArrayList<Trajectory> tracks = new ArrayList<Trajectory>();
		for(int i = 0; i < 100; i++){
			tracks.add(anomsim.generateTrajectory());
		}
		
		String[] res = null;

		String[] trueclasses = new String[6*simtracklength+1];
		int c = 0;
			Trajectory t1 = TrajectoryUtil.concactTrajectorie(anomsim.generateTrajectory(), actsim.generateTrajectory());
			
			for(int i = 0; i < simtracklength+1; i++){trueclasses[c]="SUB"; c++;}
			for(int i = 0; i < simtracklength; i++){trueclasses[c]="DIRECTED"; c++;}
			Trajectory t2 = TrajectoryUtil.concactTrajectorie(t1, anomsim.generateTrajectory());
			for(int i = 0; i < simtracklength; i++){trueclasses[c]="SUB"; c++;}
			Trajectory t3 = TrajectoryUtil.concactTrajectorie(t2, actsim.generateTrajectory());
			for(int i = 0; i < simtracklength; i++){trueclasses[c]="DIRECTED"; c++;}
			Trajectory t4 = TrajectoryUtil.concactTrajectorie(t3, freesim.generateTrajectory());
			for(int i = 0; i < simtracklength; i++){trueclasses[c]="FREE"; c++;}
			actsim = new ActiveTransportSimulator(driftspeed[2], angleVelocity, 1.0/30, 2, simtracklength);
			freesim = new FreeDiffusionSimulator(diffusioncoefficient*0.5, 1.0/30, 2, simtracklength);
			
			Trajectory t5 = TrajectoryUtil.combineTrajectory(freesim.generateTrajectory(), actsim.generateTrajectory());
	
			Trajectory t6 = TrajectoryUtil.concactTrajectorie(t4, t5);
			
			t6 = SimulationUtil.addPositionNoise(t6, sigmaPosNoise);
			for(int i = 0; i < simtracklength; i++){trueclasses[c]="DIRECTED"; c++;}
	
			// SUB, ACTIVE, SUB, ACTIVE, FREE, FLOW
			
			
			
			System.out.println("Size Track: " + t6.size());
			
			WeightedWindowedClassificationProcess wp = new WeightedWindowedClassificationProcess();
			AbstractClassifier pred = new RRFClassifierRenjin("/home/thorsten/randomForestModel.RData");
			pred.start();
			res = wp.windowedClassification(t6, pred, 30); 
			System.out.println("Res size " +  res.length);
			System.out.println("trueclasses size " +  trueclasses.length);
			printTypes(res);
			showTrack(t6,res);
			showTrack(t6,trueclasses);
			pred.stop();
			System.out.println("Complete");
			*/
		
		
		ExportImportTools io = new ExportImportTools();
		ArrayList<Trajectory> tracks = io.importTrajectoryDataFromCSV("/home/thorsten/myclassifiedtracks.csv");
		int id = 47364;
		Trajectory t = null;
		ConfinedDiffusionSimulator consim = new ConfinedDiffusionSimulator(diffusioncoefficient, timelag, 0.5, 2, 500);
		t = consim.generateTrajectory();
		/*
		for(int i = 0; i < tracks.size(); i++){
			if(tracks.get(i).getID()==id){
				t = tracks.get(i);
				break;
			}
		}
		
		*/
		System.out.println("Size: " + t.size());
		t.showTrajectory();
		RegressionDiffusionCoefficientEstimator.plotMSDLine(t, 1, t.size()/5);
		MaxDistanceBetweenTwoPositionsFeature maxdist = new MaxDistanceBetweenTwoPositionsFeature(t);
		RegressionDiffusionCoefficientEstimator regest = new RegressionDiffusionCoefficientEstimator(t, 1/timelag, 1, 2);
		
		CovarianceDiffusionCoefficientEstimator covest = new CovarianceDiffusionCoefficientEstimator(t, 1/timelag);
		
		System.out.println("D (REG): " + regest.evaluate()[0]);
		System.out.println("D (COV): " + covest.evaluate()[0]);
		ConfinedDiffusionParametersFeature conv = new ConfinedDiffusionParametersFeature(t, 1/timelag,covest,FitMethod.SIMPLEX);
		double[] convres = conv.evaluate();
		System.out.println("R: " + Math.sqrt(convres[0])/Math.sqrt(Math.PI));//+ "A1: " + conv.evaluate()[1] + " A2: " + conv.evaluate()[2]);
		System.out.println("D" + convres[1]);
		System.out.println("Max. Dist: " + maxdist.evaluate()[0]);
		
		
	//	Trajectory t = consim.generateTrajectory();
	
	}
	

	
	
	public static double[] evaluate(Trajectory t, double timelag) {
		MeanSquaredDisplacmentFeature msd = new MeanSquaredDisplacmentFeature(t, 1);
		msd.setOverlap(false);

		ArrayList<Double> xDataList = new ArrayList<Double>();
		ArrayList<Double> yDataList = new ArrayList<Double>();
		for(int i = 1; i < t.size(); i++){
			msd.setTimelag(i);
			double[] res = msd.evaluate();
			double msdvalue = res[0];
			int N = (int)res[2];
			for(int j = 0; j < N; j++){
				xDataList.add((double) i*timelag);
				yDataList.add(msdvalue);
			}
		}
		double[] xData = ArrayUtils.toPrimitive(xDataList.toArray(new Double[0]));
		double[] yData = ArrayUtils.toPrimitive(yDataList.toArray(new Double[0]));
		CurveFitter fitter = new CurveFitter(xData, yData);
		MaxDistanceBetweenTwoPositionsFeature maxdist = new MaxDistanceBetweenTwoPositionsFeature(t);
		RegressionDiffusionCoefficientEstimator regest = new RegressionDiffusionCoefficientEstimator(t, 1/timelag, 5, 5);
		double estrad = maxdist.evaluate()[0];
		double estDC = regest.evaluate()[0];
		double[] initialParams = {estrad*estrad};//,regest.evaluate()[0]};
		//fitter.doCustomFit("y=a*(1-b*exp(-4*c*"+estDC+"*x/a))", initialParams, false);
		fitter.doCustomFit("y=a*(1-exp(-4*"+estDC+"*x/a))", initialParams, false);
		double[] params = fitter.getParams();
		double[] res = {params[0]};//params[1],params[2],fitter.getFitGoodness()};
		return res;
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
				    Color.BLACK, Color.PINK, Color.ORANGE, Color.MAGENTA
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
	
	

}
