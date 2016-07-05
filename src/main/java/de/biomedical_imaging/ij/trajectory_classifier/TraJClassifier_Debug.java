/*
MIT License

Copyright (c) 2016 Thorsten Wagner

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

package de.biomedical_imaging.ij.trajectory_classifier;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.CurveFitter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

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

import de.biomedical_imaging.ij.trajectory_classifier.FeatureWorker.EVALTYPE;
import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.TrajectoryUtil;
import de.biomedical_imaging.traJ.VisualizationUtils;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.CovarianceDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.RegressionDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.features.AbstractTrajectoryFeature;
import de.biomedical_imaging.traJ.features.Asymmetry2Feature;
import de.biomedical_imaging.traJ.features.Asymmetry3Feature;
import de.biomedical_imaging.traJ.features.AsymmetryFeature;
import de.biomedical_imaging.traJ.features.ConfinedDiffusionParametersFeature;
import de.biomedical_imaging.traJ.features.EfficiencyFeature;
import de.biomedical_imaging.traJ.features.ElongationFeature;
import de.biomedical_imaging.traJ.features.FractalDimensionFeature;
import de.biomedical_imaging.traJ.features.GaussianityFeauture;
import de.biomedical_imaging.traJ.features.KurtosisFeature;
import de.biomedical_imaging.traJ.features.MSDRatioFeature;
import de.biomedical_imaging.traJ.features.MaxDistanceBetweenTwoPositionsFeature;
import de.biomedical_imaging.traJ.features.MeanSquaredDisplacmentCurvature;
import de.biomedical_imaging.traJ.features.MeanSquaredDisplacmentFeature;
import de.biomedical_imaging.traJ.features.PowerLawFeature;
import de.biomedical_imaging.traJ.features.ShortTimeLongTimeDiffusioncoefficentRatio;
import de.biomedical_imaging.traJ.features.SkewnessFeature;
import de.biomedical_imaging.traJ.features.SplineCurveDynamicsFeature;
import de.biomedical_imaging.traJ.features.StraightnessFeature;
import de.biomedical_imaging.traJ.features.TrappedProbabilityFeature;
import de.biomedical_imaging.traJ.simulation.AbstractSimulator;
import de.biomedical_imaging.traJ.simulation.ActiveTransportSimulator;
import de.biomedical_imaging.traJ.simulation.AnomalousDiffusionWMSimulation;
import de.biomedical_imaging.traJ.simulation.CentralRandomNumberGenerator;
import de.biomedical_imaging.traJ.simulation.ConfinedDiffusionSimulator;
import de.biomedical_imaging.traJ.simulation.FreeDiffusionSimulator;
import de.biomedical_imaging.traJ.simulation.SimulationUtil;
import de.biomedical_imaging.traj.math.PowerLawCurveFit.FitMethod;

public class TraJClassifier_Debug {

	public static void main(String[] args) {
		CentralRandomNumberGenerator.getInstance().setSeed(8);
		double diffusioncoefficient = 50;//9.02*Math.pow(10,-2); //[µm^2/s];
		double timelag = 1.0/30;
		double[] driftspeed = {0, 0.27,1,2.4}; // µm/s
		double angleVelocity = Math.PI/4; //rad/s
		int simtracklength = 500;
		int diffusionToNoiseRatio = 2;
		double sigmaPosNoise = Math.sqrt(2*diffusioncoefficient*timelag)/diffusionToNoiseRatio; 
		
		
		AbstractSimulator sim = new FreeDiffusionSimulator(diffusioncoefficient, timelag, 2, simtracklength);
		ArrayList<Trajectory> t = new ArrayList<Trajectory>();
		Trajectory tr = sim.generateTrajectory();
		//tr.setRelativStartTimepoint(1);
		//tr.scale(1.0/0.166);
		tr.offset(500, 500, 0);
		ArrayList<Chart> lc = new ArrayList<Chart>();
		Chart c = VisualizationUtils.getTrajectoryChart(tr);
		lc.add(c);
		//t.add(tr);
		
		
		sim = new AnomalousDiffusionWMSimulation(diffusioncoefficient, timelag, 2, 2000, 0.5);
		Trajectory tr2 = sim.generateTrajectory();
		tr2 = tr2.subList(0, 500+1);
		t.add(TrajectoryUtil.concactTrajectorie(tr, tr2));
		
		double radius_confined = Math.sqrt(-1*Math.log(0.9)*(4*diffusioncoefficient*60*timelag));
		System.out.println("Radius:" + radius_confined*1000);
		sim = new ConfinedDiffusionSimulator(diffusioncoefficient,timelag,radius_confined,2,500);
		tr = sim.generateTrajectory();
		tr.offset(250, 250, 0);
		t.add(tr);
		new ImageJ();
		IJ.getInstance().show(true);
		ImageStack is = new ImageStack(1000, 1000);
		for(int i = 0; i < 1005; i++){
			is.addSlice(new ByteProcessor(1000, 1000));
		}
		
		ImagePlus img = new ImagePlus("", is);
		img.show();
		
		TraJClassifier_ tclass = new TraJClassifier_();
		tclass.setTracksToClassify(t);
		final long timeStart = System.currentTimeMillis(); 

		tclass.run("DEBUG");
		final long timeEnd = System.currentTimeMillis(); 
	    System.out.println("Laufzeit: " + (timeEnd - timeStart) + " Millisek."); 
		
	
	}
	
	
	public static void outputFeatures(Trajectory t,double timelag){
		int numberOfSegmentsSplineFit = 7;
		int numberOfPointsForShortTimeLongTimeRatio = 2;
		AbstractTrajectoryFeature f = new ElongationFeature(t);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
		
		f  = new FractalDimensionFeature(t);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
		
		f  = new MeanSquaredDisplacmentCurvature(t);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
		
		f  = new RegressionDiffusionCoefficientEstimator(t, 1.0/timelag, 1, 3);
		PowerLawFeature f2 = new PowerLawFeature(t, 1, t.size()/3, FitMethod.SIMPLEX,0.5,f.evaluate()[0]);
		System.out.println(f.getShortName()+": " + f2.evaluate()[0]);

	//	StandardDeviationDirectionFeature sdf = new StandardDeviationDirectionFeature(t, timelagForDirectionDeviationLong);
	//	pool.submit(new FeatureWorker(sdDir, i,sdf, EVALTYPE.FIRST));
	//	if(chatty)System.out.println("SDDIR evaluated");

		f = new SplineCurveDynamicsFeature(t, numberOfSegmentsSplineFit, 1);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]/f.evaluate()[1]);
		
		f = new AsymmetryFeature(t);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);

		
		f = new Asymmetry2Feature(t);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
		
		f = new Asymmetry3Feature(t);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
		
		f = new EfficiencyFeature(t);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
		
		f = new ShortTimeLongTimeDiffusioncoefficentRatio(t, numberOfPointsForShortTimeLongTimeRatio);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
		
		f = new KurtosisFeature(t);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
		
		f = new SkewnessFeature(t);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
		
		f = new MSDRatioFeature(t, 1,5);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
		
		f = new StraightnessFeature(t);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
		
		f = new TrappedProbabilityFeature(t);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
		
		f= new GaussianityFeauture(t, 1);
		System.out.println(f.getShortName()+": " + f.evaluate()[0]);
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
