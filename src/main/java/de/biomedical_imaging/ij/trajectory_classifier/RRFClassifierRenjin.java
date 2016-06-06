package de.biomedical_imaging.ij.trajectory_classifier;

import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

import org.renjin.parser.ParseException;
import org.renjin.sexp.StringVector;

import de.biomedical_imaging.ij.trajectory_classifier.FeatureWorker.EVALTYPE;
import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.features.Asymmetry2Feature;
import de.biomedical_imaging.traJ.features.Asymmetry3Feature;
import de.biomedical_imaging.traJ.features.AsymmetryFeature;
import de.biomedical_imaging.traJ.features.EfficiencyFeature;
import de.biomedical_imaging.traJ.features.ElongationFeature;
import de.biomedical_imaging.traJ.features.FractalDimensionFeature;
import de.biomedical_imaging.traJ.features.GaussianityFeauture;
import de.biomedical_imaging.traJ.features.KurtosisFeature;
import de.biomedical_imaging.traJ.features.MSDRatioFeature;
import de.biomedical_imaging.traJ.features.MeanSquaredDisplacmentCurvature;
import de.biomedical_imaging.traJ.features.PowerLawFeature;
import de.biomedical_imaging.traJ.features.ShortTimeLongTimeDiffusioncoefficentRatio;
import de.biomedical_imaging.traJ.features.SkewnessFeature;
import de.biomedical_imaging.traJ.features.SplineCurveDynamicsFeature;
import de.biomedical_imaging.traJ.features.StraightnessFeature;
import de.biomedical_imaging.traJ.features.TrappedProbabilityFeature;

public class RRFClassifierRenjin extends AbstractClassifier  {

	private ScriptEngine engine = null;
	private boolean chatty = false;
	@Override
	public String classify(Trajectory t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void start() {
		ScriptEngineManager manager = new ScriptEngineManager();
	    // create a Renjin engine:
	    engine = manager.getEngineByName("Renjin");
	    try {
			engine.eval("library(randomForest)");
			engine.eval("library(plyr)");
			engine.eval("load(\"/home/thorsten/randomForestModel.RData\")");
		} catch (ScriptException e) {
			e.printStackTrace();
		}
		
	    // check if the engine has loaded correctly:
	    if(engine == null) {
	        throw new RuntimeException("Renjin Script Engine not found on the classpath.");
	    }
	    
		
	}

	@Override
	public void stop() {
		engine =null;
		
	}

	@Override
	public String[] classify(ArrayList<Trajectory> tracks) {
		
		int N = tracks.size();
		String[] result = new String[N];
		double[] elong = new double[N];
		double[] fd = new double[N];
		double[] msdcurvature = new double[N];
		double[] power = new double[N];
		double[] sdDir = new double[N]; 
		double[] dratio = new double[N]; 
		double[] ltStRatio = new double[N]; 
		double[] asym1 = new double[N]; 
		double[] asym2 = new double[N];
		double[] asym3 = new double[N];
		double[] efficiency = new double[N];
		double[] kurtosis = new double[N];
		double[] skewness = new double[N];
		double[] msdratio = new double[N];
		double[] straightness = new double[N];
		double[] trappedness = new double[N];
		double[] gaussianity = new double[N];
		
		int numberOfSegmentsSplineFit = 7;
		int numberOfPointsForShortTimeLongTimeRatio = 2;
		
		
		int cores = Runtime.getRuntime().availableProcessors();
		ExecutorService pool = Executors.newFixedThreadPool(cores);
		
		for(int i = 0; i < tracks.size(); i++){
			Trajectory t = tracks.get(i);

			ElongationFeature elongF = new ElongationFeature(t);
			pool.submit(new FeatureWorker(elong, i,elongF, EVALTYPE.FIRST));
			
			FractalDimensionFeature fdF = new FractalDimensionFeature(t);
			pool.submit(new FeatureWorker(fd, i,fdF, EVALTYPE.FIRST));
			
			MeanSquaredDisplacmentCurvature msdCurv = new MeanSquaredDisplacmentCurvature(t);
			pool.submit(new FeatureWorker(msdcurvature, i,msdCurv, EVALTYPE.FIRST));
			//if(chatty)System.out.println("MSDCURV evaluated");
			
			PowerLawFeature pwf = new PowerLawFeature(t, 1, t.size()-1);
			pool.submit(new FeatureWorker(power, i,pwf, EVALTYPE.FIRST));
			if(chatty)System.out.println("POWER evaluated");

		//	StandardDeviationDirectionFeature sdf = new StandardDeviationDirectionFeature(t, timelagForDirectionDeviationLong);
		//	pool.submit(new FeatureWorker(sdDir, i,sdf, EVALTYPE.FIRST));
		//	if(chatty)System.out.println("SDDIR evaluated");

			SplineCurveDynamicsFeature scdf = new SplineCurveDynamicsFeature(t, numberOfSegmentsSplineFit, 1);
			pool.submit(new FeatureWorker(dratio, i,scdf, EVALTYPE.RATIO_12));
			if(chatty)System.out.println("SCDF evaluated");

			
			
			AsymmetryFeature asymf1 = new AsymmetryFeature(t);
			pool.submit(new FeatureWorker(asym1, i,asymf1, EVALTYPE.FIRST));
		//	if(chatty)System.out.println("ASYM1 evaluated");
			
			Asymmetry2Feature asymf2 = new Asymmetry2Feature(t);
			pool.submit(new FeatureWorker(asym2, i,asymf2, EVALTYPE.FIRST));
			if(chatty)System.out.println("ASYMf2 evaluated");
			
			Asymmetry3Feature asymf3 = new Asymmetry3Feature(t);
			pool.submit(new FeatureWorker(asym3, i,asymf3, EVALTYPE.FIRST));
			if(chatty)System.out.println("ASYMf3 evaluated");
			
			EfficiencyFeature eff = new EfficiencyFeature(t);
			pool.submit(new FeatureWorker(efficiency, i,eff, EVALTYPE.FIRST));
			if(chatty)System.out.println("EFF evaluated");
			
			ShortTimeLongTimeDiffusioncoefficentRatio stltdf = new ShortTimeLongTimeDiffusioncoefficentRatio(t, numberOfPointsForShortTimeLongTimeRatio);
			pool.submit(new FeatureWorker(ltStRatio, i,stltdf, EVALTYPE.FIRST));
			if(chatty)System.out.println("STLTDF evaluated");
			
			KurtosisFeature kurtf = new KurtosisFeature(t);
			pool.submit(new FeatureWorker(kurtosis, i,kurtf, EVALTYPE.FIRST));
		//	if(chatty)System.out.println("KURT evaluated");
			
			SkewnessFeature skew = new SkewnessFeature(t);
			pool.submit(new FeatureWorker(skewness, i,skew, EVALTYPE.FIRST));
		//	if(chatty)System.out.println("SKEW evaluated");

			MSDRatioFeature msdr = new MSDRatioFeature(t, 1,10);
			pool.submit(new FeatureWorker(msdratio, i,msdr, EVALTYPE.FIRST));
			//if(chatty)System.out.println("MSDR evaluated");
			
			StraightnessFeature straight = new StraightnessFeature(t);
			pool.submit(new FeatureWorker(straightness, i,straight, EVALTYPE.FIRST));
		//	if(chatty)System.out.println("STRAIGHT evaluated");
			
			TrappedProbabilityFeature trappf = new TrappedProbabilityFeature(t);
			pool.submit(new FeatureWorker(trappedness, i,trappf, EVALTYPE.FIRST));
			if(chatty)System.out.println("TRAPP evaluated");

			GaussianityFeauture gaussf = new GaussianityFeauture(t, 1);
			pool.submit(new FeatureWorker(gaussianity, i,gaussf, EVALTYPE.FIRST));
			if(chatty)System.out.println("GAUSS evaluated");
		}
		
		pool.shutdown();
		try {
			pool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
			} catch (InterruptedException e) {
			  e.printStackTrace();
			}
		
		try {
			
			engine.put("elong", elong);
			engine.put("fd",fd);
			engine.put("msdcurvature",msdcurvature);
			engine.put("power", power);
		//	engine.put("sdDir", sdDir);
			engine.put("D.ratio", dratio);
			engine.put("LtStRatio", ltStRatio);
			engine.put("asymmetry1", asym1);
			engine.put("asymmetry2", asym2);
			engine.put("asymmetry3", asym3);
			engine.put("efficiency", efficiency);
			engine.put("kurtosis",kurtosis);
			engine.put("skewness", skewness);
			engine.put("msdratio", msdratio);
			engine.put("straightness", straightness);
			engine.put("trappedness", trappedness);
			engine.put("gaussianity", gaussianity);
			
			engine.eval("data<-data.frame(ELONG=elong,FD=fd,MSD.C=msdcurvature,"
					+ "POWER=power,SPLINE.RATIO=D.ratio,LTST.RATIO=LtStRatio,"
					+ "ASYM1=asymmetry1,MSD.R=msdratio,ASYM2=asymmetry2,ASYM3=asymmetry3,EFFICENCY=efficiency, KURT=kurtosis,"
					+ "SKEW=skewness,STRAIGHTNESS=straightness, "
					+ "TRAPPED=trappedness,GAUSS=gaussianity)");
			engine.eval("features.predict <- predict(model,data)");
			engine.eval("lvl <- levels(model$y)");
			
			StringVector res = (StringVector)engine.eval("lvl[features.predict]");
			
			result = res.toArray();
		}
		catch (ParseException e) {
		    System.out.println("R script parse error: " + e.getMessage());
		} 
		catch (ScriptException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return result;
	}

}
