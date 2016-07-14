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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

import org.renjin.parser.ParseException;
import org.renjin.sexp.DoubleVector;
import org.renjin.sexp.StringVector;

import de.biomedical_imaging.ij.trajectory_classifier.FeatureWorker.EVALTYPE;
import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.RegressionDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.features.Asymmetry3Feature;
import de.biomedical_imaging.traJ.features.EfficiencyFeature;
import de.biomedical_imaging.traJ.features.FractalDimensionFeature;
import de.biomedical_imaging.traJ.features.GaussianityFeauture;
import de.biomedical_imaging.traJ.features.KurtosisFeature;
import de.biomedical_imaging.traJ.features.MSDRatioFeature;
import de.biomedical_imaging.traJ.features.PowerLawFeature;
import de.biomedical_imaging.traJ.features.ShortTimeLongTimeDiffusioncoefficentRatio;
import de.biomedical_imaging.traJ.features.SkewnessFeature;
import de.biomedical_imaging.traJ.features.StraightnessFeature;
import de.biomedical_imaging.traJ.features.TrappedProbabilityFeature;

public class RRFClassifierRenjin extends AbstractClassifier  {

	private ScriptEngine engine = null;
	private String pathToModel;
	private double[] confindence;
	private double timelag;
	public RRFClassifierRenjin(String pathToModel,double timelag) {
		this.pathToModel = pathToModel;
		this.timelag = timelag;
	}
	
	public void setTimelag(double timelag){
		this.timelag = timelag;
	}
	
	@Override
	public String classify(Trajectory t)  {
		ArrayList<Trajectory> tracks = new ArrayList<Trajectory>();
		tracks.add(t);
		
		return classify(tracks)[0];
	}

	@Override
	public void start() {
		ScriptEngineManager manager = new ScriptEngineManager();
	    // create a Renjin engine:
	    engine = manager.getEngineByName("Renjin");
	    try {
			engine.eval("library(randomForest)");
			engine.eval("library(plyr)");
			engine.eval("load(\""+pathToModel+"\")");
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
	public String[] classify(ArrayList<Trajectory> tracks)  {
		
		int N = tracks.size();
		String[] result = new String[N];
		double[] fd = new double[N];
		int[] lengths = new int[N];
		double[] power = new double[N];
		Arrays.fill(power, -1);
		//double[] ltStRatio = new double[N]; 
		double[] asym3 = new double[N];
		double[] efficiency = new double[N];
		double[] kurtosis = new double[N];
		double[] skewness = new double[N];
		double[] msdratio = new double[N];
		double[] straightness = new double[N];
		double[] trappedness = new double[N];
		double[] gaussianity = new double[N];
		double[] pwrDCs = new double[N];
		Arrays.fill(power, -1);
		int numberOfPointsForShortTimeLongTimeRatio = 2;
		double start= System.currentTimeMillis();
		int cores = Runtime.getRuntime().availableProcessors();
		ExecutorService pool = Executors.newFixedThreadPool(cores);
		
		
		for(int i = 0; i < tracks.size(); i++){
			Trajectory t = tracks.get(i);
			
			lengths[i] = t.size();

			
			FractalDimensionFeature fdF = new FractalDimensionFeature(t);
			pool.submit(new FeatureWorker(fd, i,fdF, EVALTYPE.FIRST));
			double initDC=0;
			double initAlpha =0;
			if(i-1>0 && power[i-1]>0 && pwrDCs[i-1]>0){
				initDC = pwrDCs[i-1];
				initAlpha = power[i-1];
	
			}else{
				RegressionDiffusionCoefficientEstimator regest = new RegressionDiffusionCoefficientEstimator(t, 1.0/timelag, 1, 3);
				initDC= regest.evaluate()[0];
				initAlpha = 0.5;
			}
			
			PowerLawFeature pwf = new PowerLawFeature(t, 1, t.size()/3,initAlpha,initDC);
			pool.submit(new FeatureWorker(power, i,pwf, EVALTYPE.FIRST));
			pool.submit(new FeatureWorker(pwrDCs, i,pwf, EVALTYPE.SECOND));
			
			Asymmetry3Feature asymf3 = new Asymmetry3Feature(t);
			pool.submit(new FeatureWorker(asym3, i,asymf3, EVALTYPE.FIRST));
		
			EfficiencyFeature eff = new EfficiencyFeature(t);
			pool.submit(new FeatureWorker(efficiency, i,eff, EVALTYPE.FIRST));
			
		//	ShortTimeLongTimeDiffusioncoefficentRatio stltdf = new ShortTimeLongTimeDiffusioncoefficentRatio(t, numberOfPointsForShortTimeLongTimeRatio);
		//	pool.submit(new FeatureWorker(ltStRatio, i,stltdf, EVALTYPE.FIRST));
			
			KurtosisFeature kurtf = new KurtosisFeature(t);
			pool.submit(new FeatureWorker(kurtosis, i,kurtf, EVALTYPE.FIRST));
			
			SkewnessFeature skew = new SkewnessFeature(t);
			pool.submit(new FeatureWorker(skewness, i,skew, EVALTYPE.FIRST));

			MSDRatioFeature msdr = new MSDRatioFeature(t, 1,5);
			pool.submit(new FeatureWorker(msdratio, i,msdr, EVALTYPE.FIRST));
			
			StraightnessFeature straight = new StraightnessFeature(t);
			pool.submit(new FeatureWorker(straightness, i,straight, EVALTYPE.FIRST));
			
			TrappedProbabilityFeature trappf = new TrappedProbabilityFeature(t);
			pool.submit(new FeatureWorker(trappedness, i,trappf, EVALTYPE.FIRST));

			GaussianityFeauture gaussf = new GaussianityFeauture(t, 1);
			pool.submit(new FeatureWorker(gaussianity, i,gaussf, EVALTYPE.FIRST));
	
		}
		pool.shutdown();
	
		try {
			pool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
			} catch (InterruptedException e) {
			  e.printStackTrace();
			}

		try {
			
			engine.put("fd",fd);
			engine.put("lengths",lengths);
			engine.put("power", power);
		//	engine.put("LtStRatio", ltStRatio);
			engine.put("asymmetry3", asym3);
			engine.put("efficiency", efficiency);
			engine.put("kurtosis",kurtosis);
			engine.put("skewness", skewness);
			engine.put("msdratio", msdratio);
			engine.put("straightness", straightness);
			engine.put("trappedness", trappedness);
			engine.put("gaussianity", gaussianity);
			
			engine.eval("data<-data.frame(LENGTHS=lengths,FD=fd,"
					+ "POWER=power,"
					+ "MSD.R=msdratio,ASYM3=asymmetry3,EFFICENCY=efficiency, KURT=kurtosis,"
					+ "SKEW=skewness,STRAIGHTNESS=straightness, "
					+ "TRAPPED=trappedness,GAUSS=gaussianity)");

			engine.eval("features.predict <- predict(model,data,type=\"prob\")");
			engine.eval("fprob<-features.predict");
		
			if(tracks.size()>1){
				engine.eval("probs <- as.vector(apply(fprob[1:nrow(fprob),],1,max))");
				engine.eval("indexmax <- as.vector(apply(fprob[1:nrow(fprob),],1,which.max))");
			}
			else{
				engine.eval("probs <- max(fprob)");
				engine.eval("indexmax <- which.max(fprob)");
			}
			engine.eval("mynames <- colnames(fprob)");
			engine.eval("maxclass <- mynames[indexmax]");
			StringVector res = (StringVector)engine.eval("maxclass");
			result = res.toArray();
			DoubleVector confi = (DoubleVector)engine.eval("probs");
			confindence = confi.toDoubleArray();
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

	@Override
	public double[] getConfindence() {
		// TODO Auto-generated method stub
		return confindence;
	}

}
