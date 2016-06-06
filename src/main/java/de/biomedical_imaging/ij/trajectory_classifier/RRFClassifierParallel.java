package de.biomedical_imaging.ij.trajectory_classifier;

import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;

import de.biomedical_imaging.ij.trajectory_classifier.FeatureWorker.EVALTYPE;
import de.biomedical_imaging.ij.trajectory_classifier.r.StartRserve;
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
import de.biomedical_imaging.traJ.features.StandardDeviationDirectionFeature;
import de.biomedical_imaging.traJ.features.StraightnessFeature;
import de.biomedical_imaging.traJ.features.TrappedProbabilityFeature;

public class RRFClassifierParallel extends AbstractClassifier {
	public static boolean chatty= false;
	private RConnection c = null;
	int cores;
	
	@Override
	public void start() {
	
		//if(!StartRserve.isRserveRunning()){
			StartRserve.launchRserve("R");
			c = StartRserve.c;
			
			/*
		}else{
			try {
				c = new RConnection();
			} catch (RserveException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		*/
		try {
			c.voidEval("library(\"missForest\")");
			c.voidEval("try(library(\"randomForest\"))");
			c.voidEval("try(library(\"doSNOW\"))");
			c.voidEval("try(library(\"foreach\"))");
			c.voidEval("library(\"doParallel\")");
			cores = Runtime.getRuntime().availableProcessors();
			c.voidEval("registerDoParallel(cores="+cores+")");
			c.voidEval("load(\"/home/thorsten/randomForestModel.RData\")");
			
			c.voidEval("cl<-makeCluster("+cores+")");
			c.voidEval("registerDoSNOW(cl)");
		} catch (RserveException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	@Override
	public void stop() {
		// TODO Auto-generated method stub
		try {
			c.voidEval("stopCluster(cl)");
			c.shutdown();
			
		} catch (RserveException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	@Override
	public String classify(Trajectory t) {
		ArrayList<Trajectory> tracks = new ArrayList<Trajectory>();
		tracks.add(t);
		String[] result = classify(tracks);
		return result[0];
	}
	
	
	@Override
	public String[] classify(ArrayList<Trajectory> tracks) {
		
		/*
		 * Generate Dataset
		 */
		
		int N = tracks.size();
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
		System.out.println("Calc features");
		long startTime = System.currentTimeMillis();
		//selectedFeatures=c("TYPES","FD","TRAPPED", "EFFICENCY","POWER","SD.DIR","SPLINE.RATIO")
		for(int i = 0; i < tracks.size(); i++){
			Trajectory t = tracks.get(i);
			int timelagForDirectionDeviationLong = t.size()/20; 
			
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
		long estimatedTime = System.currentTimeMillis() - startTime;
		System.out.println("Calc features finished:" + estimatedTime);
		/*
		 * Classify
		 */

		String[] res =null;
	//	RConnection c = StartRserve.c; 
		
		try {
			
			c.assign("elong", elong);
			c.assign("fd",fd);
			c.assign("msdcurvature",msdcurvature);
			c.assign("power", power);
		//	c.assign("sdDir", sdDir);
			c.assign("D.ratio", dratio);
			c.assign("LtStRatio", ltStRatio);
			c.assign("asymmetry1", asym1);
			c.assign("asymmetry2", asym2);
			c.assign("asymmetry3", asym3);
			c.assign("efficiency", efficiency);
			c.assign("kurtosis",kurtosis);
			c.assign("skewness", skewness);
			c.assign("msdratio", msdratio);
			c.assign("straightness", straightness);
			c.assign("trappedness", trappedness);
			c.assign("gaussianity", gaussianity);
			/*
			c.voidEval("data<-data.frame(ELONG=elong,FD=fd,MSD.C=msdcurvature,"
					+ "POWER=power,SPLINE.RATIO=D.ratio,LTST.RATIO=LtStRatio,"
					+ "ASYM1=asymmetry1,MSD.R=msdratio,ASYM2=asymmetry2,ASYM3=asymmetry3,EFFICENCY=efficiency, KURT=kurtosis,"
					+ "SKEW=skewness,STRAIGHTNESS=straightness, "
					+ "TRAPPED=trappedness,GAUSS=gaussianity)");
			*/
			c.voidEval("data<-data.frame(FD=fd, KURT=kurtosis,"
					+ "POWER=power,STRAIGHTNESS=straightness,ASYM3=asymmetry3,"
					+ "MSD.R=msdratio,EFFICENCY=efficiency,LTST.RATIO=LtStRatio,"
					+ "TRAPPED=trappedness,GAUSS=gaussianity)");
			
			//c.voidEval("try(hlp<-missForest(data))");//,parallelize = \"variables\")");
			
			//c.voidEval("data<-hlp$ximp");
		
			
			
			
			c.voidEval("split_testing<-sort(rank(1:nrow(data))%%"+cores+") ");
			c.voidEval("pred<-foreach(i=unique(split_testing),"
					+ ".combine=c,.packages=c(\"randomForest\")) %dopar% {"
							+ "as.numeric(predict(model,newdata=data[split_testing==i,]))}");
			c.voidEval("lvl<-levels(model$y)");
			//c.voidEval("features.predict <- predict(model,data)");
			res = c.eval("lvl[pred]").asStrings();
		} catch (RserveException e) {
			System.out.println("Message: " + e.getMessage());
			e.printStackTrace();
			try {
				c.voidEval("save(data,file=\"/home/thorsten/baddata.Rdata\")");
				System.out.println("Bad data is saved...");
			} catch (RserveException e1) {
				System.out.println("NOT");
				e1.printStackTrace();
			}
			
			
			return null;
		} catch (REngineException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (REXPMismatchException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		return res;
	}
	
	public static void printArr(double[] arr, String title){
		System.out.print(title+":");
		System.out.print(arr[0]);
		for(int i = 1; i < arr.length; i++){
			System.out.print(","+arr[i]);
		}
		System.out.println();
		
	}

	
	
	

}
