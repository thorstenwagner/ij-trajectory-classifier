package de.biomedical_imaging.ij.trajectory_classifier.help;

import ij.IJ;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.junit.experimental.ParallelComputer;
import org.knowm.xchart.Chart;

import com.opencsv.CSVWriter;

import de.biomedical_imaging.ij.trajectory_classifier.AbstractClassifier;
import de.biomedical_imaging.ij.trajectory_classifier.RRFClassifierRenjin;
import de.biomedical_imaging.ij.trajectory_classifier.WeightedWindowedClassificationProcess;
import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.TrajectoryUtil;
import de.biomedical_imaging.traJ.VisualizationUtils;
import de.biomedical_imaging.traJ.features.BoundednessFeature;
import de.biomedical_imaging.traJ.features.MeanSquaredDisplacmentFeature;
import de.biomedical_imaging.traJ.features.TrappedProbabilityFeature;
import de.biomedical_imaging.traJ.simulation.AbstractSimulator;
import de.biomedical_imaging.traJ.simulation.ActiveTransportSimulator;
import de.biomedical_imaging.traJ.simulation.AnomalousDiffusionWMSimulation;
import de.biomedical_imaging.traJ.simulation.CentralRandomNumberGenerator;
import de.biomedical_imaging.traJ.simulation.CombinedSimulator;
import de.biomedical_imaging.traJ.simulation.ConfinedDiffusionSimulator;
import de.biomedical_imaging.traJ.simulation.FreeDiffusionSimulator;
import de.biomedical_imaging.traJ.simulation.SimulationUtil;

public class PerformanceDataGenerator {
	
	public static double diffusioncoefficient = 9.02*Math.pow(10,-2);
	public static String modelpath="";
	
	public static void main(String[] args) {
		PerformanceDataGenerator pg = new PerformanceDataGenerator();
		try {
			modelpath= pg.ExportResource("/randomForestModel.RData");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		generateNormalDiffData();
		generateActiveData();
		generateConfinedData();
		generateSubdiffusionData();
	}
	
	public static void generateActiveData(){
		CentralRandomNumberGenerator.getInstance().setSeed(10);
		int tracklength = 500;
		double dt = 1.0/30;
		double r = 1;
		double[] SNRs = {1,2,3,4,5};
		double angleVelocity = Math.PI/4;
		int N = 50;
		int w =  30;
		AbstractClassifier rrf = new RRFClassifierRenjin(modelpath,dt);
		rrf.start();
		WeightedWindowedClassificationProcess wcp = new WeightedWindowedClassificationProcess();
		double[] par = new double[15];
		double[] cor = new double[15];
		double[] sdcor = new double[15];
		double[] semcor = new double[15];
		int j  =0;
		Mean m = new Mean();
		StandardDeviation sd = new StandardDeviation();
		for (double SNR : SNRs) {
		r = 1;	
		
		while(r <= 15){
		
			double[] val = new double[N];
			for(int i = 0; i < N; i++){
				
				double drift = Math.sqrt(r*4*diffusioncoefficient/((w*2)*dt));
				double sigmaPosNoise = Math.sqrt(diffusioncoefficient*dt+drift*drift*dt*dt)/SNR; 
				//System.out.println("sigma " + sigmaPosNoise);
				AbstractSimulator sim1 = new ActiveTransportSimulator(drift, angleVelocity, dt, 2, tracklength);
				AbstractSimulator sim2 = new FreeDiffusionSimulator(diffusioncoefficient, dt, 2, tracklength);
				AbstractSimulator sim = new CombinedSimulator(sim1, sim2);
				Trajectory t = sim.generateTrajectory();
				t = SimulationUtil.addPositionNoise(t, sigmaPosNoise);
				String[] res = wcp.windowedClassification(t, rrf, w);
				val[i] = getCorrectNess(res, "DIRECTED/ACTIVE");
			}
			double meancorrectness = m.evaluate(val);
			double corrsd = sd.evaluate(val);
			System.out.println("r: " + r+ " Meancorr: " + meancorrectness);
			par[j] = r;
			cor[j] = meancorrectness;
			sdcor[j] = corrsd;
			semcor[j] = corrsd/Math.sqrt(N);
			j++;
			r+=1;
		}
		
		exportCSV("/home/thorsten/perform_active_SNR_"+SNR+"_N_"+N+".csv", "r", par, cor,sdcor,semcor);
		}
		rrf.stop();
	}
	
	public static void generateSubdiffusionData(){
		CentralRandomNumberGenerator.getInstance().setSeed(10);
		int tracklength = 500;
		double dt = 1.0/30;
		double[] alphas = {0.9,0.8,0.7,0.5};
		double SNR = 5;
		int N = 80;
		AbstractClassifier rrf = new RRFClassifierRenjin(modelpath,dt);
		rrf.start();
		WeightedWindowedClassificationProcess wcp = new WeightedWindowedClassificationProcess();
		double[] par = new double[17];
		double[] cor = new double[17];
		double[] sdcor = new double[17];
		double[] semcor = new double[17];
		int j  =0;
		Mean m = new Mean();
		StandardDeviation sd = new StandardDeviation();
		
		for (double alpha : alphas) {
			System.out.println("Start alpha: " + alpha);
			j=0;
			AbstractSimulator sim1 = new AnomalousDiffusionWMSimulation(
					diffusioncoefficient, dt, 2, 2000, alpha);
			for (int w = 15; w < 100; w += 5) {
				System.out.println("Start w: " + w);
				double[] val = new double[N];
				double sigmaPosNoise = Math.sqrt(diffusioncoefficient * dt)
						/ SNR;

				// AbstractSimulator sim1 = new
				// FreeDiffusionSimulator(diffusioncoefficient, 1.0/30, 2,
				// 2000);
				for (int i = 0; i < N; i++) {
					//System.out.println(1.0*i/(N-1)*100);
					Trajectory t = sim1.generateTrajectory();
					t = t.subList(0, tracklength);
					t = SimulationUtil.addPositionNoise(t, sigmaPosNoise);
					String[] res = wcp.windowedClassification(t, rrf, w);
					if (i == 5) {
						System.out.println("res " + res[250]);
					}
					val[i] = getCorrectNess(res, "SUBDIFFUSION");
				}
				double meancorrectness = m.evaluate(val);
				double corrsd = sd.evaluate(val);
				System.out.println("w: " + w*2 + "Alpha: " + alpha + " SNR: "
						+ SNR + " MEAN: " + meancorrectness + " SD: " + corrsd
						/ Math.sqrt(N));
				par[j] = 2*w;
				cor[j] = meancorrectness;
				sdcor[j] = corrsd;
				semcor[j] = corrsd/Math.sqrt(N);
				j++;
			}

			
			exportCSV("/home/thorsten/perform_subdiffusion_SNR_" + SNR + "_N_"
					+ N + "_alpha_" + alpha + "_.csv", "r", par, cor, sdcor,semcor);
		}
		rrf.stop();
	}

	public static void generateConfinedData() {
		CentralRandomNumberGenerator.getInstance().setSeed(10);
		int tracklength = 500;
		double dt = 1.0 / 30;
		double B = 0.5;
		double[] SNRs = {1,2,3,4,5};

		int N = 200;
		int w = 30;
		AbstractClassifier rrf = new RRFClassifierRenjin(modelpath, dt);
		rrf.start();
		WeightedWindowedClassificationProcess wcp = new WeightedWindowedClassificationProcess();
		double[] par = new double[30];
		double[] cor = new double[30];
		double[] sdcor = new double[30];
		double[] semcor = new double[30];
		int j = 0;
		Mean m = new Mean();
		StandardDeviation sd = new StandardDeviation();
		for (double SNR : SNRs) {
		B= 0.5;	
		
		while (B <= 4) {

			double[] val = new double[N];
			double sigmaPosNoise = Math.sqrt(diffusioncoefficient * dt) / SNR;
			double radius = Math.sqrt(BoundednessFeature.a(2 * w)
					* diffusioncoefficient * dt / (4 * B));
			AbstractSimulator sim1 = new ConfinedDiffusionSimulator(
					diffusioncoefficient, dt, radius, 2, tracklength);
			for (int i = 0; i < N; i++) {

				Trajectory t = sim1.generateTrajectory();
				t = SimulationUtil.addPositionNoise(t, sigmaPosNoise);
				String[] res = wcp.windowedClassification(t, rrf, w);
				val[i] = getCorrectNess(res, "CONFINED");
			}
			double meancorrectness = m.evaluate(val);
			double corrsd = sd.evaluate(val);
			System.out.println("B: " + B + " R: " + radius + " MEAN: "
					+ meancorrectness + " SD: " + corrsd);
			par[j] = B;
			cor[j] = meancorrectness;
			sdcor[j] = corrsd;
			semcor[j] = corrsd/Math.sqrt(N);
			j++;
			B += 0.2;
		}
		
		exportCSV("/home/thorsten/perform_confined_SNR_" + SNR + "_N_" + N
				+ "_B_" + B + "_.csv", "r", par, cor, sdcor,semcor);
		}
		rrf.stop();
	}

	public static void generateNormalDiffData() {
		CentralRandomNumberGenerator.getInstance().setSeed(10);
		int tracklength = 500;
		double dt = 1.0/30;
		double[] SNRs = {1,2,3,4,5};
		int N = 200;

		AbstractClassifier rrf = new RRFClassifierRenjin(modelpath, dt);
		rrf.start();
		WeightedWindowedClassificationProcess wcp = new WeightedWindowedClassificationProcess();
		double[] par = new double[16];
		double[] cor = new double[16];
		double[] sdcor = new double[16];
		double[] semcor = new double[16];
		Mean m = new Mean();
		StandardDeviation sd = new StandardDeviation();
		for (double SNR : SNRs) {
			System.out.println("SNR: " + SNR);
			int w = 15;
			int j = 0;
			while (w <= 90) {

				double[] val = new double[N];
				for (int i = 0; i < N; i++) {

					double sigmaPosNoise = Math.sqrt(diffusioncoefficient * dt)
							/ SNR;
					AbstractSimulator sim = new FreeDiffusionSimulator(
							diffusioncoefficient, dt, 2, tracklength);
					Trajectory t = sim.generateTrajectory();
					t = SimulationUtil.addPositionNoise(t, sigmaPosNoise);
					String[] res = wcp.windowedClassification(t, rrf, w);
					val[i] = getCorrectNess(res, "NORM. DIFFUSION");
				}
				double meancorrectness = m.evaluate(val);
				double corrsd = sd.evaluate(val);
				System.out.println("w: " + w + " Meancorr: " + meancorrectness + "SEM: " + corrsd/Math.sqrt(N));
				par[j] = 2 * w;
				cor[j] = meancorrectness;
				sdcor[j] = corrsd;
				semcor[j] = corrsd/Math.sqrt(N);
				j++;
				w += 5;
			}

			exportCSV("/home/thorsten/perform_normal_SNR_" + SNR + "_N_" + N
					+ ".csv", "w", par, cor, sdcor,semcor);

		}
		rrf.stop();
	}
	
	public static void exportCSV(String path, String parameter, double[] par, double[] cor, double[] corsd,double[] semcor){
		String[] nextLine = null;
		
			CSVWriter writer;
			try {
				writer = new CSVWriter(new FileWriter(path, false), '\t');
				nextLine = new String[]{parameter,"correctness","sd","sem"};
				writer.writeNext(nextLine);
				
				for(int i = 0; i < cor.length; i++){
					double parv = par[i];
					double corv = cor[i];
					double sd = corsd[i];
					double sem = semcor[i];
					nextLine = new String[]{""+parv,""+corv,""+sd};
					writer.writeNext(nextLine);
					
				}
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		
	}
	
	public static double getCorrectNess(String[] actual, String expclass){
		int corr = 0;
		for(int i = 0; i < actual.length; i++){
			if(actual[i].equals(expclass)){
				corr++;
			}
		}
		return 1.0*corr/actual.length;
	}
	
	 public String ExportResource(String resourceName) throws Exception {
	        InputStream stream = null;
	        OutputStream resStreamOut = null;
	        String tmpFolder;
	        try {
	            stream = this.getClass().getResourceAsStream(resourceName);//note that each / is a directory down in the "jar tree" been the jar the root of the tree
	            if(stream == null) {
	            	IJ.error("Cannot get resource \"" + resourceName + "\" from Jar file.");
	                throw new Exception("Cannot get resource \"" + resourceName + "\" from Jar file.");
	            }

	            int readBytes;
	            byte[] buffer = new byte[4096];
	            File folderDir = new File(IJ.getDirectory("temp")+"/.trajclassifier");

	            // if the directory does not exist, create it
	            if (!folderDir.exists()) {
	            	folderDir.mkdir();
	            }
	            tmpFolder = folderDir.getPath().replace('\\', '/');
	            resStreamOut = new FileOutputStream(tmpFolder + resourceName);
	            while ((readBytes = stream.read(buffer)) > 0) {
	                resStreamOut.write(buffer, 0, readBytes);
	            }
	        } catch (Exception ex) {
	        	IJ.error(ex.getMessage());
	            throw ex;
	        } finally {
	            stream.close();
	            resStreamOut.close();
	        }

	        return tmpFolder + resourceName;
	    }


}
