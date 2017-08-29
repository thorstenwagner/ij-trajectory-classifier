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
import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Font;
import java.awt.font.FontRenderContext;
import java.awt.geom.AffineTransform;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;
import org.omg.CORBA.OMGVMCID;

import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.AbstractDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.CovarianceDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.RegressionDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.drift.StaticDriftCalculator;
import de.biomedical_imaging.traJ.drift.StaticDriftCorrector;
import de.biomedical_imaging.traJ.features.AbstractTrajectoryFeature;
import de.biomedical_imaging.traJ.features.ActiveTransportParametersFeature;
import de.biomedical_imaging.traJ.features.Asymmetry3Feature;
import de.biomedical_imaging.traJ.features.CenterOfGravityFeature;
import de.biomedical_imaging.traJ.features.ConfinedDiffusionParametersFeature;
import de.biomedical_imaging.traJ.features.EfficiencyFeature;
import de.biomedical_imaging.traJ.features.FractalDimensionFeature;
import de.biomedical_imaging.traJ.features.GaussianityFeauture;
import de.biomedical_imaging.traJ.features.KurtosisFeature;
import de.biomedical_imaging.traJ.features.MSDRatioFeature;
import de.biomedical_imaging.traJ.features.MeanSpeedFeature;
import de.biomedical_imaging.traJ.features.PowerLawFeature;
import de.biomedical_imaging.traJ.features.StraightnessFeature;
import de.biomedical_imaging.traJ.features.TrappedProbabilityFeature;
import ij.IJ;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.io.OpenDialog;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;

public class TraJClassifier_ implements PlugIn {

	private double timelag;
	private double minTrackLength;
	private int windowSizeClassification;
	private int minSegmentLength;
	private double pixelsize;
	private int resampleRate;
	private boolean showID;
	private boolean showOverviewClasses;
	private boolean removeGlobalDrift;
	private boolean useReducedModelConfinedMotion;
	private int ommittedTrajectories;
	private ArrayList<Subtrajectory> classifiedSegments;
	private ArrayList<Trajectory> tracksToClassify;
	//private ArrayList<Trajectory> tracks 
	private static TraJClassifier_ instance;
	ArrayList<Trajectory> parentTrajectories;
	String documentPath;
	public TraJClassifier_() {
		minTrackLength=160;
		windowSizeClassification=60;
		minSegmentLength=30;
		pixelsize=0.166;
		timelag=1.0/30;
		resampleRate=1;
		showID = true;
		ommittedTrajectories=0;
		instance = this;
	}
	
	public static TraJClassifier_ getInstance(){
		if(instance == null){
			instance = new TraJClassifier_();
		}
		return instance;
	}
	
	public void run(String arg) {
		
		/*
		 * Export model!
		 */
		String modelpath="";
		try {
			modelpath= ExportResource("/randomForestModel.RData");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		/*
		 * Import Data
		 */
		if(!arg.contains("DEBUG")){
			OpenDialog open = new OpenDialog("Choose the TrackMate xml file");
			documentPath = open.getPath();
			if(documentPath==null){
				return;
			}
			TrackMateImporter tMateImport = new TrackMateImporter();
			tracksToClassify = tMateImport.importTrackMateXML(documentPath);
		}
		int maxNumberOfPositions = 0;
		for(int i = 0; i < tracksToClassify.size(); i++){
			if(tracksToClassify.get(i).size()> maxNumberOfPositions){
				maxNumberOfPositions = tracksToClassify.get(i).size();
			}
		}
		boolean visualize = true;
		
		if(WindowManager.getCurrentImage()==null || (IJ.getImage().getNFrames()<maxNumberOfPositions && IJ.getImage().getNSlices()<maxNumberOfPositions)){
		
				IJ.showMessage("Your image does not contain enough frames for visualization."
						+ "Therefore visualization will be deactivated. \n "
						+ "For visualization please open the corresponding image stack to the trackmate file");
				visualize = false;
			
		}

		/*
		 * GUI
		 */
		
		String classifierVersion = getClass().getPackage().getImplementationVersion();
		
		if(!arg.contains("NOGUI")){
			//Load previous settings
			
			minTrackLength = Prefs.get("trajclass.minTrackLength", 160);
			windowSizeClassification = (int) Prefs.get("trajclass.windowSize", 60);
			minSegmentLength = (int) Prefs.get("trajclass.minSegmentLength", 60);
			resampleRate = (int) Prefs.get("trajclass.ResampleRate", 1);
			pixelsize = Prefs.get("trajclass.pixelsize", 0.166);
			timelag = 1.0/Prefs.get("trajclass.framerate",30);
			useReducedModelConfinedMotion = Prefs.get("trajclass.useReducedModelConfinedMotion", false);
			showID = Prefs.get("trajclass.showID", true);
			showOverviewClasses = Prefs.get("trajclass.showOverviewClasses", true);
			removeGlobalDrift = Prefs.get("trajclass.removeGlobalDrift", false);
			
			//Show GUI
			GenericDialog gd = new GenericDialog("TraJectory Classification ("+classifierVersion+")");
		
			gd.addSlider("Min. tracklength", 10, 1000, minTrackLength);
			gd.addSlider("Windowsize (positions)", 10, 1000, windowSizeClassification);
			gd.addSlider("Min. segment length",10,1000,minSegmentLength);
			gd.addNumericField("Resample rate*", resampleRate, 0);
			gd.addNumericField("Pixelsize (µm)**", pixelsize, 4);
			gd.addNumericField("Framerate (FPS)", 1/timelag, 0);
			gd.addCheckbox("Use reduced model confined motion", useReducedModelConfinedMotion);
			gd.addCheckbox("Show_IDs", showID);
			gd.addCheckbox("Show_overview classes", showOverviewClasses);
			gd.addCheckbox("Remove_global_drift", removeGlobalDrift);
			gd.addMessage("* The ratio of window size / resample rate have to be at least 30.");
			gd.addMessage("** Set to zero if the imported data is already correctly scaled.");
			gd.addHelp("http://imagej.net/TraJClassifier");
			gd.addDialogListener(new DialogListener() {
				
				@Override
				public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
					minTrackLength = gd.getNextNumber();
					windowSizeClassification = (int) (gd.getNextNumber());
					minSegmentLength = (int)gd.getNextNumber();
					resampleRate = (int)gd.getNextNumber();
					pixelsize = gd.getNextNumber();
					timelag = 1/gd.getNextNumber();
					useReducedModelConfinedMotion = gd.getNextBoolean();
					showID = gd.getNextBoolean();
					showOverviewClasses = gd.getNextBoolean();
					removeGlobalDrift = gd.getNextBoolean();
					boolean valid = (resampleRate == 1 || windowSizeClassification/resampleRate>=30) && 
							(minTrackLength >= windowSizeClassification) &&
							(minTrackLength>=10) && (minSegmentLength>=10) && (windowSizeClassification>=10);
					
					return valid;
				}
			});
			gd.showDialog();
			if(gd.wasCanceled()){
				return;
			}
			minTrackLength = gd.getNextNumber();
			windowSizeClassification = (int) (gd.getNextNumber()/2);
			minSegmentLength = (int)gd.getNextNumber();
			resampleRate = (int)gd.getNextNumber();
			pixelsize = gd.getNextNumber();
			timelag = 1/gd.getNextNumber();
			useReducedModelConfinedMotion = gd.getNextBoolean();
			showID = gd.getNextBoolean();
			showOverviewClasses = gd.getNextBoolean();
			removeGlobalDrift = gd.getNextBoolean();
			// Save settings
			Prefs.set("trajclass.minTrackLength", minTrackLength);
			Prefs.set("trajclass.windowSize", windowSizeClassification*2);
			Prefs.set("trajclass.minSegmentLength", minSegmentLength);
			Prefs.set("trajclass.ResampleRate", resampleRate);
			Prefs.set("trajclass.pixelsize", pixelsize);
			Prefs.set("trajclass.framerate", 1/timelag);
			Prefs.set("trajclass.showID", showID);
			Prefs.set("trajclass.showOverviewClasses", showOverviewClasses);
			Prefs.set("trajclass.removeGlobalDrift", removeGlobalDrift);
			Prefs.set("trajclass.useReducedModelConfinedMotion", useReducedModelConfinedMotion);
		}

		
		/*
		 * Scaling
		 */
		if(pixelsize>0.000001){
			for(int i = 0; i < tracksToClassify.size(); i++){
				tracksToClassify.get(i).scale(pixelsize);
			}
		}
		
		/*
		 * Trajectories which to often change its position are not suitable for the classifier. At least for 50% percent
		 * of the trajectories, the position have to change.
		 */
		
		for(int i = 0; i < tracksToClassify.size(); i++){
			Trajectory t = tracksToClassify.get(i);
			int changesCounter=0;
			for(int j = 1; j < t.size(); j++){
				if(t.get(j).distance(t.get(j-1)) > Math.pow(10,-12)){
					changesCounter++;
				}
			}
			if(1.0*changesCounter/t.size() < 0.5){
				tracksToClassify.remove(i);
				ommittedTrajectories++;
				i--;
			}
		}
		
		
	
		/*
		 * Classification, Segmentation & Visualization
		 */

		HashMap<String, Color> mapTypeToColor = new HashMap<String, Color>();
		mapTypeToColor.put("DIRECTED/ACTIVE", Color.MAGENTA);
		mapTypeToColor.put("NORM. DIFFUSION", Color.RED);
		mapTypeToColor.put("CONFINED", Color.YELLOW);
		mapTypeToColor.put("SUBDIFFUSION", Color.GREEN);

		/*
		 *  Remove tracks which are too short
		 */
		parentTrajectories = new ArrayList<Trajectory>();
		for (Trajectory track : tracksToClassify) {
			if(track.size() > minTrackLength){
				parentTrajectories.add(track);
			}
		}

		/*
		 * Remove drift
		 */
		StaticDriftCalculator<Trajectory> dcalc = new StaticDriftCalculator<Trajectory>();
		if(removeGlobalDrift){
			double[] drft = dcalc.calculateDrift(parentTrajectories);
			StaticDriftCorrector dcorr = new StaticDriftCorrector(drft);
			parentTrajectories = dcorr.removeDrift(parentTrajectories);
		}
 		
		/*
		 * Classification and segmentation
		 */
		
		classifiedSegments = classifyAndSegment(parentTrajectories, modelpath, windowSizeClassification, minSegmentLength, 10, resampleRate);


		/*
		 * Visualization 
		 */
		if (visualize) {
			// Trajectories
			Overlay ov = new Overlay();
			for (int i = 0; i < classifiedSegments.size(); i++) {
				Subtrajectory tr = classifiedSegments.get(i);

				ArrayList<Roi> prois = null;
				if (pixelsize > 0.000001) {
					prois = VisualizationUtils
							.generateVisualizationRoisFromTrack(tr,
									mapTypeToColor.get(tr.getType()), showID,
									pixelsize);

				} else {
					prois = VisualizationUtils
							.generateVisualizationRoisFromTrack(tr,
									mapTypeToColor.get(tr.getType()), showID,IJ.getImage().getCalibration().pixelWidth);
				}
				for (Roi r : prois) {
					ov.add(r);
				}
			}

			// Classes
			if (showOverviewClasses) {
				Set<String> classes = mapTypeToColor.keySet();

				Iterator<String> it = classes.iterator();
				int y = 5;
				
				float fsize = 20;
				AffineTransform affinetransform = new AffineTransform();     
				FontRenderContext frc = new FontRenderContext(affinetransform,true,true); 
				int width = (int)IJ.getImage().getProcessor().getFont().getStringBounds("Norm. Diffusion", frc).getWidth();
				Font f = IJ.getImage().getProcessor().getFont();
				while(1.0*width/IJ.getImage().getWidth() > 0.08){
					fsize--;
					f = f.deriveFont(fsize);
					width = (int)f.getStringBounds("Norm. Diffusion", frc).getWidth();
				}
				
				TextRoi.setFont("TimesRoman", (int)fsize, Font.PLAIN);
				while (it.hasNext()) {
					String type = it.next();
					TextRoi troi = new TextRoi(5, y, type);
					troi.setFillColor(Color.DARK_GRAY);
					troi.setStrokeColor(mapTypeToColor.get(type));
					ov.add(troi);
					y = y + 20;

				}
			}
			
			IJ.getImage().setOverlay(ov);
			IJ.getImage().updateAndRepaintWindow();
		}
		

		/*
		 * Fill results table
		 */
		
		
		HashMap<String, TraJResultsTable> rtables = new HashMap<String, TraJResultsTable>();
		rtables.put("DIRECTED/ACTIVE", new TraJResultsTable());
		rtables.put("NORM. DIFFUSION", new TraJResultsTable());
		rtables.put("SUBDIFFUSION", new TraJResultsTable());
		rtables.put("CONFINED", new TraJResultsTable());

		double sumConf = 0;
		for (int i = 0; i < classifiedSegments.size(); i++) {
				IJ.showProgress(i, classifiedSegments.size());
				Subtrajectory t = classifiedSegments.get(i);
				TraJResultsTable rt = rtables.get(t.getType());
				if(rt==null){
					IJ.log("Type: " + t.getType());
					ExportImportTools eit = new ExportImportTools();
					ArrayList<Trajectory> hlp = new ArrayList<Trajectory>();
					eit.exportTrajectoryDataAsCSV(hlp, "/home/thorsten/bad.csv");
					IJ.log(t.toString());
					
				}
				rt.incrementCounter();
				rt.addValue("PARENT-ID", t.getParent().getID());
				rt.addValue("ID", t.getID());
				rt.addValue("LENGTH", t.size());
				rt.addValue("START", t.getRelativeStartTimepoint());
				rt.addValue("END", t.getRelativeStartTimepoint()+t.size()-1);
				rt.addValue("CLASS", t.getType());

				AbstractTrajectoryFeature dcEstim=null;
				double dc =0;
				DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.ENGLISH);
				NumberFormat formatter = new DecimalFormat("0.###E0",otherSymbols);
				double[] res;
				double goodness=0;
				double alphaGoodness=0;
				switch (t.getType()) {
				case "DIRECTED/ACTIVE":
					
					ActiveTransportParametersFeature apf = new ActiveTransportParametersFeature(t, timelag);
					
					//dcEstim = new RegressionDiffusionCoefficientEstimator(t,1/timelag,1,t.size()-1);
					//dc = dcEstim.evaluate()[0];
					res = apf.evaluate();
					rt.addValue("(FIT) D", formatter.format(res[0]));
					rt.addValue("(FIT) Velocity", res[1]);
					goodness = res[2];
					break;
				case "NORM. DIFFUSION":
					
					dcEstim = new RegressionDiffusionCoefficientEstimator(t, 1/timelag, 1, t.size()/3);
					res = dcEstim.evaluate();
					dc = res[0];
					rt.addValue("(FIT) D", formatter.format(dc));
					goodness = res[3];
					
					break;
				case "CONFINED":
					AbstractDiffusionCoefficientEstimator dcEst = new RegressionDiffusionCoefficientEstimator(t,1/timelag,1,3);
					ConfinedDiffusionParametersFeature confp = new ConfinedDiffusionParametersFeature(t,timelag,useReducedModelConfinedMotion,dcEst);
					double[] p = confp.evaluate();
					dc = p[1];
					if(useReducedModelConfinedMotion){
						rt.addValue("(FIT) CONF. RADIUS", Math.sqrt(p[0]));
						rt.addValue("(FIT) D", formatter.format(p[1]));
						goodness = p[2];
					}else{
						rt.addValue("(FIT) CONF. RADIUS", Math.sqrt(p[0]));
						rt.addValue("(FIT) A [CONF. SHAPE]", p[2]);
						rt.addValue("(FIT) B (CONF SHAPE)", p[3]);
						rt.addValue("(FIT) D", formatter.format(p[1]));
						goodness = p[4];
					}
					
	
					break;
				case "SUBDIFFUSION":
					PowerLawFeature pwf = new PowerLawFeature(t, 1/timelag,1, t.size()/3);
					res = pwf.evaluate();
					dc = res[1];
					
					rt.addValue("(FIT) D", formatter.format(dc));
		
					goodness = res[2];
					break;
				default:
					break;
				}
				
				AbstractTrajectoryFeature pwf = new PowerLawFeature(t, 1/timelag,1, t.size()/3);
				res = pwf.evaluate();
				double alpha = res[0];
				alphaGoodness = res[2];
				
				
				AbstractTrajectoryFeature f = new CenterOfGravityFeature(t);
				double cog_x = f.evaluate()[0];
				double cog_y = f.evaluate()[1];
				rt.addValue("X (COG)", cog_x);
				rt.addValue("Y (COG)", cog_y);
				
				if(t.getType().equals("NONE")==false){
					FractalDimensionFeature fdf = new FractalDimensionFeature(t);
					double v = fdf.evaluate()[0];;
					rt.addValue("FRACT. DIM.", v);
					
					TrappedProbabilityFeature trapped = new TrappedProbabilityFeature(t);
					v = trapped.evaluate()[0];
					rt.addValue("TRAPPEDNESS", v);
					
					EfficiencyFeature eff = new EfficiencyFeature(t);
					v=eff.evaluate()[0];
					rt.addValue("EFFICENCY", v);
					
					StraightnessFeature straight = new StraightnessFeature(t);
					v= straight.evaluate()[0];
					rt.addValue("STRAIGHTNESS", v);
					
					MeanSpeedFeature msfeature = new MeanSpeedFeature(t, timelag);
					v = msfeature.evaluate()[1];
					rt.addValue("SPEED", v);
					
					KurtosisFeature kurt = new KurtosisFeature(t);
					v = kurt.evaluate()[0];
					rt.addValue("KURTOSIS", v);
					
				//	AbstractTrajectoryFeature pwf = new PowerLawFeature(t, 1, t.size()/3);
				//	res = pwf.evaluate();
				//	v = res[0];
					rt.addValue("(FIT) ALPHA", alpha);
					
					GaussianityFeauture gauss = new GaussianityFeauture(t, 1);
					v = gauss.evaluate()[0];
					rt.addValue("GAUSSIANITY", v);
					
					Asymmetry3Feature asym3 = new Asymmetry3Feature(t);
					v = asym3.evaluate()[0];
					rt.addValue("Asymmetry", v);
					
					MSDRatioFeature msdratio = new MSDRatioFeature(t, 1, 5);
					v = msdratio.evaluate()[0];
					rt.addValue("MSDRatio", v);
					
					CovarianceDiffusionCoefficientEstimator cest = new CovarianceDiffusionCoefficientEstimator(t, 1/timelag);
					res = cest.evaluate();
					rt.addValue("Loc. noise_sigma", (res[1]+res[2])/2);
					rt.addValue("Fit Goodness", goodness);
					rt.addValue("Alpha Fit Goodness", alphaGoodness);
					double conf = t.getConfidence();
					sumConf+=conf;
					rt.addValue("Confidence", conf);
					
				}

		}
		
		/*
		 * Fill parents results table
		 */
		Iterator<String> rtIt = rtables.keySet().iterator();
		
		ResultsTable parents = new TraJResultsTable(true);
		
		for(int i = 0; i < parentTrajectories.size(); i++){
			parents.incrementCounter();
			Trajectory t = parentTrajectories.get(i);
			parents.addValue("ID", t.getID());
			parents.addValue("LENGTH", t.size());
			parents.addValue("START", t.getRelativeStartTimepoint());
			parents.addValue("END", t.getRelativeStartTimepoint()+t.size()-1);
			int subPosCount =0;
			int subSegCount =0;
			int normPosCount =0;
			int normSegCount = 0;
			int directedPosCount =0;
			int directSegCount = 0;
			int confPosCount=0;
			int confSegCount=0;

			ArrayList<Subtrajectory> sameParent = Subtrajectory.getTracksWithSameParant(classifiedSegments, t.getID());
			for (Subtrajectory sub : sameParent) {
				switch (sub.getType()) {
				case "DIRECTED/ACTIVE":
					directedPosCount+=sub.size();
					directSegCount++;
					break;
				case "NORM. DIFFUSION":
					normPosCount+=sub.size();
					normSegCount++;
					break;
				case "CONFINED":
					confPosCount+=sub.size();
					confSegCount++;
					break;
				case "SUBDIFFUSION":
					subPosCount+=sub.size();
					subSegCount++;
					break;
				default:
					break;
				}
			}
			parents.addValue("#SEG_NORM", normSegCount);
			parents.addValue("#POS_NORM", normPosCount);
			parents.addValue("#SEG_SUB", subSegCount);
			parents.addValue("#POS_SUB", subPosCount);
			parents.addValue("#SEG_CONF", confSegCount);
			parents.addValue("#POS_CONF", confPosCount);
			parents.addValue("#SEG_DIRECTED", directSegCount);
			parents.addValue("#POS_DIRECTED", directedPosCount);
		}
		
		
		String trajVersion = Trajectory.class.getPackage().getImplementationVersion();
		double[] drift = dcalc.calculateDrift(parentTrajectories);
		ResultsTable overall = new ResultsTable();
		overall.incrementCounter();
		overall.addValue("Mean confindence", sumConf/classifiedSegments.size());
		overall.addValue("Drift x", drift[0]);
		overall.addValue("Drift y", drift[1]);
		overall.addValue("Omitted segments", ommittedTrajectories);
		overall.addValue("Min. track length", minTrackLength);
		overall.addValue("Window size", windowSizeClassification*2);
		overall.addValue("Min. segment length", minSegmentLength);
		overall.addValue("Resamplerate", resampleRate);
		overall.addValue("Pixelsize", pixelsize);
		overall.addValue("Framerate", 1/timelag);
		overall.addValue("Reduced conf. model", Boolean.toString(useReducedModelConfinedMotion));
		overall.addValue("Remove global drift", Boolean.toString(removeGlobalDrift));
		overall.addValue("TraJclassifier version", classifierVersion);
		overall.addValue("TraJ library version", trajVersion);
		overall.addValue("TrackMate file:", documentPath);
		
		overall.show("Settings & Miscellaneous");
		
		
		// show tables
		parents.show("Parents");
		rtIt = rtables.keySet().iterator();
		while(rtIt.hasNext()){
			String rt = rtIt.next();
			rtables.get(rt).show(rt + " trajectories");
		}
		
	
	}
	
	public ArrayList<Subtrajectory> classifyAndSegment(Trajectory trackToClassify, String modelpath, int windowSizeClassification, int minSegmentLength, int modeFilterLength, int resampleRate){
		ArrayList<Trajectory> help = new ArrayList<Trajectory>();
		help.add(trackToClassify);
		return classifyAndSegment(help, modelpath, windowSizeClassification, minSegmentLength, modeFilterLength, resampleRate);
	}
	
	public ArrayList<Subtrajectory> classifyAndSegment(ArrayList<Trajectory> tracksToClassify, String modelpath, int windowSizeClassification, int minSegmentLength, int modeFilterLength, int resampleRate){
		ArrayList<Subtrajectory> classified = new ArrayList<Subtrajectory>();
		int j = 0;
		RRFClassifierRenjin rrf = new RRFClassifierRenjin(modelpath,resampleRate*timelag);
		rrf.start();
		WeightedWindowedClassificationProcess wcp = new WeightedWindowedClassificationProcess();
		int subidcounter = 1;
		for (Trajectory track : tracksToClassify) {
			j++;
			IJ.showProgress(j, tracksToClassify.size());
			Trajectory mTrack = track;

			String[] classes = wcp.windowedClassification(mTrack, rrf, windowSizeClassification,resampleRate);
			
				double[] classConfidence = wcp.getPositionConfidence();
				//Moving mode
				classes = movingMode(classes, modeFilterLength);
				double sumConf = 0;
				int Nconf = 0;
				Subtrajectory tr = new Subtrajectory(track,2);

				tr.setID(subidcounter);
				subidcounter++;
				tr.add(track.get(0).x, track.get(0).y, 0);
				sumConf += classConfidence[0];
				Nconf++;
				String prevCls = classes[0];
				int start = track.getRelativeStartTimepoint();
				tr.setRelativStartTimepoint(start);
				tr.setType(prevCls);

				for(int i = 1; i < classes.length; i++){
					if(prevCls == classes[i]){
						tr.add(track.get(i).x, track.get(i).y,0);
						sumConf += classConfidence[i];
						Nconf++;
					}else{;
					tr.setConfidence(sumConf/Nconf);
					classified.add(tr);
					tr = new Subtrajectory(track,2);
					tr.setID(subidcounter);
					subidcounter++;
					tr.setRelativStartTimepoint(start+i);
					tr.add(track.get(i).x, track.get(i).y,0);
					sumConf = classConfidence[i];
					Nconf = 1;
					prevCls = classes[i];
					tr.setType(prevCls);
					}
				}
				tr.setConfidence(sumConf/Nconf);
				classified.add(tr);
				sumConf = 0;
				Nconf = 0;
			
			
		}
		rrf.stop();
		
		/*
		 * FILTER
		 */
		
		
		
		//Remove segments smaller then the minimum segment length
		for(int i = 0; i < classified.size(); i++){
			if(classified.get(i).size()<minSegmentLength){
				classified.remove(i);
				i--;
			}
		}
		return classified;
	}
	
	public double getTimelag(){
		return timelag;
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
			medTypes[i] = types.get(i);
		}
		for(int i = types.size()- n; i < types.size(); i++){
			medTypes[i] = types.get(i);
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
	
	public ArrayList<Subtrajectory> getClassifiedTrajectories(){
		return classifiedSegments;
	}
	
	public ArrayList<Trajectory> getParentTrajectories(){
		return parentTrajectories;
	}
	
	
	 /**
     * Export a resource embedded into a Jar file to the local file path.
     *
     * @param resourceName ie.: "/SmartLibrary.dll"
     * @return The path to the exported resource
     * @throws Exception
     */
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

    public void setTracksToClassify(ArrayList<Trajectory> t) {
    	tracksToClassify = t;
    }
    

	public double getMinTrackLength() {
		return minTrackLength;
	}

	public void setMinTrackLength(double minTrackLength) {
		this.minTrackLength = minTrackLength;
	}

	public double getPixelsize() {
		return pixelsize;
	}

	public void setPixelsize(double pixelsize) {
		this.pixelsize = pixelsize;
	}

	public boolean isShowID() {
		return showID;
	}

	public void setShowID(boolean showID) {
		this.showID = showID;
	}

	public int getWindowSizeClassification() {
		return windowSizeClassification;
	}
	
	public boolean isUseReducedModelConfinedMotion(){
		return useReducedModelConfinedMotion;
	}

	public void setTimelag(double timelag) {
		this.timelag = timelag;
	}

	
	
	public void setWindowSizeClassification(int windowSizeClassification){
		this.windowSizeClassification = windowSizeClassification;
	}




}
