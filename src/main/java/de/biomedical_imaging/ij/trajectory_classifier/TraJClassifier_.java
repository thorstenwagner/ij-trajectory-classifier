package de.biomedical_imaging.ij.trajectory_classifier;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;

import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.AbstractDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.CovarianceDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.features.AbstractTrajectoryFeature;
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
import de.biomedical_imaging.traJ.features.ShortTimeLongTimeDiffusioncoefficentRatio;
import de.biomedical_imaging.traJ.features.StraightnessFeature;
import de.biomedical_imaging.traJ.features.TrappedProbabilityFeature;
import de.biomedical_imaging.traj.math.ConfinedDiffusionMSDCurveFit.FitMethod;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.io.OpenDialog;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.process.FloatPolygon;

public class TraJClassifier_ implements PlugIn {

	private boolean useR = false;
	private double timelag;
	
	private ArrayList<Subtrajectory> classifiedTrajectories;
	private ArrayList<Trajectory> tracksToClassify;
	//private ArrayList<Trajectory> tracks 
	private static TraJClassifier_ instance;
	
	public TraJClassifier_() {
		instance = this;
	}
	
	public static TraJClassifier_ getInstance(){
		return instance;
	}
	
	public void run(String arg) {
		
		/*
		 * Export model!
		 */
		String modelpath="";
		try {
			modelpath= ExportResource("/randomForestModel.RData");
			System.out.println("P: " + modelpath);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		/*
		 * GUI
		 */
		GenericDialog gd = new GenericDialog("Parameters Classification");
	
		gd.addSlider("Min. tracklength", 1, 1000, 160);
		gd.addSlider("Min. segment size", 0, 1000, 60);
		gd.addSlider("Windowsize", 1, 200, 60);
		gd.addSlider("Filtersize", 0, 200, 0);
		gd.addNumericField("Min. diffusion coeffcient (µm^2 / s)", 0, 0);
		gd.addNumericField("Pixelsize (µm)*", 0.166, 4);
		gd.addNumericField("Framerate", 30, 0);
		gd.addCheckbox("Use parallized R**", false);
		gd.addMessage("* Set to zero if the imported data is already correctly scaled.");
		gd.addMessage("** R has to be installed and correctly configured.");
		gd.addHelp("www.MY HELP PAGE.de");
		gd.showDialog();
		double minlength = gd.getNextNumber();
		double minseglength = gd.getNextNumber();
		int windowSizeClassification = (int) (gd.getNextNumber()/2);
		int modeSize = (int) (gd.getNextNumber()/2);
		double minDiffusionCoefficient = gd.getNextNumber();
		double pixelsize = gd.getNextNumber();
		timelag = 1/gd.getNextNumber();
		useR = gd.getNextBoolean();

		/*
		 * Import Data
		 */
		if(!arg.contains("DEBUG")){
			OpenDialog open = new OpenDialog("Choose the TrackMate xml file");
			String filepath = open.getPath();
			TrackMateImporter tMateImport = new TrackMateImporter();
			tracksToClassify = tMateImport.importTrackMateXML(filepath);
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
		 * Classification, Segmentation & Visualization
		 */
		AbstractClassifier rrf = null;
		if(useR){
			rrf = new RRFClassifierParallel(modelpath);
		}else{
			rrf = new RRFClassifierRenjin(modelpath);
		}
		
		rrf.start();
		WeightedWindowedClassificationProcess wcp = new WeightedWindowedClassificationProcess();
		
		
		int j = 0;

		HashMap<String, Color> mapTypeToColor = new HashMap<String, Color>();
		mapTypeToColor.put("DIRECTED/ACTIVE", Color.MAGENTA);
		mapTypeToColor.put("NORM. DIFFUSION", Color.RED);
		mapTypeToColor.put("CONFINED", Color.YELLOW);
		mapTypeToColor.put("SUBDIFFUSION", Color.GREEN);
		mapTypeToColor.put("STALLED", Color.ORANGE);
		mapTypeToColor.put("NONE", Color.LIGHT_GRAY);

		/*
		 *  Remove tracks which are too short
		 */
		ArrayList<Trajectory> minLengthTracks = new ArrayList<Trajectory>();
		for (Trajectory track : tracksToClassify) {
			if(track.size() > minlength){
				minLengthTracks.add(track);
			}
		}
		
		/*
		 * Classification and segmentation
		 */
		classifiedTrajectories = new ArrayList<Subtrajectory>();
		for (Trajectory track : minLengthTracks) {
			j++;
			IJ.showProgress(j, minLengthTracks.size());

			String[] classes = wcp.windowedClassification(track, rrf, windowSizeClassification);
			
			if(modeSize>0){
				classes = movingMode(classes, modeSize);
			}

			Subtrajectory tr = new Subtrajectory(track,2);
			tr.add(track.get(0).x, track.get(0).y, 0);
			
			String prevCls = classes[0];
			int start = track.getRelativeStartTimepoint();
			tr.setRelativStartTimepoint(start);
			tr.setType(prevCls);
			
			for(int i = 1; i < classes.length; i++){
				//System.out.println("i " + i + " class: "+classes[i]);
				if(prevCls == classes[i]){
					tr.add(track.get(i).x, track.get(i).y,0);
				}else{;
					
					classifiedTrajectories.add(tr);
					tr = new Subtrajectory(track,2);
					tr.setRelativStartTimepoint(start+i);
					tr.add(track.get(i).x, track.get(i).y,0);
					prevCls = classes[i];
					tr.setType(prevCls);
				}
			}
			classifiedTrajectories.add(tr);
		}
		rrf.stop();
		
		/*
		 * FILTER
		 */
		
		//Remove tracks shorter than min. length
		for(int i = 0; i < classifiedTrajectories.size(); i++){
			if(classifiedTrajectories.get(i).size()<minseglength){
				classifiedTrajectories.remove(i);
				i--;
			}
		}
	
		//Segments smaller than the window size seems to be very error prone. Set them to "None".
		for (Trajectory trajectory : classifiedTrajectories) {
			if(trajectory.size()<windowSizeClassification*2){
				trajectory.setType("NONE");
			}
		}
		
		//Trajectories with a negative short time diffusion coefficient are considered as stalled particles
		for (Trajectory trajectory : classifiedTrajectories) {
			CovarianceDiffusionCoefficientEstimator regest = new CovarianceDiffusionCoefficientEstimator(trajectory, 1/timelag);
			if(regest.evaluate()[0]<minDiffusionCoefficient){
				trajectory.setType("STALLED");
			}
		}
		
		/*
		 * Visualization
		 */
		
		Overlay ov = new Overlay(); 
		for(int i = 0; i < classifiedTrajectories.size(); i++){
			Trajectory tr =  classifiedTrajectories.get(i);
			ArrayList<Roi> prois = null;
			if(pixelsize>0.000001){
				prois = VisualizationUtils.generateVisualizationRoisFromTrack(tr, mapTypeToColor.get(tr.getType()),pixelsize);
				
			}else{
				prois = VisualizationUtils.generateVisualizationRoisFromTrack(tr, mapTypeToColor.get(tr.getType()));
			}
			for (Roi r : prois) {
				ov.add(r);
			}
		}
		
		Set<String> classes = mapTypeToColor.keySet();

		Iterator<String> it = classes.iterator();
		int y = 5;
		TextRoi.setFont("TimesRoman", 12, Font.PLAIN);
		while(it.hasNext()){
			String type = it.next();
			TextRoi troi = new TextRoi(5, y, type);

			troi.setStrokeColor(mapTypeToColor.get(type));
			ov.add(troi);
			y = y + 10;

		}

		IJ.getImage().setOverlay(ov);
		IJ.getImage().updateAndRepaintWindow();
		
		/*
		 * Export classified trajectories
		 */
		IJ.log("Export");
		ExportImportTools iotools = new ExportImportTools();
		iotools.exportTrajectoryDataAsCSV(classifiedTrajectories, "/home/thorsten/myclassifiedtracks.csv");
		
		/*
		 * Fill results table
		 */
		ResultsTable parents = new ResultsTable();
		for(int i = 0; i < minLengthTracks.size(); i++){
			parents.incrementCounter();
			Trajectory t = minLengthTracks.get(i);
			parents.addValue("ID", t.getID());
			parents.addValue("LENGTH", t.size());
			parents.addValue("START", t.getRelativeStartTimepoint());
			parents.addValue("END", t.getRelativeStartTimepoint()+t.size()-1);
		}
		
		HashMap<String, TraJResultsTable> rtables = new HashMap<>();
		rtables.put("DIRECTED/ACTIVE", new TraJResultsTable());
		rtables.put("NORM. DIFFUSION", new TraJResultsTable());
		rtables.put("SUBDIFFUSION", new TraJResultsTable());
		rtables.put("CONFINED", new TraJResultsTable());
		rtables.put("STALLED", new TraJResultsTable());
		rtables.put("NONE", new TraJResultsTable());
		
		
		IJ.log("Fill results table");
		for (int i = 0; i < classifiedTrajectories.size(); i++) {
				IJ.showProgress(i, classifiedTrajectories.size());
				Subtrajectory t = classifiedTrajectories.get(i);
				System.out.println("Type: " + t.getType());
				TraJResultsTable rt = rtables.get(t.getType());
				if(rt == null){
				System.out.println("null!");
				}
				rt.incrementCounter();
				rt.addValue("ID", t.getID());
				IJ.log("ID: " + t.getID());
				rt.addValue("PARENT-ID", t.getParent().getID());
				rt.addValue("LENGTH", t.size());
				rt.addValue("START", t.getRelativeStartTimepoint());
				rt.addValue("END", t.getRelativeStartTimepoint()+t.size()-1);
				rt.addValue("CLASS", t.getType());
				//CovarianceDiffusionCoefficientEstimator cov = new CovarianceDiffusionCoefficientEstimator(t, fps)
				//rt.addValue("D (SHORT)", value);
				CovarianceDiffusionCoefficientEstimator covest=null;
				switch (t.getType()) {
				case "DIRECTED/ACTIVE":
					covest = new CovarianceDiffusionCoefficientEstimator(t, 1/timelag);
					rt.addValue("D", String.format("%6.3e",covest.evaluate()[0]));
					break;
				case "NORM. DIFFUSION":
					covest = new CovarianceDiffusionCoefficientEstimator(t, 1/timelag);
					rt.addValue("D", String.format("%6.3e",covest.evaluate()[0]));
					break;
				case "CONFINED":
					AbstractDiffusionCoefficientEstimator dcEst = new CovarianceDiffusionCoefficientEstimator(t,1/timelag);
					ConfinedDiffusionParametersFeature confp = new ConfinedDiffusionParametersFeature(t,timelag,dcEst,FitMethod.SIMPLEX);
			
					double[] p = confp.evaluate();
					rt.addValue("CONF. SIZE", p[0]);
					covest = new CovarianceDiffusionCoefficientEstimator(t, 1/timelag);
					rt.addValue("D", String.format("%6.3e",covest.evaluate()[0]));
					break;
					
				case "SUBDIFFUSION":
					PowerLawFeature pwf = new PowerLawFeature(t, 1, t.size()/3);
					double res[] = pwf.evaluate();
					rt.addValue("D", String.format("%6.3e",res[1]));
					break;
				case "STALLED":
					break;
				case "NONE":
					break;
				default:
					break;
				}
				AbstractTrajectoryFeature f = new CenterOfGravityFeature(t);
				double cog_x = f.evaluate()[0];
				double cog_y = f.evaluate()[1];
				rt.addValue("X (COG)", cog_x);
				rt.addValue("Y (COG)", cog_y);
				
				
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
				
				AbstractTrajectoryFeature pwf = new PowerLawFeature(t, 1, t.size()/3);
				double[] res = pwf.evaluate();
				v = res[0];
				rt.addValue("ALPHA", v);
				
				GaussianityFeauture gauss = new GaussianityFeauture(t, 1);
				v = gauss.evaluate()[0];
				rt.addValue("GAUSSIANITY", v);
				
				Asymmetry3Feature asym3 = new Asymmetry3Feature(t);
				v = asym3.evaluate()[0];
				rt.addValue("Asymmetry", v);
				
				MSDRatioFeature msdratio = new MSDRatioFeature(t, 1, 5);
				v = msdratio.evaluate()[0];
				rt.addValue("MSDRatio", v);
				
				
				ShortTimeLongTimeDiffusioncoefficentRatio ltstratio = new ShortTimeLongTimeDiffusioncoefficentRatio(t, 2);
				v = ltstratio.evaluate()[0];
				rt.addValue("LTSTRatio", v);

				
		}
		//rt.show("Results");
		parents.show("Parents");
		Iterator<String> rtIt = rtables.keySet().iterator();
		while(rtIt.hasNext()){
			String rt = rtIt.next();
			rtables.get(rt).show(rt + " trajectories");
		}
		
		IJ.log("Done");
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
	
	public ArrayList<Subtrajectory> getClassifiedTrajectories(){
		return classifiedTrajectories;
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




}
