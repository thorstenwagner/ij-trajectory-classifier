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
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;

import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.AbstractDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.CovarianceDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.RegressionDiffusionCoefficientEstimator;
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
import ij.IJ;
import ij.Prefs;
import ij.WindowManager;
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
	private double minDiffusionCoefficient;
	private double pixelsize;
	private boolean showID;
	private boolean showOverviewClasses;
	private ArrayList<Subtrajectory> classifiedTrajectories;
	private ArrayList<Trajectory> tracksToClassify;
	//private ArrayList<Trajectory> tracks 
	private static TraJClassifier_ instance;
	
	
	public TraJClassifier_() {
		minTrackLength=160;
		windowSizeClassification=60;
		minDiffusionCoefficient=0;
		pixelsize=0.166;
		timelag=1.0/30;
		showID = true;
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
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		/*
		 * Import Data
		 */
		if(!arg.contains("DEBUG")){
			OpenDialog open = new OpenDialog("Choose the TrackMate xml file");
			String filepath = open.getPath();
			TrackMateImporter tMateImport = new TrackMateImporter();
			tracksToClassify = tMateImport.importTrackMateXML(filepath);
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
		
		if(!arg.contains("NOGUI")){
			//Load previous settings
			minTrackLength = Prefs.get("trajclass.minTrackLength", 160);
			windowSizeClassification = (int) Prefs.get("trajclass.windowSize", 60);
			minDiffusionCoefficient = Prefs.get("trajclass.minDC", 0);
			pixelsize = Prefs.get("trajclass.pixelsize", 0.166);
			timelag = 1.0/Prefs.get("trajclass.framerate",30);
			showID = Prefs.getBoolean("trajclass.showID", true);
			showOverviewClasses = Prefs.getBoolean("trajclass.showOverviewClasses", true);
			
			//Show GUI
			GenericDialog gd = new GenericDialog("Parameters Classification");
		
			gd.addSlider("Min. tracklength", 1, 1000, minTrackLength);
			System.out.println("window: " + windowSizeClassification);
			gd.addSlider("Windowsize", 1, 200, windowSizeClassification);
			gd.addNumericField("Min. diffusion coeffcient (µm^2 / s)", minDiffusionCoefficient, 0);
			gd.addNumericField("Pixelsize (µm)*", pixelsize, 4);
			gd.addNumericField("Framerate", 1/timelag, 0);
			gd.addCheckbox("Show IDs", showID);
			gd.addCheckbox("Show overview classes", showOverviewClasses);
			gd.addMessage("* Set to zero if the imported data is already correctly scaled.");
			gd.addHelp("http://forum.imagej.net");
			gd.showDialog();
			minTrackLength = gd.getNextNumber();
			windowSizeClassification = (int) (gd.getNextNumber()/2);
			minDiffusionCoefficient = gd.getNextNumber();
			pixelsize = gd.getNextNumber();
			timelag = 1/gd.getNextNumber();
			showID = gd.getNextBoolean();
			showOverviewClasses = gd.getNextBoolean();
			
			// Save settings
			Prefs.set("trajclass.minTrackLength", minTrackLength);
			Prefs.set("trajclass.windowSize", windowSizeClassification*2);
			Prefs.set("trajclass.minDC", minDiffusionCoefficient);
			Prefs.set("trajclass.pixelsize", pixelsize);
			Prefs.set("trajclass.framerate", 1/timelag);
			Prefs.set("trajclass.showID", showID);
			Prefs.set("trajclass.showOverviewClasses", showOverviewClasses);
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
		AbstractClassifier rrf = new RRFClassifierRenjin(modelpath);
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
			if(track.size() > minTrackLength){
				minLengthTracks.add(track);
			}
		}
		
		/*
		 * Classification and segmentation
		 */
		classifiedTrajectories = new ArrayList<Subtrajectory>();
		int subidcounter = 1;
		for (Trajectory track : minLengthTracks) {
			j++;
			IJ.showProgress(j, minLengthTracks.size());

			String[] classes = wcp.windowedClassification(track, rrf, windowSizeClassification);
			
			Subtrajectory tr = new Subtrajectory(track,2);
			tr.setID(subidcounter);
			subidcounter++;
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
					tr.setID(subidcounter);
					subidcounter++;
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
		
		//Segments have to contain at least 30 steps
		for(int i = 0; i < classifiedTrajectories.size(); i++){
			if(classifiedTrajectories.get(i).size()<30){
				
			}
		}
	
		//Segments smaller than the window size seems to be very error prone. Set them to "None".
		for (Trajectory trajectory : classifiedTrajectories) {
			if(trajectory.size()<windowSizeClassification*2){
				trajectory.setType("NONE");
			}
		}
		
		/*
		 * Visualization 
		 */
		if (visualize) {
			// Trajectories
			Overlay ov = new Overlay();
			for (int i = 0; i < classifiedTrajectories.size(); i++) {
				Subtrajectory tr = classifiedTrajectories.get(i);

				ArrayList<Roi> prois = null;
				if (pixelsize > 0.000001) {
					prois = VisualizationUtils
							.generateVisualizationRoisFromTrack(tr,
									mapTypeToColor.get(tr.getType()), showID,
									pixelsize);

				} else {
					prois = VisualizationUtils
							.generateVisualizationRoisFromTrack(tr,
									mapTypeToColor.get(tr.getType()), showID);
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
				TextRoi.setFont("TimesRoman", 12, Font.PLAIN);

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
		 * Export classified trajectories
		 */
		//IJ.log("Export");
		//ExportImportTools iotools = new ExportImportTools();
		//iotools.exportTrajectoryDataAsCSV(classifiedTrajectories, "/home/thorsten/myclassifiedtracks.csv");
		
		/*
		 * Fill results table
		 */
		
		
		HashMap<String, TraJResultsTable> rtables = new HashMap<String, TraJResultsTable>();
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
				TraJResultsTable rt = rtables.get(t.getType());
				
				rt.incrementCounter();
				rt.addValue("ID", t.getID());
				rt.addValue("PARENT-ID", t.getParent().getID());
				rt.addValue("LENGTH", t.size());
				rt.addValue("START", t.getRelativeStartTimepoint());
				rt.addValue("END", t.getRelativeStartTimepoint()+t.size()-1);
				rt.addValue("CLASS", t.getType());

				AbstractTrajectoryFeature dcEstim=null;
				double dc =0;
				switch (t.getType()) {
				case "DIRECTED/ACTIVE":
					dcEstim = new RegressionDiffusionCoefficientEstimator(t,1/timelag,1,t.size()-1);
					dc = dcEstim.evaluate()[0];
					rt.addValue("D", String.format("%6.3e",dc));
					
					break;
				case "NORM. DIFFUSION":
					dcEstim = new CovarianceDiffusionCoefficientEstimator(t, 1/timelag);
					dc = dcEstim.evaluate()[0];
					rt.addValue("D", String.format("%6.3e",dc));
					break;
				case "CONFINED":
					AbstractDiffusionCoefficientEstimator dcEst = new RegressionDiffusionCoefficientEstimator(t,1/timelag,1,3);
					ConfinedDiffusionParametersFeature confp = new ConfinedDiffusionParametersFeature(t,timelag,dcEst);
					double[] p = confp.evaluate();
					dc = p[1];
					rt.addValue("CONF. SIZE", p[0]);
					rt.addValue("A (CONF SHAPE)", p[2]);
					rt.addValue("B (CONF SHAPE)", p[3]);
					rt.addValue("D", String.format("%6.3e",p[1]));
					break;
				case "SUBDIFFUSION":
					PowerLawFeature pwf = new PowerLawFeature(t, 1, t.size()/3);
					double res[] = pwf.evaluate();
					dc = res[1];
					rt.addValue("D", String.format("%6.3e",dc));
					break;
				case "NONE":
					break;
				default:
					break;
				}
				
				if(dc<minDiffusionCoefficient){
					t.setType("STALLED");
				}
				
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

		}
		
		ResultsTable parents = new ResultsTable();
		for(int i = 0; i < minLengthTracks.size(); i++){
			parents.incrementCounter();
			Trajectory t = minLengthTracks.get(i);
			parents.addValue("ID", t.getID());
			parents.addValue("LENGTH", t.size());
			parents.addValue("START", t.getRelativeStartTimepoint());
			parents.addValue("END", t.getRelativeStartTimepoint()+t.size()-1);
			int subCount =0;
			int normCount =0;
			int directedCount =0;
			int confCount=0;
			int stalledCount=0;
			int noneCount =0;
			ArrayList<Subtrajectory> sameParent = Subtrajectory.getTracksWithSameParant(classifiedTrajectories, t.getID());
			for (Subtrajectory sub : sameParent) {
				switch (sub.getType()) {
				case "DIRECTED/ACTIVE":
					directedCount+=sub.size();
					break;
				case "NORM. DIFFUSION":
					normCount+=sub.size();
					break;
				case "CONFINED":
					confCount+=sub.size();
					break;
				case "SUBDIFFUSION":
					subCount+=sub.size();
					break;
				case "STALLED":
					stalledCount+=sub.size();
					break;
				case "NONE":
					noneCount+=sub.size();
					break;
				default:
					break;
				}
			}
			parents.addValue("#POS_NORM", normCount);
			parents.addValue("#POS_SUB", subCount);
			parents.addValue("#POS_CONF", confCount);
			parents.addValue("#POS_DIRECTED", directedCount);
			parents.addValue("#POS_STALLED", stalledCount);
			parents.addValue("#POS_NONE", noneCount);
		}
		
		
		// show tables
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
    

	public double getMinTrackLength() {
		return minTrackLength;
	}

	public void setMinTrackLength(double minTrackLength) {
		this.minTrackLength = minTrackLength;
	}
	
	public double getMinDiffusionCoefficient() {
		return minDiffusionCoefficient;
	}

	public void setMinDiffusionCoefficient(double minDiffusionCoefficient) {
		this.minDiffusionCoefficient = minDiffusionCoefficient;
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

	public void setTimelag(double timelag) {
		this.timelag = timelag;
	}

	
	
	public void setWindowSizeClassification(int windowSizeClassification){
		this.windowSizeClassification = windowSizeClassification;
	}




}
