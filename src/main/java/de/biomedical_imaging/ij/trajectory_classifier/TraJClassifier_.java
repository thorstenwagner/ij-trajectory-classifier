package de.biomedical_imaging.ij.trajectory_classifier;
import java.awt.Color;
import java.awt.Font;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.TextRoi;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;

public class TraJClassifier_ implements PlugIn {


	
	@Override
	public void run(String arg) {
		
		/*
		 * GUI
		 */
		GenericDialog gd = new GenericDialog("Parameters Classification");
		gd.addSlider("Min. tracklength", 160, 1000, 160);
		gd.addNumericField("FPS", 30, 0);
		gd.showDialog();
		double minlength = gd.getNextNumber();
		double fps = 1.0/gd.getNextNumber();
		/*
		 * Gut funktioniert hat bisher wsc = 30, ms = 80
		 */
		int windowSizeClassification = 30; // -> 2*n+1=121
		int modeSize = 80;	// -> 2*n+1 = 81
		
		/*
		 * Import Data
		 */

		OpenDialog open = new OpenDialog("Choose the TrackMate xml file");
		String filepath = open.getPath();
		
		
		TrackMateImporter tMateImport = new TrackMateImporter();
		ArrayList<Track> tracks = tMateImport.importTrackMateXML(filepath);
		
		
		/*
		 * Classification & Visualization
		 */
		AbstractClassifier rrf = new RRFClassifierParallel(fps);
		rrf.start();
		WindowedClassificationProcess wcp = new WindowedClassificationProcess();
		Overlay ov = new Overlay(); 
		int j = 0;

		
		HashMap<String, Color> mapTypeToColor = new HashMap<String, Color>();
		mapTypeToColor.put("DIRECTED/ACTIVE", Color.YELLOW);
		mapTypeToColor.put("NORM. DIFFUSION", Color.RED);
		mapTypeToColor.put("SUBDIFFUSION", Color.GREEN);
		mapTypeToColor.put("NONE", Color.LIGHT_GRAY);

		
		for (Track track : tracks) {
			ArrayList<PolygonRoi> prArr = new ArrayList<PolygonRoi>();
			j++;
			if(track.size() > minlength){
				IJ.log("Now at: " + j + " / " + tracks.size());
		
				String[] classes = wcp.windowedClassification(track, rrf, windowSizeClassification);
				classes = movingMode(classes, modeSize);
				
			
			
	
			
				FloatPolygon p = new FloatPolygon();
				p.addPoint(track.get(0).x, track.get(0).y);
				String prevCls = classes[0];
				for(int i = 1; i < classes.length; i++){
					if(prevCls == classes[i]){
						p.addPoint(track.get(i).x, track.get(i).y);
						PolygonRoi pr = new PolygonRoi(p,PolygonRoi.POLYLINE);
						pr.setStrokeColor(mapTypeToColor.get(prevCls));
						pr.setPosition(track.getRelativeStartTimepoint()+i);
						for (PolygonRoi roi : prArr) {
							PolygonRoi roi2 = new PolygonRoi(roi.getFloatPolygon(), PolygonRoi.POLYLINE);
							roi2.setPosition(track.getRelativeStartTimepoint()+i);
							roi2.setStrokeColor(roi.getStrokeColor());
							ov.add(roi2);
							
						}
						ov.add(pr);
						
					}else{
						PolygonRoi pr = new PolygonRoi(p,PolygonRoi.POLYLINE);
						pr.setStrokeColor(mapTypeToColor.get(prevCls));
						pr.setPosition(track.getRelativeStartTimepoint()+i-1);
						prArr.add(pr);
						p = new FloatPolygon();
						p.addPoint(track.get(i).x, track.get(i).y);
						prevCls = classes[i];
						
					}
				}
			}
		}
		rrf.stop();
		IJ.log("Plot!");
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
