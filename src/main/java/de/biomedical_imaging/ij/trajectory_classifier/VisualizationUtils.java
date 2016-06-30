package de.biomedical_imaging.ij.trajectory_classifier;

import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.process.FloatPolygon;

import java.awt.Color;
import java.awt.Font;
import java.util.ArrayList;

import de.biomedical_imaging.traJ.Trajectory;

public class VisualizationUtils {
	
	public static ArrayList<Roi> generateVisualizationRoisFromTrack(Trajectory t, Color c){
		return generateVisualizationRoisFromTrack(t, c, 1);
	}
	
	public static ArrayList<Roi> generateVisualizationRoisFromTrack(Trajectory t, Color c, double pixelsize){
		ArrayList<Roi> proi = new ArrayList<Roi>();
		FloatPolygon p = new FloatPolygon();
		double sumx = 0;
		double sumy = 0;
		TextRoi.setFont("TimesRoman", 8, Font.PLAIN);
		for(int i = 0; i < t.size(); i++){
			System.out.println("x " + t.get(i).x);
			sumx += t.get(i).x/pixelsize;
			sumy += t.get(i).y/pixelsize;
			p.addPoint(t.get(i).x/pixelsize, t.get(i).y/pixelsize);
			PolygonRoi pr = new PolygonRoi(new FloatPolygon(p.xpoints, p.ypoints,i+1), PolygonRoi.POLYLINE);
			pr.setStrokeColor(c);
			pr.setPosition(t.getRelativeStartTimepoint()+i+1);
			proi.add(pr);
			
			TextRoi troi = new TextRoi(sumx/(i+1), sumy/(i+1),""+t.getID());
			troi.setPosition(t.getRelativeStartTimepoint()+i+1);
			troi.setStrokeColor(c);
			
			proi.add(troi);
		}
		return proi;
	}

}
