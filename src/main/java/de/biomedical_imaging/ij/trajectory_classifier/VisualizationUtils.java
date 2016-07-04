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
	
	public static ArrayList<Roi> generateVisualizationRoisFromTrack(Subtrajectory t, Color c){
		return generateVisualizationRoisFromTrack(t, c, 1);
	}
	
	public static ArrayList<Roi> generateVisualizationRoisFromTrack(Subtrajectory t, Color c, double pixelsize){
		ArrayList<Roi> proi = new ArrayList<Roi>();
		FloatPolygon p = new FloatPolygon();
		double sumx = 0;
		double sumy = 0;
		TextRoi.setFont("TimesRoman", 8, Font.PLAIN);
		for(int i = 0; i < t.getParent().size(); i++){
			int to = t.size();
			if(i< t.size()){
				sumx += t.get(i).x/pixelsize;
				sumy += t.get(i).y/pixelsize;
				p.addPoint(t.get(i).x/pixelsize, t.get(i).y/pixelsize);
				
				to = i+1;
			}
			
			PolygonRoi pr = new PolygonRoi(new FloatPolygon(p.xpoints, p.ypoints,to), PolygonRoi.POLYLINE);
			pr.setStrokeColor(c);
			pr.setPosition(t.getRelativeStartTimepoint()+i+1);
			proi.add(pr);
			
			TextRoi troi = new TextRoi(sumx/to, sumy/to," "+t.getID()+" ");
			troi.setPosition(t.getRelativeStartTimepoint()+i+1);
			troi.setFillColor(Color.BLACK);
			troi.setStrokeColor(c);
			troi.setAntialiased(true);
			
			proi.add(troi);
		}
		return proi;
	}

}
