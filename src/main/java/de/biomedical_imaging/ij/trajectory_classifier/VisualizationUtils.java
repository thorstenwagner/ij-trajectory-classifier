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

import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.process.FloatPolygon;

import java.awt.Color;
import java.awt.Font;
import java.util.ArrayList;


public class VisualizationUtils {
	
	public static ArrayList<Roi> generateVisualizationRoisFromTrack(Subtrajectory t, Color c, boolean showID){
		return generateVisualizationRoisFromTrack(t, c, showID,1);
	}
	
	public static ArrayList<Roi> generateVisualizationRoisFromTrack(Subtrajectory t, Color c, boolean showID, double pixelsize){
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
			
			if(showID){
				long parentID = t.getParent().getID();
				TextRoi troi = new TextRoi(sumx/to, sumy/to," "+parentID+":"+t.getID()+" ");
				troi.setPosition(t.getRelativeStartTimepoint()+i+1);
				troi.setFillColor(Color.BLACK);
				troi.setStrokeColor(c);
				troi.setAntialiased(true);
				proi.add(troi);
			}
		}
		return proi;
	}

}
