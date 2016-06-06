package de.biomedical_imaging.ij.trajectory_classifier;

import ij.IJ;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;


public class TrackMateImporter {
	
	public ArrayList<Track> importTrackMateXML(String path){
		
		ArrayList<Track> trajectories = new ArrayList<Track>();
		
		/*
		 * 2. Load xml file
		 */
		Document doc = null;
		try {
			File fXmlFile = new File(path);
			DocumentBuilderFactory dbFactory = DocumentBuilderFactory
					.newInstance();
			DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
			doc = dBuilder.parse(fXmlFile);

		} catch (SAXException e) {
			// TODO Auto-generated catch block
			IJ.log("" + e.getStackTrace());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			IJ.log("" + e.getStackTrace());
		} catch (ParserConfigurationException e) {
			// TODO Auto-generated catch block
			IJ.log("" + e.getStackTrace());
		}
		
		/*
		 * 3. Build datastructure
		 */
	
		NodeList nTracks = doc.getElementsByTagName("particle");

		for (int i = 0; i < nTracks.getLength(); i++) {
			Track t = new Track(2);
			Node track = nTracks.item(i);
			boolean firstPosition = true;
			NodeList nSteps = track.getChildNodes();
			for (int j = 0; j < nSteps.getLength(); j++) {
				Node step = nSteps.item(j);
				if (step.getNodeName() == "detection") {
					Element e = (Element) step;
					double x = Double.parseDouble(e.getAttribute("x"));
					double y = Double.parseDouble(e.getAttribute("y"));
					int time = Integer.parseInt(e.getAttribute("t"));
					if(firstPosition){
						t.setRelativStartTimepoint(time);
						firstPosition = false;
					}
					t.add(x,y,0);
				}

			}
			trajectories.add(t);
		}

		return trajectories;
	}

}
