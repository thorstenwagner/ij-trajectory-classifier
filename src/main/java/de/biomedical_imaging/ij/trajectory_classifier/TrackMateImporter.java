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

import de.biomedical_imaging.traJ.Trajectory;


public class TrackMateImporter {
	
	public ArrayList<Trajectory> importTrackMateXML(String path){
		
		ArrayList<Trajectory> trajectories = new ArrayList<Trajectory>();
		Trajectory.restIDCounter();
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
			Trajectory t = new Trajectory(2);
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
