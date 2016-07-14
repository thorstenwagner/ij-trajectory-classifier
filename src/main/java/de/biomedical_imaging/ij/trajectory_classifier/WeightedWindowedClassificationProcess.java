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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.stat.Frequency;
import org.apache.commons.math3.stat.descriptive.moment.Mean;

import de.biomedical_imaging.traJ.Trajectory;

public class WeightedWindowedClassificationProcess {
	
	public String[] windowedClassification(Trajectory t, AbstractClassifier c, int n){
		int windowsize = 2*n+1;
		
		int increment = 1;
		ArrayList<Trajectory> tracks = new ArrayList<Trajectory>();
		for(int i = 0; i < (t.size()-windowsize+increment);i=i+increment){
			Trajectory sub = t.subList(i, i+windowsize-1);
			tracks.add(sub);
		}
		
	
		String[] res = c.classify(tracks);
		double[] confidence = c.getConfindence();
		
		
		
		
		
		String[] types = applyWeightening(res, confidence, n, t.size());
		return types;
	}
	
	protected String[] applyWeightening(String[] res, double[] confidence, int n, int tracklength){
		
		String[] types = new String[tracklength];
		for(int i = 0; i < types.length; i++){
			types[i] = "NONE";
		}
		
		double start= System.currentTimeMillis();
		//Build mapping
				ArrayList<String> restypes = new ArrayList<String>();
				for(int i =0; i < res.length; i++){
					restypes.add(res[i]);
				}
				restypes.add("NONE");
				HashSet<String> uniqueTypes = new HashSet<String>();
				uniqueTypes.addAll(restypes);
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
				
				ArrayList<Double[]> weightes = new  ArrayList<Double[]>();
				for(int i = 0; i < tracklength; i++){
					Double[] h = new Double[key];
					Arrays.fill(h, new Double(0));
					weightes.add(h);
				}
				
	
				
				for(int i = 0; i < res.length; i++){
					for(int j = i; j < i+2*n+1;j++){
						int typ = mapTypeToInt.get(res[i]);


						weightes.get(j)[typ]=weightes.get(j)[typ]+confidence[i];
					}
				}

				for(int i = 0; i < types.length; i++){
					if(weightes.get(i).length>0){
					//	int[] help = arrListToArrayInt(weightes.get(i));
					//	Arrays.sort(help);
						//int[] modes = getMode(help,0,help.length);
				
				
						int mode1 = getHighest(weightes.get(i));
						//if(v!=(int)modes[0]){
						//	System.out.println("!!! v: " + v + " modes" + modes[0]);
						//}
						String mode = mapIntToType.get(mode1);
						types[i] = mode;
				
					}
					else{
						types[i] = "NONE";
					}
				}
				
				return types;
	}
	
	private int getHighest(Double[] weightes){
		double max = 0;
		int maxindex = 0;
		for(int i = 0; i < weightes.length;i++){
			if(weightes[i]>max){
				max = weightes[i];
				maxindex = i;
			}
		}
		return maxindex;
		
	}
	
	private int[] getMode(int[] values, final int begin, final int length) {
        // Add the values to the frequency table
        Frequency freq = new Frequency();
        for (int i = begin; i < begin + length; i++) {
            final int value = values[i];
                freq.addValue(Integer.valueOf(value));
        }
        List<Comparable<?>> list = freq.getMode();
        // Convert the list to an array of primitive double
        int[] modes = new int[2];
        modes[0] = ((Long) list.get(0)).intValue();
      //  modes[1] = freq.getCount(Integer.valueOf(modes[0]))/(1.0*length);
        return modes;
    }
	
	public int findPopular(int[] a) {

	    if (a == null || a.length == 0)
	        return 0;

	    Arrays.sort(a);

	    int previous = a[0];
	    int popular = a[0];
	    int count = 1;
	    int maxCount = 1;

	    for (int i = 1; i < a.length; i++) {
	        if (a[i] == previous)
	            count++;
	        else {
	            if (count > maxCount) {
	                popular = a[i-1];
	                maxCount = count;
	            }
	            previous = a[i];
	            count = 1;
	        }
	    }

	    return count > maxCount ? a[a.length-1] : popular;

	}
	
	public double[] arrListToArray(ArrayList<Integer> l){
		double[] a = new double[l.size()];
		for(int i = 0; i < l.size(); i++){
			a[i] = l.get(i).intValue();
		}
		return a;
	}
	
	public int[] arrListToArrayInt(ArrayList<Integer> l){
		int[] a = new int[l.size()];
		for(int i = 0; i < l.size(); i++){
			a[i] = l.get(i).intValue();
		}
		return a;
	}

}
