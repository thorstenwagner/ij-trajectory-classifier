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

import de.biomedical_imaging.traJ.features.AbstractTrajectoryFeature;

public class FeatureWorker extends Thread {
	enum EVALTYPE{
		FIRST,RATIO_01,RATIO_10,RATIO_12;
	}
	double[] result;
	AbstractTrajectoryFeature c;
	EVALTYPE ev;
	int resIndex;
	
	public FeatureWorker(double[] result, int resIndex, AbstractTrajectoryFeature c, EVALTYPE ev) {
		this.result = result;
		this.c = c;
		this.ev =ev;
		this.resIndex = resIndex;
	}
	
	@Override
	public void run() {
		double[] res ;
			switch (ev) {
			case FIRST:
				result[resIndex] = c.evaluate()[0];
				break;
			case RATIO_01:
				res = c.evaluate();
				result[resIndex] = res[0]/res[1];
				break;
			case RATIO_10:
				res = c.evaluate();
				result[resIndex] = res[1]/res[0];
				break;
			case RATIO_12:
				res = c.evaluate();
				result[resIndex] = res[1]/res[2];
				break;

			default:
				break;
			}
			
		}

}
