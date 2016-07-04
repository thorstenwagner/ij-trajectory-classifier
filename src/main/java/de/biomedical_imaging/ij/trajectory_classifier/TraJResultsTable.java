package de.biomedical_imaging.ij.trajectory_classifier;

import java.awt.Color;
import java.awt.Menu;
import java.awt.MenuItem;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.knowm.xchart.Chart;

import de.biomedical_imaging.traJ.Trajectory;
import de.biomedical_imaging.traJ.TrajectoryUtil;
import de.biomedical_imaging.traJ.VisualizationUtils;
import de.biomedical_imaging.traJ.DiffusionCoefficientEstimator.RegressionDiffusionCoefficientEstimator;
import de.biomedical_imaging.traJ.features.ConfinedDiffusionParametersFeature;
import de.biomedical_imaging.traJ.features.PowerLawFeature;
import ij.IJ;
import ij.WindowManager;
import ij.measure.ResultsTable;
import ij.text.TextPanel;

public class TraJResultsTable extends ResultsTable {
	
	@Override
	public void show(String windowTitle) {
		// TODO Auto-generated method stub
		super.show(windowTitle);
		final String name = windowTitle;
		final ResultsTable table = this;
		
		
		/*
		 * Plot trajectory
		 */
		MenuItem plotSelectedTrajectory = new MenuItem("Plot trajectory");
		plotSelectedTrajectory.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
			//	WindowManager.getFrame(windowTitle)
				//WindowManager.getFrame(windowTitle)
				
				if (WindowManager.getFrame(name).getComponent(0) instanceof TextPanel){
					TextPanel p = (TextPanel) WindowManager.getFrame(name).getComponent(0);
				
					int selectionStart = p.getSelectionStart();
					int selectionEnd = p.getSelectionEnd();
					ArrayList<Chart> charts = new ArrayList<Chart>();
					if(selectionStart>=0 && selectionStart==selectionEnd){
						int id = (int) table.getValue("ID", selectionStart);
						ArrayList<? extends Trajectory> cTracks = TraJClassifier_.getInstance().getClassifiedTrajectories();
						Trajectory t = TrajectoryUtil.getTrajectoryByID(cTracks, id);
						Chart c = VisualizationUtils.getTrajectoryChart("Trajectory with ID " + id,t);
						charts.add(c);
						if(t.getType().equals("SUBDIFFUSION")){
							PowerLawFeature pwf = new PowerLawFeature(t, 1, t.size()/3);
							double[] res = pwf.evaluate();
							c = VisualizationUtils.getMSDLineWithPowerModelChart(t, 1, t.size()/3, TraJClassifier_.getInstance().getTimelag(), res[0], res[1]);
							charts.add(c);
							

						}else if(t.getType().equals("CONFINED")){
							ConfinedDiffusionParametersFeature cfeature = new ConfinedDiffusionParametersFeature(t, TraJClassifier_.getInstance().getTimelag());
							double[] res= cfeature.evaluate();
							c = VisualizationUtils.getMSDLineWithConfinedModelChart(t, 1, t.size()/3, TraJClassifier_.getInstance().getTimelag(), res[0], res[2], res[3], res[1]);
							charts.add(c);
						}else if(t.getType().equals("NORM. DIFFUSION")){
							RegressionDiffusionCoefficientEstimator regest = new RegressionDiffusionCoefficientEstimator(t, 1/TraJClassifier_.getInstance().getTimelag(), 1, t.size()/3);
							double[] res =regest.evaluate();
							double dc = res[0];
							double intercept = res[2];
							System.out.println("RES2 " + res[2]);
							c = VisualizationUtils.getMSDLineWithFreeModelChart(t, 1, t.size()/3, TraJClassifier_.getInstance().getTimelag(), dc, intercept);
							charts.add(c);
						}else{
							c = VisualizationUtils.getMSDLineChart(t, 1, t.size()/3);
							charts.add(c);
						}
						VisualizationUtils.plotCharts(charts);
					}else if( selectionStart!=selectionEnd){
						IJ.showMessage("Plot of multiple trajectories is not possible");
					}
					else{
						IJ.showMessage("No trajectory selected");
					}
				}
				
			}
		});
		
		/*
		 * Export trajectory(s)
		 */
		MenuItem exportTrajectories = new MenuItem("Export trajetorie(s)");
		
		exportTrajectories.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				if (WindowManager.getFrame(name).getComponent(0) instanceof TextPanel){
					TextPanel p = (TextPanel) WindowManager.getFrame(name).getComponent(0);
				
					int selectionStart = p.getSelectionStart();
					int selectionEnd = p.getSelectionEnd();
					ArrayList<Trajectory> selectedTrajectories = new ArrayList<Trajectory>();
					for( int i = selectionEnd; i <= selectionEnd; i++){
						int id = (int) table.getValue("ID", selectionStart);
						ArrayList<? extends Trajectory> cTracks = TraJClassifier_.getInstance().getClassifiedTrajectories();
						Trajectory t = TrajectoryUtil.getTrajectoryByID(cTracks, id);
						selectedTrajectories.add(t);
					}
					
					JFileChooser chooser=new JFileChooser();
					chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
					chooser.addChoosableFileFilter(new FileNameExtensionFilter("Comma seperated value (.csv)", "csv"));
					chooser.setAcceptAllFileFilterUsed(false);
					int c =chooser.showSaveDialog(null);
					if(c == JFileChooser.APPROVE_OPTION){
						String path=chooser.getSelectedFile().getAbsolutePath();
						if(!path.substring(path.length()-3, path.length()).equals("csv")){
							path += ".csv";
						}
				
						ExportImportTools eit = new ExportImportTools();
						eit.exportTrajectoryDataAsCSV(selectedTrajectories, path);
					}
					
				}
				
				
				
			}
		});
		
		Menu traJ = new Menu("TraJ");
		traJ.add(plotSelectedTrajectory);
		traJ.add(exportTrajectories);
				//ResultsTable.getResultsWindow().getMenuBar().
		//Hinzufügen von Export Funktionen
		WindowManager.getFrame(windowTitle).getMenuBar().add(traJ);
	}

}
