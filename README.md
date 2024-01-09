# relaxingreplication
Replication archive for the paper "Relaxing the Exclusion Restriction in Shift-Share Instrumental Variable Estimation"

# R
00Main.R - Loads packages and runs programs.
01Functions_repo.R - Contains functions needed to run simulations.
02Simulations.R - Runs simulations for single regressors and creates graphs.
03SimulationsMultiple.R - Simulations from Online Appendix, for multiple regressors.
04BP-Shares.R - ALasso and CIM selection when using Anderson-Rubin test for main analysis.
05BP-Shifts.R - ALasso and CIM selection when using Anderson-Rubin test for auxiliary analysis (in online appendix) where multiple shift-share IVs available.
06Visualizations.R - Creates visualization shown in the paper. 
07VisualizationsMultiple.R - Creates visualizations and calculates upper bounds of number of invalid IVs for multiple regressor case (mostly in Appendix).

# Stata
00Programs-Final.do - Includes ssada and sscim programs needed in applications.
01Analysis-Shares-Final.do - Runs main analysis for single and multiple regressor case and creates tables.
Shifts*.do - Creates the shift-share variables to be used in 02Analysis-Shifts.do. Data sources for raw data can be found in the text. 
02Analysis-Shifts.do - Runs the additional results for multiple shift-share variables, to be found in the Appendix.
sscim.ado - Program that runs sscim. Description of syntax can be found at the end of the online appendix on journal homepage.