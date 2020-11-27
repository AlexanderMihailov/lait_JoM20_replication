readme_lait_JoM20_replictaion.txt

Default path: [F:\archive\comp\ on the external drive used, but can be anything] lait_JoM20_replication\
First created: 17 July 2015, by Alexander Mihailov
Last modified: 16 January 2020, by Alexander Mihailov
(to be printed out in Landscape orientation, 110-character length of text observed below)

This readme file describes the folder structure in the zip archive with replication files for McKnight,
Mihailov and Pompa Rangel, "What do Latin American inflation targeters care about? A comparative Bayesian
estimation of central preferences", Journal of Macroeconomics, Volume 63, March 2020, Article 103188.
__________________________________________________________________________________________________________

Main folder: \lait_JoM20_replication - contains 5 subfolders (as described below)

\chain - contains MATLAB files used in computing a range of diagnostics for the MCMC convergence

\data -	our dataset in *.xlsx (Excel spreadsheet) format: quarterly time series for the 5 LAIT economies:
	Brazil (bra), Chile (chi), Colombia (col), Mexico (mex) and Peru (per)
 
\func - contains MATLAB functions and code that solve the model and calculate marginal likelihood

\output - contains (in further respective subfolders) the output from running the programs, for each
	country and model version

	\bra
		\... [4 model version subfolders]
	\chi	
		\... [4 model version subfolders]
	\col
		\... [4 model version subfolders]
	\mex	
		\... [4 model version subfolders]		
	\per	
		\... [5 model version subfolders]

\redpolicy - contains MATLAB functions and code that solves for monetary policy under discretion

The main folder also contains 4 MATLAB programs (model versions, as in the paper) for each country.
These programs are annotated in much detail so that reading through explains what each does and how.

The main folder contains as well some other MATLAB files, and laitJoMReplicationGuide.pdf provides
more detail on these and on how the code works and what the components of the respective file (and
folder) names mean, taking as illustration the case of Mexico.
