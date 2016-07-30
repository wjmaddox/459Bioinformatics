# EECS 459 Bioinformatics for Systems Biology Final Project
Title: Finding Important Proteins in Exocytosis and Endocytosis Protein Networks
Professor Mehmet Koyuturk, Spring 2015 

Project was performed with Dong Liu who provided biological domain expertise.

Project goal was to find rich-clubs in proteins in exocytosis and endocytosis protein interaction networks. Utilized two different measures of measuring and finding rich-clubs 
* FPC scoring: 
	Pei Wang, Xinghuo Yu & Jinhu Lu. Identification and Evolution of Structurally Dominant Nodes in Protein-Protein Interaction Networks. IEEE Trans. Biomed. Circuits Syst. 8, 87–97 (2014).
* Rich-club: 
	Rombach, M. P., Porter, M. A., Fowler, J. H. & Mucha, P. J. Core-Periphery Structure in Networks. SIAM J. Appl. Math. 74, 167–190 (2014).

This is a slightly strange set-up for the Git repository because all of the validation runs were done on Case's HPC cluster.
So the folder Server is directly from what was uploaded to the cluster in order to be able to run these validation runs in a reasonable amount of time.

Everything else was on a Windows 8 laptop with an Intel i5 processor.

The report is available upon request.

##Structure 
###Code folder

GeneratePPI.R - Merges and creates a protein protein interaction network (igraph) including both endocytosis and exocytosis from data downloaded from Gene Ontology and BIOGRID databases.

GeneratePartialPPI.R -  Generates endocytosis and exocytosis only PPI networks.

ExoEndoFPC.R - Performs FPC scoring on all three networks.

RCCFns.R - Helper fucntions to perform the rich club analysis as described in Rombach et al.

RCCExample.R - Test code of both analyses on a very simple graph.

RCCAll.R - Runs rich club analysis in parallel on a cluster environment.

MakePlots.R - Generates all plots from the analysis.

###Main Folder
The rest of the files are either RData encodings of the output from the analysis, or eps files containing the images used in our report.




