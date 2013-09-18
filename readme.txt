RSA normalization Project

----------------------------------------------------------------------------------------------------------------

Table of Contents:

I. Important Programs
	A. Mine PDB information
		1. parse_alignment.py		
		2. get_PDB.py
	B. Theoretical Model Construction
		2. Geometry.py
		3. PeptideBuilder.py
	C. Theoretical data generation
		1. iterateThroughModels.py
		3. DSSPData.py
	D. Data Analysis
		1. max_bins_with_population_restriction.py
		2. max_bins_with_population_restriction_theoretical.py
		3. SeperateOverAndUnderRSA1.py
		4. EmpVCalc_get_Diff_with_pop_restriction.py
		5. max_bin_all_data.py
		6. getRoseRSA.py
		7. make_ALLOWED_GeoFiles.py
		8. get_ALLOWED_bins.py
		9. get_CORE_bins.py
		10. get_GENEROUS_bins.py
	E. R scripts
		1. getMaximumValues.r
		2. getMeanSA.r
		3. CorrelatonTableNewScales.r
		4. get_population_cut_offs.r
		4. FigureScripts
			a. barGRSA.r
			b. makeRSAdistribution.r
			c. makePlotsWithPopRestriction.r
			d. makeRamaPlot.r
			e. makeEmpCalVpop.r
			f. makeNormCorPlot.r
			g. makeALLOWEDBinnedRamaPlot.r
			h. makeCOREBinnedRamaPlot.r
			i. makeGENEROUSBinnedRamaPlot.r
			j. getAngles.r
	F. Misc.
		1. editHydroScales.py
		2. runAll.py

		
II. Data Files
	A. Xxx_geo
	B. AnglesIteratedThroughAgainXXX
	C. XXX_SA_Over/Under, Xxx_Rose_RSA
	D. XXX_max_bins_all
	E. EmpericalVCalculated_diff_pop_nonZeroed_with_pop_restriction_XXX
	F. NormalizationValuesByPercentDataCoverage
	G. NormalizationValuesByPercentDataCoverageAndGenerous.txt
	H. Hydrophobicity_Scales_Updated.txt
	I. Wolfden.txt, rose.txt, Kite_Doolittle.txt, Fauchere.txt, Wimley.txt, Moon.txt, Radzicka.txt, MacCallum.txt 
	J. cullpdb_pc30_res1.8_R0.25_d130607_chains4961.gz
	K. Allowed, Core, and Generous Bins
	

----------------------------------------------------------------------------------------------------------------

I. Important Programs
	This section is dedicated to all python, R, and other scripts used to obtain or analyze data.
	
	A. Mine PDB Information
		To mine the protein structures I parsed the "cullpdb_pc30_res1.8_R0.25_d130607_chains4961.gz" 
		file in "parse_alignment.py".
		
		1. parse_alignment.py
			This is a python script that takes in a list of PDB and Chain ID's (we used the 
			"cullpdb_pc30_res1.8_R0.25_d130607_chains4961.gz" file from the Dunbrak lab).
			This program creates the "Xxx_geo" files (ie. Ala_geo, Asn_geo, etc...).
			This program also only outputs information from non-chain-terminating residues, which
			is determined by the peptide bond length between two residues and non-ambigous neighbors.
			This program is mainly a frame to process the output from the real PDB parser "get_PDB"
			Any files that is either corrupted or did not exist in the PDB database is reported to
			the "Error_report_bond_length.txt".
			The two numbers that are printed to the terminal are for testing.
	
		2. get_PDB.py
			This is a python program that takes in a PDB file name and Chain ID name, downloads
			the PDB file to the "structures" folder than extracts all the information from the
			PDB file. This program needs the PDB parser from Biopython and DSSP and the proper
			parser for it. The program outputs lists of information. This program contains many
			functions inorder to properly mine the data and needs
			DSSPData.py to work.
			
	B. Theoretical Model Construction
		"Geometry" and "PeptideBuilder" are both used to build any more information about their 
		functionality is discussed in. 
		
		1. Geometry.py
			This is a library of Amino Acid Geometry objects. By reading in the one-letter Amino
			Acid abbreviations, it can create the correct geometric parameters to construct a
			protein residue. There are 20 classes (one for each amino acid) and one function
			to create the geometry object. Some of these classes have an inputRotamers function
			that takes in a list of integers in order to change the Amino Acid's rotamers. The
			generateRandomRotamers method is used in "iterateThroughModels"

		2. PeptideBuilder.py
			This takes in geometry objects to construct amino acid chains. There are 20 methods
			to construct the amino acids and has the calculate coordinates method in "makeStructure"
			The program also contains a makeStructure method that takes in a string of Amino Acids,
			Phi list of float, Psi list of float numbers. It contains two add residue methods, an
			initialize_residue method, and a makeExtended Structure method. It also has a output
			structure method. This program creates pdb files of your name choice

	C. Theoretical data generation
		These scripts were used to iterate through the models phi and psi conformations.
		"iterateThroughModels" is a script that uses Geometry and Peptide Builder to build the
		phi and psi conformation, and if possible, all rotamer conformations.
	
		1. iterateThroughModels.py
			The program need the "Geometry" and "PeptideBuilder" to create the phi and psi
			conformation. This has a nice DSSP method and methods to iterate through the psi, psi,
			and chi angles of a residue. This program creates "AnglesIteratedThroughAgainXXX"
			
		2. DSSPData.py
			Parser object to read the output of DSSP program

	D. Data Analysis
		After obtaining the information in the "Xxx_geo" files and rotating through the theoretical
		models, we wrote scripts to parse and analyze the data.
		
		1. max_bins_with_population_restriction.py
			OBSOLETE, replaced by I.D.5
			This program takes in the "Xxx_geo" data file and bins the data into 5-degree by 5-
			degree Phi and Psi coordinates and put in the max SA found for that bin and the number
			of data points in that bin. The input of the program requires an all-caps three letter
			abbreviation of which amino acid you want to look at. This outputs the "XXX_max_
			emperical_bins".

		2. max_bins_with_population_restriction_theoretical.py
			OBSOLETE, replaced by I.D.5
			This program takes in the "XXX_max_emperical_bins_pop_restriction" file and the 
			"AnglesIteratedThroughAgainXXX" file to make the "XXX_max_theoretical_bin_Again" file.
			The command line argument is the all-caps three letter abbreviation of which amino acid 
			you want to look at. This file makes sure that the theoretical data is binned and
			treated in parallel with the empirical data. This program was not used in the final
			write up of the paper, but is a good script nonetheless.

		3. SeperateOverAndUnderRSA1.py
			This program takes in the "Xxx_geo" data file and creates two files "XXX_SA_Over/Under"
			This program is fairly basic, but it's a good basis to parse the "Xxx_geo" files. The
			command line argument is the all-caps three letter abbreviation of which amino acid you
			want to look at.

		4. EmpVCalc_get_Diff_with_pop_restriction.py
			This program takes in the "XXX_geo" file and the "AnglesIteratedThroughAgain" file to
			create the "EmpericalVCalculated_diff_pop_nonZeroed_with_pop_restriction_XXX" file.
			The command line argument is the all-caps three letter abbreviation of which amino acid 
			you want to look at. This program just compares the bins from both "XXX_max_bins" files.

		5. max_bin_all_data.py
			This program takes in the "XXX_geo" files and the "AnglesIteratedThroughAgain" files to
			create the "XXX_max_all_bins" file. The command line argument is the all-caps three letter 
			abbreviation of which amino acid you want to look at. This program compresses and bins 
			SA information from both files into one.

		6. getRoseRSA.py
			This program main use is to make the "Xxx_Rose_RSA" data files, which contain the
			normalized RSA values using the normalization constants from the Rose paper. The 
			input for this program is the "Xxx_geo" files.
			
		7. make_ALLOWED_GeoFiles.py
			This program parses the "XXX_geo" files by the information obtained from the data 
			generated by the "max_bin_all_data.py." The program takes in the "XXX_max_bin_all",
			"XXX_geo", "NormalizationValuesByPercentDataCoverageAndGenerous.txt", to make new
			data files based on being in an ALLOWED region of the Ramachandran plot. This program
			outputs the "XXX_ALLOWED_geo" files and the "AnglesIteratedThroughAgain_ALLOWED_XXX" 
			data files.
			
		8. get_ALLOWED_bins.py
			This program bins the data from the "make_ALLOWED_GeoFiles.py" script.
			
		9. get_CORE_bins.py
			This program bins the data from "XXX_geo" files and "AnglesIteratedThroughAgainXXX" files
			based on the regions defined as CORE in the "NormalizationValuesByPercentDataCoverageAndGenerous.txt"
			output.
		
		10. get_GENEROUS_bins.py
			This program bins the data from "XXX_geo" files and "AnglesIteratedThroughAgainXXX" files
			based on the regions defined as GENEROUS in the 
			"NormalizationValuesByPercentDataCoverageAndGenerous.txt" output.
			
	E. R scripts
		The R-scripts were mainly used to make figures and to do simple things that would have been a 
		bit more complicated to do in python.

		1. getMaximumValues.r
			OBSOLETE
			Obsolete but useful script to look at the data from another perspective
			This gets the maximum SA values from both "XXX_max_bins" data files. It outputs the
			"NormalizationValues.txt" in a csv format.

		2. getMeanSA.r
			This gets the mean RSA, median RSA, square root mean RSA, box-cox transformed mean RSA,
			fraction of 100% buried residues, and fraction of 95% residues (for theoretical 
			normalization values only). For all the average estimates, the script uses both 
			empirical and theoretical normalization values from ALLOWED regions. The output file is a 
			csv file called "Hydrophobicity_Scales_updated.txt". It also has optional 
			"MeanHydrophobicityScales.txt" and "BuriedHydrophobicityScales.txt" output. The script 
			requires the "NormalizationValuesByPercentDataCoverageAndGenerous.txt" file or the 
			"NormalizationValuesByPercentDataCoverage.txt".

		3. CorrelatonTable.r
			This makes the correlation table of all scales in "Hydrophobicity_Scales_updated.txt"
			with the Wolfenden, Kyte, Radzicka, MacCallum, Moon, Wimley, Fauchere, and Rose scales. 
			This require all files to run. This performs the pearson correlation test. 
			
		4. get_population_cut_offs.r
			This script obtains the normalization (MAX RSA) of both the empirical and theoretical data.
			This script is unique in that it needs to be run twice and commented out. First run, is to obtain
			the ALLOWED and CORE bin angle cut-offs. After the ALLOWED and CORE bins are obtained, the program
			is ran again with the uncommented out lines. The GENEROUS bins rely on the areas of the ALLOWED 
			bins and are calculated from the "XXX_ALLOWED_geo" and "AnglesIteratedThroughAgain_ALLOWED_XXX". 

		5. Figure Scripts
			These scripts were used to make the figures in R. most have a pdf/png alternate code in 
			their scripts that are commented out.

			a. barGRSA.r
				Using the "Xxx_geo" files and the "Xxx_Rose_RSA" files, the script makes the 
				"BarGraphRSA.pdf" figure. There is an optional png script at the bottom of this
				file, it is commented out currently.
 
			b. makeRSAdistribution.r
				This script to makes the RSA distribution of Alanine. However, with a bit of tweaking 
				it can make all Amino Acids. This makes the "Alanine_RSA_distribution.pdf" 
				This uses the hist function in R.

			c. makePlotsWithPopRestriction.r
				SEMI-OBSOLETE
				This file makes the best figure. It makes the "XXX_Rama_HSV.svg" figure (the 
				figure that compares the Theoretical and Empirical results. This figure needs
				to be edited in Inkscape. It requires both "XXX_all_bins" files. The file requires
				that you specify "code" which is the three-letter amino acid abbreviation for the
				amino acid you want to make the figure for. This uses the image function to make
				the figures.

			d. makeRamaPlot.r
				This file makes the Ramanchandran plot of all RSA with the Miller normalization 
				values, where the data point are put into two catagories RSA>1 and RSA<=1. This
				script requires the "XXX_SA_Over" and "XXX_SA_Under." The file reads out as 
				"XXX_RamaPlotRSA.pdf". To specify which XXX you want to make, assign the variable
				'code' to the three-letter abbreviation.
		
			e. makeEmpCalVpop.r
				This script makes the plot where I map the population of data point in the
				empirical bins to the difference between the Theoretical and Empirical maximums
				for each bin. The script needs the 
				"EmpericalVCalculated_diff_pop_nonZeroed_with_pop_restriction_XXX" file and makes 
				"XXX_DifferenceVPopulation" where you have to specify or create the variable 
				'code' in order to get the plot for your amino acid of interest.		

			f. makeNormCorPlot.r
				This script makes the correlation (3x3) plot. This script needs the "Wolfden.txt",
				"Rose.txt", and "Kite_Doolittle.txt". It also needs the 
				"Hydrophobicity_Scales_updated.txt" file. It outputs the 
				'NormalizedCorrelations.pdf.' This file takes the pearson correlation. 
				
			g. makeALLOWEDBinnedRamaPlot.r
				Same idea as II.5.c but with the "XXX_max_bins_ALLOWED"
				
			h. makeCOREBinnedRamaPlot.r
				Same idea as II.5.c but with the "XXX_max_bins_CORE"
				
			i. makeGENEROUSBinnedRamaPlot.r
				Same idea as II.5.c but with the "XXX_max_bins_GENEROUS"
				
			j. getAngles.R
			
				FOR DARIA to fill out
				
				
II. Data Files
	Descriptions how some of the data files look like.

	A. Xxx_geo files (RSA-normalization-values\GeoFiles\
		These files hold all bond angles, bond length, dihedral angles, SA, RSA (Miller), neighboring 
		residues, and secondary structure from the PDB files I mined from the 			
		"cullpdb_pc30_res1.8_R0.25_d130607_chains4961.gz" file. This is tab delimited. I read this into R
		with the read.delim command
		
		Needs: "get_PDB.py", "cullpdb_pc30_res1.8_R0.25_d130607_chains4961.gz", DSSP program, DSSPData.py
		generated from program: parse_alignment

	B. AnglesIteratedThroughAgainXXX (RSA-normalization-values\AnglesIteratedThroughAgain)
		These files are the 1 degree discrete rotations of all phi and psi conformations. I did not record
		the rotamer conformations that gave the high SA. This file has three column: "SA\tPhi\tPsi". It is
		tab delimited.

		Needs: "Geometry.py", "PeptideBuilder.py", DSSP
		generated from program: iterateThroughModels.py

	C. XXX_SA_Over/Under, Xxx_Rose_RSA (RSA-normalization-values\SA_Over_Under\Over, ...\Under)
		These files are pretty self explanitory. The "XXX_SA_Over/Under" files have the phi and psi values
		of all the residues with RSA >1 and RSA <=1, respectively. The file is tab delimited. Xxx_Rose_RSA
		has the normalized SA values using the Rose normalization values. It has two columns that are tab
		delimited. The AA and its neighbors followed by the RSA (with Rose constant)

		Needs: "Xxx_geo"
		generated from program: "SeperateOverAndUnderRSA1.py" or "getRoseRSA.py"

	D. XXX_max_bins_all
		These files are the binned empirical and theoretical results. This file has four columns. 
		"Phi \t Psi \t max_obs_SA \t max_theo_SA \t obs_bin_pop". This is also tab delimited.

		Needs: "Xxx_geo"
		generated from program: "max_bin_all_data.py"

	E. EmpericalVCalculated_diff_pop_nonZeroed_with_pop_restriction_XXX (RSA-normalization-values\EmpiricalVTheoretical)
		These files are used to make the population difference figures. The file has two columns the 
		SA_difference for a bin and the population of the bin. This is also tab delimited. 

		Needs: "XXX_max_bins_all"
		generated from program: "EmpVCalc_get_Diff_with_pop_restriction.py"

	F. NormalizationValuesByPercentDataCoverage
		This is a tab delimited file is the output of the R script getMaximumValues.r. This is commonly
		used in every R script to normalize things.

		Needs: "Xxx_geo" and "AnglesIteratedThroughAgainXXX"
		generated from program: "get_population_cut_offs.r"
	
	G. NormalizationValuesByPercentDataCoverageAndGenerous.txt
		This is a tab delimited file is the output of the R script getMaximumValues.r. This is commonly
		used in every R script to normalize things.

		Needs: "Xxx_ALLOWED_geo" and "AnglesIteratedThroughAgain_ALLOWED_XXX"
		generated from program: "get_population_cut_offs.r"

	H. Hydrophobicity_Scales_Updated.txt
		This is a tab delimited file that contains the Empirical and Theoretical normalized RSA mean, 
		median, square root mean, box-cox mean, and the the percent buried residues. 
		
		Needs: "NormalizationValuesByPercentDataCoverageAndGenerous.txt", "Xxx_geo"
		generated from program: "getMeanSA.r"
	
	I. Wolfden.txt, rose.txt, Kite_Doolittle.txt, Fauchere.txt, Wimley.txt, Moon.txt, Radzicka.txt, MacCallum.txt 
		This is a tab delimited file that contains the hydrophobic values from each of the respective
		papers.

	J. cullpdb_pc30_res1.8_R0.25_d130607_chains4961.gz
		In the list below, the resolution and percent identity cutoffs are given in each filename. E.g., for cullpdb_pc20_res1.8_R0.25_d130517_chains3211, the percentage identity cutoff is 20%, the resolution cutoff is 1.8 angstroms, and the R-factor cutoff is 0.25. The list was generated on May 22, 2013. The number of chains in the list is 3211

	K. Allowed, Core, and Generous Bins
		Files used to make figures and to estimate resonable Ramachandran angle cut-offs
