RSA normalization Project

----------------------------------------------------------------------------------------------------------------

Table of Contents:

I. Important Programs
	A. Mine PDB information
		1. parse_alignment.py		
		2. get_PDB.py
	B. Theoretical Model Construction
		1. makeStructure.py
		2. Geometry.py
		3. PeptideBuilder.py
		4. makeFigure.py
	C. Theoretical data generation
		1. iterateThroughModels.py
		2. iterateThroughPhiPsi.py
	D. Data Analysis
		1. max_bins_with_population_restriction.py
		2. SeperateOverAndUnderRSA1.py
		3. max_bins_with_population_theoretical.py
		4. EmpVCalc_get_Diff_with_pop_restriction.py
		5. getRoseRSA.py
	E. R scripts
		1. getMaximumValues.r
		2. getMeanSA.r
		3. CorrelatonTable.r
		4. FigureScripts
			a. barGRSA.r
			b. getRSAdistribution.r
			c. makePlotsWithPopRestriction.r
			d. makeRamaPlot.r
			e. makeEmpCalVpop.r
			f. makeNormCorPlot.r
II. Data Files
	A. Xxx_geo files
	B. AnglesIteratedThroughAgainXXX
	C. XXX_SA_Over/Under, Xxx_Rose_RSA
	D. XXX_max_emperical_bins_pop_restriction
	E. XXX_max_theoretical_bins_Again
	F. EmpericalVCalculated_diff_pop_nonZeroed_with_pop_restriction_XXX
	G. NormalizationValues.txt
	H. Hydrophobicity_Scales_Updated.txt
	I. Wolfden.txt, rose.txt, Kite_Doolittle.txt
	J. cullpdb_pc20_res1.8_R0.25_d130517_chains3211

----------------------------------------------------------------------------------------------------------------

I. Important Programs
	This section is dedicated to all python, R, and other scripts used to obtain or analyze data.
	
	A. Mine PDB Information
		To mine the protein structures I parsed the "cullpdb_pc20_res1.8_R0.25_d110520_chains2908" 
		file in "parse_alignment.py".
		
		1. parse_alignment.py
			This is a python script that takes in a list of PDB and Chain ID's (I used the 
			"cullpdb_pc20_res1.8_R0.25_d110520_chains2908" file from the Dunbrak lab).
			This program creates the "Xxx_geo" files (ie. Ala_geo, Asn_geo, etc...).
			This program also only outputs information from non-chain-terminating residues, which
			is determined by the peptide bond length between two residues.
			This program is mainly a frame to process the output from the real PDB parser "get_PDB"
			Any files that is either corrupted or did not exist in the PDB database is reported to
			the "Error_report_bond_length.txt".
			The two numbers that are printed to the terminal are for testing.
	
		2. get_PDB.py
			This is a python program that takes in a PDB file name and Chain ID name, downloads
			the PDB file to the "structures" folder than extracts all the information from the
			PDB file. This program needs the PDB parser from Biopython and DSSP and the proper
			parser for it. The program outputs lists of information. This program contains many
			functions inorder to properly mine the data, most of which I barely understand. needs
			DSSP to work. I made DSSP a command,s you will have to modify this.

	B. Theoretical Model Construction
		I approached Model Construction through a couple of different ways. "makeStructure" is quite
		limited in that it can only make Tripeptides, however, it is useful in make figures in which
		the tripeptide is at the origin. "Geometry" and "PeptideBuilder" are both used to build any
		structure. "makeFigure" is an example script that I used to make figures and use the libraries.
		
		1. makeStructure.py
			This program is not flexible, meaning changing protein parameters is a bit harder to
			do. It makes tripeptide in a unique way, however. It put's the amino acid at interest
			at the origin and then builds the glycine atoms to the right and left of the origin.
			The main function is the calculate coordinates method and the other methods simply
			build off this method. There is a construction method for each amino acid.

		2. Geometry.py
			This is a library of Amino Acid Geometry objects. By reading in the one-letter Amino
			Acid abbreviations, it can create the correct geometric parameters to construct a
			protein residue. There are 20 classes (one for each amino acid) and one function
			to create the geometry object. Some of these classes have an inputRotamers function
			that takes in a list of integers in order to change the Amino Acid's rotamers. The
			generateRandomRotamers method is a nice method that I used in "iterateThroughModels"

		3. PeptideBuilder.py
			This takes in geometry objects to construct amino acid chains. There are 20 methods
			to construct the amino acids and has the calculate coordinates method in "makeStructure"
			The program also contains a makeStructure method that takes in a string of Amino Acids,
			Phi list of float, Psi list of float numbers. It contains two add residue methods, an
			initialize_residue method, and a makeExtended Structure method. It also has a output
			structure method. This program creates pdb files of your name choice

		4. makeFigure.py
			This was a script I used to make figures. It's a good example of how to use 
			PeptideBuilder

	C. Theoretical data generation
		After each model was built, I needed to iterate through the models phi and psi conformations.
		"iterateThroughModels" is a script that uses Geometry and Peptide Builder to build the
		phi and psi conformation, and if possible, the rotamer conformations. "iterateThroughPhiPsi"
		is less sophisticated in that it can only iterate through the phi and psi conformations. It
		uses, however, a nice mathematical rotation matrix to iterate through the phi and psi 
		conformations.
	
		1. iterateThroughModels.py
			The program need the "Geometry" and "PeptideBuilder" to create the phi and psi
			conformation. This has a nice DSSP method and methods to iterate through the psi, psi,
			and chi angles of a residue. This program creates "AnglesIteratedThroughAgainXXX"
			
		2. iterateThroughPhiPsi.py	
			The program has some nice testing output, in that you can print out the protein
			geometries of the tripeptide to ensure that everything loooks right. It also has 
			a cool rotation matrix operation that rotates the neighboring residues. 			
			This program creates "AnglesIteratedThroughXXX"

	D. Data Analysis
		After obtaining the information in the "Xxx_geo" files and rotating through the theoretical
		models, I wrote a bunch of scripts to help me analyze the data.
		
		1. max_bins_with_population_restriction.py
			This program takes in the "Xxx_geo" data file and bins the data into 5-degree by 5-
			degree Phi and Psi coordinates and put in the max SA found for that bin and the number
			of data points in that bin. The input of the program requires an all-caps three letter
			abbreviation of which amino acid you want to look at. This outputs the "XXX_max_
			emperical_bins_pop_restriction".
	
		2. SeperateOverAndUnderRSA1.py
			This program takes in the "Xxx_geo" data file and creates two files "XXX_SA_Over/Under"
			This program is fairly basic, but it's a good basis to parse the "Xxx_geo" files. The
			command line argument is the all-caps three letter abbreviation of which amino acid you
			want to look at.

		3. max_bins_with_population_theoretical.py
			This program takes in the "XXX_max_emperical_bins_pop_restriction" file and the 
			"AnglesIteratedThroughAgainXXX" file to make the "XXX_max_theoretical_bin_Again" file.
			The command line argument is the all-caps three letter abbreviation of which amino acid 
			you want to look at. This file makes sure that the theoretical data is binned and
			treated in parallel with the empirical data.

		4. EmpVCalc_get_Diff_with_pop_restriction.py
			This program takes in the "XXX_geo" file and the "AnglesIteratedThroughAgain" file to
			create the "EmpericalVCalculated_diff_pop_nonZeroed_with_pop_restriction_XXX" file.
			The command line argument is the all-caps three letter abbreviation of which amino acid 
			you want to look at. This program just compares the bins from both "XXX_max_bins" files.

		5. getRoseRSA.py
			This program main use is to make the "Xxx_Rose_RSA" data files, which contain the
			normalized RSA values using the normalization constants from the Rose paper. The 
			input for this program is the "Xxx_geo" files.
			
	E. R scripts
		The R-scripts were mainly used to make figures and to do simple things that would have been a 
		bit more complicated to do in python.

		1. getMaximumValues.r
			This gets the maximum SA values from both "XXX_max_bins" data files. It outputs the
			"NormalizationValues.txt" in a csv format.

		2. getMeanSA.r
			This gets the mean RSA, median RSA, square root mean RSA, box-cox transformed mean RSA,
			fraction of 100% buried residues, and fraction of 95% residues (for theoretical 
			normalization values only). For all the average estimates, the script uses both 
			empirical and theoretical normalization values. The output file is a csv file called
			"Hydrophobicity_Scales_updated.txt". It also has optional "MeanHydrophobicityScales.txt"
			and "BuriedHydrophobicityScales.txt" output. The script requires the
			"NormalizationValues.txt" file.

		3. CorrelatonTable.r
			This makes the correlation table of all scales in "Hydrophobicity_Scales_updated.txt"
			with the Wolfenden, Kyte, and Rose scales. This require all four files to run. This
			performs the pearson correlation test. Needs the "Wolfden.txt", "Rose.txt", and
			"Kite_Doolittle.txt" files. 

		4. Figure Scripts
			These scripts were used to make the figures in R. most have a pdf/png alternate code in 
			their scripts that are commented out.

			a. barGRSA.r
				Using the "Xxx_geo" files and the "Xxx_Rose_RSA" files, the script makes the 
				"BarGraphRSA.pdf" figure. There is an optional png script at the bottom of this
				file, it is commented out currently.
 
			b. getRSAdistribution.r
				This file can use somework. I only made the script to make the RSA distribution
				Alanine. However, with a bit of tweaking it can make all Amino Acids. There is 
				a script however (III.D) that can make all the RSA distributions for each amino
				acid. This makes the "Ala_dist.pdf" This uses the hist function in R.

			c. makePlotsWithPopRestriction.r
				This file makes the best figure. It makes the "XXX_Rama_HSV.svg" figure (the 
				figure that compares the Theoretical and Empirical results. This figure needs
				to be edited in Inkscape. It requires both "XXX_max_bins" files. The file requires
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
II. Data Files
	Here I jsut describe how some of the data files look like.

	A. Xxx_geo files
		These files hold all bond angles, bond length, dihedral angles, SA, RSA (Miller), neighboring 
		residues, and secondary structure from the PDB files I mined from the 			
		"cullpdb_pc20_res1.8_R0.25_d110520_chains2908" file. This is tab delimited. I read this into R
		with the read.delim command
		
		Needs: "get_PDB.py", "cullpdb_pc20_res1.8_R0.25_d110520_chains2908", DSSP 
		generated from program: parse_alignment

	B. AnglesIteratedThroughAgainXXX
		These files are the 1 degree discrete rotations of all phi and psi conformations. I did not record
		the rotamer conformations that gave the high SA. This file has three column: "SA\tPhi\tPsi". It is
		tab delimited.

		Needs: "Geometry.py", "PeptideBuilder.py", DSSP
		generated from program: iterateThroughModels.py

	C. XXX_SA_Over/Under, Xxx_Rose_RSA
		These files are pretty self explanitory. The "XXX_SA_Over/Under" files have the phi and psi values
		of all the residues with RSA >1 and RSA <=1, respectively. The file is tab delimited. Xxx_Rose_RSA
		has the normalized SA values using the Rose normalization values. It has two columns that are tab
		delimited. The AA and its neighbors followed by the RSA (with Rose constant)

		Needs: "Xxx_geo"
		generated from program: "SeperateOverAndUnderRSA1.py" or "getRoseRSA.py"

	D. XXX_max_emperical_bins_pop_restriction
		These files are the binned empirical results. I know empirical is spelled wrong. Sorry. This file
		has four columns. "Phi \t Psi \t maxSA \t population". This is also tab delimited.

		Needs: "Xxx_geo"
		generated from program: "max_bins_with_population_restriction.py"

	E. XXX_max_theoretical_bins_Again
		These files are the binned theoretical results that mirroed the binning and population 
		restrictions set in the "XXX_max_empirical" file. The file is tab delimited with "Phi \t
		Psi \t maxSA"

		Needs: "Xxx_geo", "AnglesIteratedThroughAgain"
		generated from program: "max_bins_with_population_theoretical.py"

	F. EmpericalVCalculated_diff_pop_nonZeroed_with_pop_restriction_XXX
		These files are used to make the population difference figures. The file has two columns the 
		SA_difference for a bin and the population of the bin. This is also tab delimited. 

		Needs: "XXX_max_emperical_bins_pop_restriction" and "XXX_max_theoretical_bins_Again"
		generated from program: "max_bins_with_population_theoretical.py"

	G. NormalizationValues.txt
		This is a csv file that reads is the output of the R script getMaximumValues.r. This is commonly
		used in every R script to normalize things. I use the read.csv("NormalizationValues.txt", 
		row.names=1) command in R to read this into a data structure. The header is "names, Theoretical,
		Empirical"

		Needs: "Xxx_geo" and "AnglesIteratedThroughAgainXXX"
		generated from program: "getMaximumValues.r"

	H. Hydrophobicity_Scales_Updated.txt
		This is a csv file that contains the Empirical and Theoretical normalized RSA mean, median, square
		root mean, box-cox mean, and the the percent buried residues. This is read into R as 
		read.csv("Hydrophobicity_Scales_Updated.txt", row.names=1)

		Needs: "NormalizationValues.txt", "Xxx_geo"
		generated from program: "getMeanSA.r"
	
	I. Wolfden.txt, rose.txt, Kite_Doolittle.txt
		This is a tab delimited file that contains the hydrophobic values from each of the respective
		papers.

	J. cullpdb_pc20_res1.8_R0.25_d130517_chains3211
		In the list below, the resolution and percent identity cutoffs are given in each filename. E.g., for cullpdb_pc20_res1.8_R0.25_d130517_chains3211, the percentage identity cutoff is 20%, the resolution cutoff is 1.8 angstroms, and the R-factor cutoff is 0.25. The list was generated on May 22, 2013. The number of chains in the list is 3211

