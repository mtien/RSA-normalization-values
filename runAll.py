#!usr/bin/python

import os

AA=["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu",
    "Lys", "Met", "Pro", "Phe", "Ser", "Thr", "Trp", "Tyr", "Val"]

##script_line="python getRoseRSA.py"
##os.system(script_line)

for aa in AA:
    ##script_line="python SeperateOverAndUnderRSA1.py " + aa
    ##os.system(script_line)

    ##script_line_max="python max_bins_with_population_restriction.py " + aa
    ##os.system(script_line_max)

    ##script_line_max="python max_bins_with_population_restriction_theoretical.py " + aa
    ##os.system(script_line_max)

    ##script_line_difference="python EmpVCalc_get_Diff_with_pop_restriction.py " + aa
    ##os.system(script_line_difference)

    script_line_max="python max_bin_all_data.py " + aa
    os.system(script_line_max)

    print aa
