#!/usr/bin/python

"""
#
# Keep only single mutation SNPs, Remove all deletions (*) and SNPs with >1 ALT from a variants table (converted from a vcf using GATK).
#
# Usage: clean_variants.py -i <genotype_table> -o <output_name>
#
# Where:
# genotype_table = table generated by the VariantsToTable tool from GATK (https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php)
# output_name = name to save the resulting file
#
# Options: -h for usage help
#
"""


import re, sys, getopt


# Check for the arguments, open the inputs and print useful help messages

try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:")
except getopt.GetoptError:
    print '\n', '####     Invalid use     ####', '\n'
    print 'Usage: clean_variants.py -i <genotype_table> -o <output_name>'
    print 'For help use clean_variants.py -h'
    sys.exit(99)
    

for opt, arg in opts:
    if opt == '-h':
        print '\n', 'Remove all deletions (*) and SNPs with >1 ALT from a variants table (converted from a vcf using GATK). ', '\n'
        print 'Usage: clean_variants.py -i <genotype_table> -o <output_name>', '\n'
        print 'Where: genotype_table = table generated by the VariantsToTable tool from GATK (https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php)'
        print 'output_name = name to save the resulting file'
        print 'Options: -h for help'
        
        sys.exit()
        
    elif len(arg) >= 2:
        if opt in ("-i"):
            raw_variants = open(arg)
        if opt in ("-o"):
            output = open(arg,"w")


    else:
        assert False, "unhandled option"
        
        
        
### Create useful variables
h = 0
alt_pos = 0

lst_ind_pos = []

ind_gen = ""
a1 = ""

already = False
duplicated = False

### Cleaning the genotype file

for line in raw_variants:
	# Get position of the individuals GT and ALT
	if h == 0:
		header = line.split("\t") 
		alt_pos = header.index("ALT")
		ref_pos = header.index("REF")
		
		for col in header:
			if col.endswith(".GT"):
				ind_pos = header.index(col)
				lst_ind_pos.append(ind_pos)
		
		output.write(line) 
		h = 1
		

	# Check if there is a deletion in line
	elif "*" not in line:
		leave = False
		
		info = line.split("\n")
		info = info[0]
		info = info.split("\t") 
		alt = info[alt_pos]
		alt_len = alt.split(",")
		ref = info[ref_pos]
		
		# Check if all possible alterations are SNPs
		if (len(ref) != 1):
			continue
		
		for variants in alt_len:
			alt_n = variants
			if (len(alt_n) != 1) :
				leave = True
				continue
		
		if leave is True:
			continue
		
		# Check if there are more than 2 ALT
		elif len(alt) == 1 : # if it is only one ALT and it's a single mutation
			output.write(line) 
			continue
		else:
			alt = alt.split(",")			
			duplicated = False
			already = False
			
			for pos in lst_ind_pos: # For all individuals in the file (for each line)
				ind_gen = info[pos]
				for i in alt: # For all possible ALT markers
					if (i in ind_gen) and (already is False) :
						already = True
						a1 = i
					elif (i in ind_gen) and (already is True) and ( i != a1) :
						duplicated = True # mark if more than two ALT are found
						
			if duplicated is False: # if only one alternative is found 
				info[alt_pos] = a1
				new_line = "\t".join(info)
				output.write(new_line) 
					

output.close()	 
	
	
	
	
	
	
	
