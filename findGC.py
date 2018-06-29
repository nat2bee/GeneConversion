#!/usr/bin/python

"""
#
# Find gene conversion (of SNP type variants) using a phased file from "PedPhase.py". Analyzes must be performed by chrm. 
# Graphs are generated using R, so it needs to be installed.
#
# Usage: findGC.py -L <ind1,ind2,indN> -p <pedigree_info> -i <PedPhase_final.txt> -o <output_prefix>
#
# Where:
# ind1,ind2,indN = list of IDs for the individuals in which you want to search for gene conversion events (as in the PedPhase_final.txt), i.e. events happened in the parental meiose of these individuals but will be found in them. Minimun of one ID is needed. Individual and its parents need to be phased in PedPhase_final.txt
# pedigree_info = a file containing the pedigree information for all individuals of interest (format file as described in http://zzz.bwh.harvard.edu/plink/data.shtml#ped)
# PedPhase_final.txt = final output generated with "PedPhase.py"
# output_prefix = prefix name to save the output(s)
#
#
# Options: -h for usage help
#          -goff By default the program confirms the gene conversion events searching for the transmision of the haplotype to the individual offspring. If the offspring is not available in the dataset, or if you just desire to remove this feature use this option
#		   -e <event_size> Only consider as possible gene conversion events of size equal to or smaller than this number (i.e. tract size of the gene conversion event). Deafult = 3000bp.
#		   -w <window_size> Number of events to include before and after the start-end of the gene conversion event in the picture file. (E.g. '-w 50' will include 50 marker before the start of the event and 50 after). Deafult = 20.
#
#
# Outputs: One file per individual named <output_prefix>_<id>_all_phase.txt, and a .pdf file per individual containing figures illustrating the region where the gene conversion ocurred. 
#			In both outputs 0 means not informative site (i.e. unphased); 1 means allele received from the paternal chrm of the parent; 2 means allele received from the maternal chrm of the parent
#			Columns in *_all_phase.txt means CHROM	POS 	MALE_GAM_PHASE_INFO		FEMALE_GAM_PHASE_INFO
# 			Where: CHROM and POS, as in PedPhase_final.txt; MALE_GAM_PHASE_INFO and FEMALE_GAM_PHASE_INFO, is the code (as explained above) for the haplotype in each position considering the paternal or maternal gamete received, respectively. (usefull for graphical visualization of the event)
#
"""

###########

## Create default variables

goff_analyses = True
event_size = 3000
window_size = 20
link_info = False
link_info_GO = False

papa = False
mama = False


#############


## Import necessary packages

import re, sys, getopt

import os



#########################################


# Check for the arguments, open the inputs and print useful help messages

try:
    opts, args = getopt.getopt(sys.argv[1:],"he:w:L:p:i:o:",["goff", "link", "linkG"])
except getopt.GetoptError:
    print '\n', '####     Invalid use     ####', '\n'
    print 'Usage: findGC.py -L <ind1,ind2,indN> -p <pedigree_info> -i <PedPhase_final.txt> -o <output_prefix>'
    print 'For help use findGC.py -h'
    sys.exit(99)
    

for opt, arg in opts:
    if opt == '-h':
        print '\n', 'Find gene conversion (of SNP type variants) using a phased file from "PedPhase.py". Analyzes must be performed by chrm. '
        print 'Graphs are generated using R, so it needs to be installed.', '\n'
        print 'Usage: findGC.py -L <ind1,ind2,indN> -p <pedigree_info> -i <PedPhase_final.txt> -o <output_prefix>', '\n'
        print 'Where: ind1,ind2,indN = list of IDs for the individuals in which you want to search for gene conversion events (as in the PedPhase_final.txt), i.e. events happened in the parental meiose of these individuals but will be found in them. Minimun of one ID is needed. Individual and its parents need to be phased in PedPhase_final.txt'
        print 'pedigree_info = a file containing the pedigree information for all individuals of interest (format file as described in http://zzz.bwh.harvard.edu/plink/data.shtml#ped)'
        print 'PedPhase_final.txt = final output generated with "PedPhase.py"'
        print 'output_prefix = prefix name to save the output', '\n'
        print 'Options: -h for help (print this menu)'
        print '--goff By default the program confirms the gene conversion events searching for the transmision of the haplotype to the individual offspring. If the offspring is not available in the dataset, or if you just desire to remove this feature use this option'
        print '--link By default the program will not correct non-informative markers in between informative ones in the probands (e.g. 101 -> 111), but if this opition is used this will be performed before the discovery of events'
        print '--linkG Like --link but to correct non-info markers in the grand offspring if the transmission is analyzed.'
        print '-e <event_size> Only consider as possible gene conversion events of size equal to or smaller than this number (i.e. tract size of the gene conversion event). Deafult = 3000 (base pairs).'
        print '-w <window_size> Number of markers to include before and after the start-end of the gene conversion event in the picture file. (E.g. "-w 50" will include 50 marker before the start of the event and 50 after). Deafult = 20.'
        sys.exit()
        
    elif opt in ("-L"):
        probands = arg
        probands = probands.split(",")
    elif opt in ("-p"):
        ped = arg
    elif opt in ("-i"):
        phase = arg
    elif opt in ("-o"):
        out_prefix = arg
    elif opt in ("--goff"):
        goff_analyses = False
    elif opt in ("--link"):
        link_info = True
    elif opt in ("--linkG"):
        link_info_GO = True
    elif opt in ("-e"):
        event_size = int(arg)
    elif opt in ("-w"):
        window_size = int(arg)

    else:
        assert False, "unhandled option"
        


##########################################

## Define functions

############################################

## Function to find possible gene conversion events

def find_GC(prob, event_size, ped, phase, out_prefix):
	
	## Variables in the function
	father = ""
	mother = ""
	chrm = ""
	prob_phase = "" 
	father_phase = ""
	mother_phase = ""
	prob_phase2 = "" 
	father_phase2 = ""
	mother_phase2 = ""
	position = ""
	
	h = 0
	prob_pos = 0
	prob_Mgam = 0
	prob_Fgam = 0

	father_Mgam = 0
	father_Fgam = 0

	mother_Mgam = 0
	mother_Fgam = 0
	
	pos = 0
	
	start = 0
	end = 0
	n_markers = 0
	start2 = 0
	end2 = 0
	n_markers2 = 0
	
	not_confirmed_father = 0
	not_confirmed_mother = 0
	
	MG = 0
	FG = 0
	
	started = False
	started2 = False
	
	dictGAM = dict()
	dictGAM2 = dict()

	
	
	## get the parents ID
	pedigree = open(ped)
	for entry in pedigree:
		if entry.startswith("#"): # skip the header
			continue
		entry = entry.split("\t")
		id = str(entry[1])
		if id == prob : 
			father = str(entry[2])
			mother = str(entry[3])
		else:
			continue
	
	
	# start output for the graph file
	outname = out_prefix + prob + "_uncor_4graph.txt"
	out1 = open(outname,"w")

	#useful info
	print "Determining gametes received from the parents.\n"
	
	# Compare proband and parents haplotypes for each position and gamete. Outuput code (0,1,2)
	final_file = open(phase)	
	for line in final_file: 
		line = line.split("\n")
		line = line[0]
		line = line.split("\t") 
		
		# get the position of the proband and of the parent in the header of the phase file
		if h == 0: 
			prob_pos = line.index(prob)
			father_pos = line.index(father)
			mother_pos = line.index(mother)
			
			h = 1
			continue

		# Compare each position
		else:
			chrm = line[0]
			pos = line[1]
			prob_phase = line[prob_pos] # To compare with the paternal gametes
			prob_phase2 = line[prob_pos] # To compare with maternal gametes
			
			## Compare the haplotype in each position

			if ":" in prob_phase:  # if the proband is phased
				prob_phase = prob_phase.split(":")
				prob_phase2 = prob_phase2.split(":")

				# compare the parental haplotype with the father
				father_phase = line[father_pos]
				if ":" not in father_phase:  # if the father is not phased, site is not informative
					MG = 0
					dictGAM[pos] = 0

				if ":" in father_phase:
					father_phase = father_phase.split(":")
					father_gam1 = father_phase[0]
					father_gam2 = father_phase[1]
					prob_phase = prob_phase[0]

					
					## Check which gamete received from father

					if (prob_phase == father_gam1) and (prob_phase != father_gam2):
						MG = 1
						dictGAM[pos] = 1
						
					elif (prob_phase == father_gam2) and (prob_phase != father_gam1):
						MG = 2
						dictGAM[pos] = 2
						
					else:
						MG = 0
						dictGAM[pos] = 0
						
					
					
				# compare the parental haplotype with the mother
				mother_phase = line[mother_pos]
				if ":" not in mother_phase:  # if the father is not phased, site is not informative
					FG = 0
					dictGAM2[pos] = 0

				if ":" in mother_phase:
					mother_phase = mother_phase.split(":")
					mother_gam1 = mother_phase[0]
					mother_gam2 = mother_phase[1]
					prob_phase2 = prob_phase2[1]

					
					## Check which gamete received from father

					if (prob_phase2 == mother_gam1) and (prob_phase2 != mother_gam2):
						FG = 1
						dictGAM2[pos] = 1
						
					elif (prob_phase2 == mother_gam2) and (prob_phase2 != mother_gam1):
						FG = 2
						dictGAM2[pos] = 2
						
					else:
						FG = 0
						dictGAM2[pos] = 0
						
			else: # if it is not phased in proband save as not informative site
				MG = 0
				dictGAM[pos] = 0
				FG = 0
				dictGAM2[pos] = 0

				
			# Save the output for graph
			out1.write(str(chrm) + "\t" + str(pos) + "\t" + str(MG) + "\t" + str(FG) + "\n")
	
	
	# Close output
	out1.close()
		
######################	
	keys = sorted(map(int,dictGAM.keys()))
	keys2 = sorted(map(int,dictGAM2.keys()))
	
	## If the correction by flaking markers of the non-informative sites is requested 
	if link_info is True :
		#useful info
		print "Review non-info sites based on flaking markers."
	
		# start output for the graph file
		outname = out_prefix + prob + "_cor_4graph.txt"
		out2 = open(outname,"w")
	
		## Update non-info positions based on flanking markers
		## For the paternal chromossome
		t = 0
		c = 10000
	
		for position in dictGAM:
			MG = dictGAM[position]
			first_key = str(keys[0])
			last_key = str(keys[(len(keys)-1)])
		
			if ((position != first_key) and (position != last_key)) and MG == 0 : # first/ last value in the dictionary, don't correct. If the site is non-info, see if can correct it
				
				# Get the flanking markers
				ind = keys.index(int(position))
				next_key = str(keys[ind+1])
				next_value = dictGAM[next_key]
			
				previous_key = str(keys[ind-1])
				previous_value = dictGAM[previous_key]
				
				t= t+1
				if t == c:
					print "Total markers to check:", len(keys) , "\nChecked:", t
					c = c + 10000

				# Check if next value is info or not, if not keep until one info is found
				if next_value == 0:
					n_ind = ind+1
					while (next_value == 0) : # stop if one info site is found
						if next_key == last_key: # stop if the last position is reached
							break
						n_ind = n_ind+1
						next_key = str(keys[n_ind])
						next_value = dictGAM[next_key]
					
				# Check if previous value is info or not, if not keep until one info is found
				if previous_value == 0:
					p_ind = ind-1
					while (previous_value == 0) :  # stop if one info site is found
						if previous_key == first_key: # stop if the first position is reached
							break
						p_ind = p_ind-1
						previous_key = str(keys[p_ind])
						previous_value = dictGAM[previous_key]

			
				# Compare with flaking markers
				if (previous_value == 1) and (next_value == 1):
					dictGAM[position] = 1

				elif (previous_value == 2) and (next_value == 2):
					dictGAM[position] = 2
	
				else: # not possible to correct if flaking markers are non-info
					dictGAM[position] = 0
			
		print "Paternal chrm corrected."
		
		## For the maternal chromossome
		
		t = 0
		c = 10000
	
		for position in dictGAM2:
			FG = dictGAM2[position]
			first_key = str(keys2[0])
			last_key = str(keys2[(len(keys2)-1)])
		
			if ((position != first_key) and (position != last_key)) and FG == 0 : # first/ last value in the dictionary, don't correct. If the site is non-info, see if can correct it

				# Get the flanking markers
				ind = keys2.index(int(position))

				next_key = str(keys2[ind+1])
				next_value = dictGAM2[next_key]
			
				previous_key = str(keys2[ind-1])
				previous_value = dictGAM2[previous_key]
				
				t= t+1
				if t == c:
					print "Total markers to check:", len(keys2) , "\nChecked:", t
					c = c + 10000
			
			
				# Check if next value is info or not, if not keep until one info is found
				if next_value == 0:
					n_ind = ind+1
					while (next_value == 0) :
						if next_key == last_key:
							break
						n_ind = n_ind+1
						next_key = str(keys2[n_ind])
						next_value = dictGAM2[next_key]
					
				# Check if previous value is info or not, if not keep until one info is found
				if previous_value == 0:
					p_ind = ind-1
					while (previous_value == 0) :
						if previous_key == first_key:
							break
						p_ind = p_ind-1
						previous_key = str(keys2[p_ind])
						previous_value = dictGAM2[previous_key]

			
				# Compare with flaking markers
				if (previous_value == 1) and (next_value == 1):
					dictGAM2[position] = 1
		
				elif (previous_value == 2) and (next_value == 2):
					dictGAM2[position] = 2
	
				else: # not possible to correct if flaking markers are non-info
					dictGAM2[position] = 0
			
	
			
				# Save the output for graph
				out2.write(str(chrm) + "\t" + str(position) + "\t" + str(dictGAM[position]) + "\t" + str(dictGAM2[position]) + "\n")
	
	
		# Close output
		out2.close()
		
		print "Maternal chrm corrected.\n"

#####################	
	## Find the events with change of markers in the dataset
	
	#useful info
	print "Looking for possible events...\n"
	
	## variables
	her1 = 0
	her2 = 0
	
	t_events = 0
	s = 0
	
	start_info = False
	started = False
	
	lst_markers = list()
	
	## Possible gene conversion events in the paternal gamete
	for position in keys:
		her2 = dictGAM[str(position)]
	
		
		if start_info is False: # for the first position
			her1 = her2
			start_info = True
			continue
		
		# No event
		if (her2 != 0) and (her1 != 0) and (her2 == her1) and (started is False):
			continue
			
		elif (her2 != her1) and (her1 == 0) and (started is False):
			her1 = her2
			continue
		
		# If there was a new change in markers	(start of the event)
		elif (her1 != 0) and (her2 != 0) and (her2 != her1) and (started is False): # if had a change and marker is info
			# define start site
			her_test = 0
			i = 1
			while her_test == 0:
				start = keys.index(position) - i
				start = int(keys[start])
				her_test = dictGAM[str(start)]
				i = i+1
			#print position, her1, her2, start
			started = True
			lst_markers.append(position)
			n_markers = n_markers + 1 # start markers counts
			her1 = her2
			continue
			
		
		# Event continues
		elif (her2 == her1) and (her2 != 0) and (started is True):
			n_markers = n_markers + 1 # continue markers counts
			lst_markers.append(position)
				
		
		# Event ends
		elif (her2 != her1) and (her1 != 0) and (her2 != 0) and (started is True):
			#print position, her1, her2, start, (int(position) - start), event_size
			her1 = her2
			end = int(position)
			started = False
			
			
			if (end - start) <= event_size:
				t_events = t_events + 1
				# start output for the unconfirmed GC if needed
				if s == 0:
					outname = out_prefix + prob + "_unconfirmedGC.txt"
					out3 = open(outname,"w")
					out3.write("CHROM\tstart\tend\tn_markers\tpositions\tproband\tparent\n")
					s = 1	
				
				not_confirmed_father = not_confirmed_father + 1 # count the number of events
				print "Possible GC event from paternal chrm in individual", prob, ".\nEvent start:end", start, ":", end, "\nSize=", (end - start), "\tNumber of info markers=", n_markers, "\nNot confirmed by offspring transmition yet"
								
				# Save the output of the unconfirmed GC
				if len(lst_markers) > 1:
					m_positions = ",".join(str(e) for e in lst_markers)
				elif len(lst_markers) <= 1:
					m_positions = str(lst_markers[0])
				out3.write(str(chrm) + "\t" + str(start) + "\t" + str(end) + "\t" + str(n_markers) + "\t" + m_positions + "\t" + str(prob) + "\t" + str(father) + "\n")
			
			# Restart the markers counts		
			n_markers = 0
			lst_markers = list()
			continue
			
			

	
	print "\nAnalyses completed in the paternal chromossome."
	print "Total events found:", t_events, "\n"
	
	## In the maternal gamete
	start_info = False
	n_markers = 0
	started = False
	her1 = 0
	her2 = 0
	t_events = 0
	lst_markers = list()
	
	for position in keys2:
		her2 = dictGAM2[str(position)]
		
		if start_info is False: # for the first position
			her1 = her2
			start_info = True
			continue

		# No event
		if (her2 != 0) and (her1 != 0) and (her2 == her1) and (started is False):
			continue
			
		elif (her2 != her1) and (her1 == 0) and (started is False):
			her1 = her2
			continue

		
		# If there was a new change in markers	(start of the event)
		if (her1 != 0) and (her2 != 0) and (her2 != her1) and (started is False): # if had a change and marker is info
			
			# define start site
			her_test = 0
			i = 1
			while her_test == 0:
				start = keys2.index(position) - i
				start = int(keys2[start])
				her_test = dictGAM2[str(start)]
				i = i+1
			
			started = True
			lst_markers.append(position)
			n_markers = n_markers + 1 # start markers counts
			her1 = her2
			continue
			
		
		# Event continues
		elif (her2 == her1) and (her2 != 0) and (started is True):
			n_markers = n_markers + 1 # continue markers counts
			lst_markers.append(position)
			

		# Event ends
		elif (her2 != her1) and (her1 != 0) and (her2 != 0) and (started is True):
			her1 = her2
			end = int(position)
			started = False
			
			if (end - start) <= event_size:
				t_events = t_events + 1
				# start output for the unconfirmed GC if needed
				if s == 0:
					outname = out_prefix + prob + "_unconfirmedGC.txt"
					out3 = open(outname,"w")
					out3.write("CHROM\tstart\tend\tn_markers\tpositions\tproband\tparent\n")
					s = 1	
					
				not_confirmed_father = not_confirmed_father + 1 # count the number of events
				print "Possible GC event from maternal chrm in individual", prob, ".\nEvent start:end", start, ":", end, "\nSize=", (end - start), "\tNumber of info markers=", n_markers, "\nNot confirmed by offspring transmition yet"
								
				# Save the output of the unconfirmed GC
				if len(lst_markers) > 1:
					m_positions = ",".join(str(e) for e in lst_markers)
				elif len(lst_markers) <= 1:
					m_positions = str(lst_markers[0])
				out3.write(str(chrm) + "\t" + str(start) + "\t" + str(end) + "\t" + str(n_markers) + "\t" + m_positions + "\t" + str(prob) + "\t" + str(mother) + "\n")
			
			# Restart the markers counts		
			n_markers = 0
			lst_markers = list()
			

	# Close output if open
	if s == 1:
		out3.close()
				
	print "\nAnalyses completed in the maternal chromossome."
	print "Total events found:", t_events, "\n"

##########################################################################
			
# function to save the output graph file for the grand-offspring and also creates a dictionary with positions comparison	

def check_transmission(phase_file, offspring, proband, p_sex, p_parental):
	off_comp = dict()
	gam = ""
	off_gam = ""
	proband_gam = ""
	chrm = ""
	
	# Discover if the event happend in the male/ female chrm of the proband
	pedigree = open(ped)
	for entry in pedigree:
		if entry.startswith("#"): # skip the header
			continue
		
		entry = entry.split("\t")
		p = str(entry[1])
		f = str(entry[2])
		m = str(entry[3])
			
		if (f == p_parental) and (p == proband): 
			gam = "male"
		
		if (m == p_parental) and (p == proband): 
			gam = "female"
	
	
	h = 0

	# Compare proband and parents haplotypes for each position and gamete. Outuput code (0,1,2)
	final_file = open(phase_file)	
	for line in final_file: 
		line = line.split("\n")
		line = line[0]
		line = line.split("\t") 
		
		# get the position of the proband and of the parent in the header of the phase file
		if h == 0: 
			off_pos = line.index(offspring)
			proband_pos = line.index(proband)
			
			h = 1
			continue

		# Compare each position and save in a dictionary/ output
		else:
			chrm = line[0]
			pos = line[1]
			off_phase = line[off_pos] 
			proband_phase = line[proband_pos] 
			
			if (":" in proband_phase) and (":" in off_phase): # must be phased in both to compare
				off_phase = off_phase.split(":")
				proband_phase = proband_phase.split(":")
				
				# Get the right gamete of the offspring
				if p_sex == "male":
					off_gam = off_phase[0]
				if p_sex == "female":
					off_gam = off_phase[1]
				
				# Get the right gamete of the proband
				if gam == "male":
					proband_gam = proband_phase[0]
				if gam == "female":
					proband_gam = proband_phase[1]
				
				# Compare the SNPs (0 =non-info; 1= equal to the proband; 2= different from proband)
				
				if off_gam == proband_gam:
					off_comp[pos] = 1
					
				elif off_gam != proband_gam:
					off_comp[pos] = 2
			
			else: # if it is not phase
				off_comp[pos] = 0

	# Check if the chromosome of the event might have been transmited
	if 1 in off_comp.values():
		if link_info_GO is False: # if the link correction is not necessary
			
			# Print output for graph
			
			# start output for the graph file
			outname = out_prefix + offspring + "_gam_" + p_parental + "_4graph.txt"
			out5 = open(outname,"w")
			
			for site in off_comp:
				out5.write(str(chrm) + "\t" + str(site) + "\t" + str(off_comp[site]) + "\n")
			
			# Close output
			out5.close()
			
			# Return dict
			return off_comp
			
		if link_info_GO is True:
		
			#useful info
			print "Review non-info sites based on flaking markers for the offspring:", offspring
			
			# First correct the positions before finishing
			keys = sorted(map(int,off_comp.keys()))
			
			# start output for the graph file
			outname = out_prefix + offspring + "_4graph.txt"
			out5 = open(outname,"w")
			
			t = 0
			c = 10000
			
			# Correct non-info site with flaking markers
			for site in off_comp:
				marker = off_comp[site]
				first_key = str(keys[0])
				last_key = str(keys[(len(keys)-1)])
				
				if ((site != first_key) and (site != last_key)) and marker == 0 : # first/ last value in the dictionary, don't correct. If the site is non-info, see if can correct it

					# Get the flanking markers
					ind = keys.index(int(site))
				
					next_key = str(keys[ind+1])
					next_value = off_comp[next_key]
			
					previous_key = str(keys[ind-1])
					previous_value = off_comp[previous_key]
					
					
					t= t+1
					if t == c:
						print "Total markers to check:", len(keys) , "\nChecked:", t
						c = c + 10000
			
					# Check if next value is info or not, if not keep until one info is found
					if next_value == 0:
						n_ind = ind+1
						while (next_value == 0) :
							if next_key == last_key:
								break
							n_ind = n_ind+1
							next_key = str(keys[n_ind])
							next_value = off_comp[next_key]
					
					# Check if previous value is info or not, if not keep until one info is found
					if previous_value == 0:
						p_ind = ind-1
						while (previous_value == 0) :
							if previous_key == first_key:
								break
							p_ind = p_ind-1
							previous_key = str(keys[p_ind])
							previous_value = off_comp[previous_key]

			
					# Compare with flaking markers
					if (previous_value == 1) and (next_value == 1):
						off_comp[site] = 1
		
					elif (previous_value == 2) and (next_value == 2):
						off_comp[site] = 2
	
					else: # not possible to correct if flaking markers are non-info
						off_comp[site] = 0
			
	
			
				# Save the output for graph
				out5.write(str(chrm) + "\t" + str(site) + "\t" + str(off_comp[site]) + "\n")
	
	
			# Close output
			out5.close()
			
			# Return dict
			return off_comp

	else:
		return False
					

	
#########################################################################

        
############################################


## Create useful variables

cwd = os.getcwd()

lst_results = []
off_her_MPgam = []
off_her_FPgam = []


dictGAM = dict()
MP_dictGC = dict()
FP_dictGC = dict()

prob_new = ""
par = ""

MP_gam = ""
FP_gam = ""

individual = ""

father = ""
mother = ""

sex = ""

unconf_GC = ""





## Define the possible GC events in all probands in the input

for individual in probands:
	print "######################################"
	print "\nStarting analyses for:", individual, "\n" 
	GC = False
	find_GC(individual, event_size, ped, phase, out_prefix)  # Already give the final outputs if dont want to check transmission
	
	if goff_analyses is False: # don't need to analyse the transmission
		# insert code to generate the graph in R
		continue 

	# Check the transmission of the unconfirmed events
	if (goff_analyses is True):
		cwd = os.getcwd()
		mama = False
		papa = False
	
		for file in os.listdir(cwd):
			prob_file = individual + "_unconfirmedGC.txt"
			if file.endswith(prob_file):
				GC = True
				unconf_GC = open(file)

		if GC is False: # No event to check
			continue
			
		elif GC is True:

			# start output to write the confirmed events
			outname = out_prefix + individual + "_confirmedGC.txt"
			out4 = open(outname,"w")
			out4.write("CHROM\tstart\tend\tn_markers\tposition\toffspring\tproband\tgamete\n")
			
			
			## Check if the parental Proband is male or female and all its offspring
			pedigree = open(ped)
			off_list = list()

			
			for entry in pedigree:
				if entry.startswith("#"): # skip the header
					continue
		
				entry = entry.split("\t")
				id = str(entry[1])
				f = str(entry[2])
				m = str(entry[3])
		
				if (f == individual): 
					sex = "male"
					off_list.append(id)
		
				if (m == individual): 
					sex = "female"
					off_list.append(id)
			
			
			h = 0
			
			for event in unconf_GC:
				
				if h == 0: # skip header
					h = 1
					continue

				## get info
				event = event.split("\t")
				par = event[6]
				par = par.split("\n")
				par = par[0]
				chrm = event[0]
				start = int(event[1])
				end = int(event[2])
				markers = event[3]
				interest_marker = event[4]
				interest_marker = interest_marker.split(",")

				#useful info
				print "\nChecking transmission of the event from", start, "to", end, "from", par, "meiosis"
				
				# Create a dictionary comparing each position
				for off in off_list:
					off_dict = check_transmission(phase, off, individual, sex, par)
										
					if off_dict is not False:
						
						for pos in off_dict:
	
							if pos in interest_marker: # i.e. marker in inside the event
								comp = off_dict[str(pos)]
								comp_start = off_dict[str(start)]
								comp_end = off_dict[str(end)]
								
								#print start, end, off
								print comp_start, comp_end, comp
							
								if (comp != 2) and (comp_start != 2) and (comp_end != 2): # event marker is confirmed (1) or non-informative (0)
									print "Event marker possibly transmited to", off
									# write in the confirmed output
									out4.write(str(chrm) + "\t" + str(start) + "\t" + str(end) + "\t" + str(markers) + "\t" + str(pos) + "\t" + str(off) + "\t" + str(individual) + "\t" + str(par) + "\n")



			out4.close ()






"""

	



## Print graphic output of the found events

for individual in probands:

	if (goff_analyses is False):
		outname = out_prefix + individual + "_unconfirmedGC.txt"
		
	# print the unconfirmed events (w_size - antes/ depois)

	# use
	
	
					for file in os.listdir(cwd):
						if file.endswith("_unconfirmedGC.txt"):
							unconf_GC = open(file)
	
	outname = out_prefix + off + "_4graph.txt"
	
	for individual in probands:
	
	# options:
	w_size; "_unconfirmedGC.txt",  out_prefix
	
	# save one pdf/ individual; one graph per unconfirmed

	
	
	
	if (goff_analyses is True):
		outname = out_prefix + individual + "_confirmedGC.txt"
		
	# print the confirmed events (w_size - antes/ depois)

	# use
	
					
					for file in os.listdir(cwd):
						if file.endswith("_confirmedGC.txt"):
							conf_GC = open(file)
	
	outname = out_prefix + off + "_4graph.txt"
	
	for individual in probands
		get offs with the same event (same gam, same start)

	
	# options:
	w_size; "_confirmedGC.txt", out_prefix
	
	# save one pdf/ individual; one graph per event with ind. and all offs
	
				
"""



        