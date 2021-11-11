#!/usr/bin/env python

from itertools import combinations
from HaploFunct import *
from GFF_lib import *
from FASTA_lib import *


gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


############################################################################################################################################################
# AGP file format, tab delimited, 1-base coordinates #######################################################################################################
############################################################################################################################################################
#	##agp-version	2
#	# 0                1          2          3        4            5                   6                    7                  8
#	# Obj_Name	       Obj_start  Obj_End	 PartNum  Compnt_Type  CompntId/GapLength  CompntStart_GapType  CompntEnd_Linkage  Orientation_LinkageEvidence
#	Super-Scaffold_9   1	      2          1	      N  	       2	               scaffold	            yes	               map
#	Super-Scaffold_9   3          2581380    2	      W	           SMEL_DRAFT_562	   1	                2581378	           -
#	Super-Scaffold_9   2581381    2716079    3	      N	           134699	           scaffold	            yes	               map
############################################################################################################################################################
#	0 - Obj_Name	- Object Out_sequence name (options.out + "_{1,2}_" + Target_sequence)
#	1 - Obj_start 	- Object start on Out_sequence (start form 1, = 1 + prev_Obj_End)
#	2 - Obj_End		- Object end on Out_sequence ( start + CompntEnd - CompntStart | start + GapLength - 1)
#	3 - PartNum		- Progressive element count (start 1)
#	4 - Compnt_Type 				- (W -> sequence) 						| (U -> Gap unknown size) [not use (N -> Gap known size) ]
#	5 - CompntId/GapLength			- query_seq name 						| gap_length
#	6 - CompntStart_GapType			- query_seq start [ == 1]				| gap_type [== scaffold]
#	7 - CompntEnd_Linkage			- query_seq end [ == query_seq(len)]	| Linkage [== yes]
#	8 - Orientation_LinkageEvidence	- strand								| Kind of evidence linking [align_genus]
############################################################################################################################################################

### Query list element format
# Gap:	[61252	,	(0:61252)		,	gap		,	61252	,	0]
#		[length	,	(T_start:Tstop)	, 	"gap" 	, 	length 	, 	0]
# Object:	[b40-14.HS_iter4_seq3733|+	,	(61252:6463804)	,	93612:7595148	,	-6402552			,	4526208]
#			[ID|strand					,	(T_start:Tstop)	,	Q_start:Q_stop	,	-(alignment length)	,	matches]


def read_agp( agp_file ) :
	agp_db = {}
	# agp_db[chr][start]=agp_line.split("\t")
	for agp_line in open(agp_file) :
		if agp_line[0] == "#" : continue
		el = agp_line.rstrip().split("\t")
		if el[0] not in agp_db:
			agp_db[el[0]] = {}
		agp_db[el[0]][int(el[1])] = el

	return agp_db


def write_agp( agp_db , filename ) :
	outfile = open(filename , 'w')

	for seq_id in sorted(agp_db.keys()) :
		for start in sorted(agp_db[seq_id].keys()) :
			print >> outfile , "\t".join( [ str(x) for x in agp_db[seq_id][start] ] )

	outfile.close()

	return filename


def regions_from_agp( agp_db , direction = "new" ) :
	regions = {}

	if direction == "new" :
		for chr in sorted(agp_db.keys()) :
			regions[chr] = {}
			for start in sorted(agp_db[chr].keys()) :
				Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = agp_db[chr][start]
				regions[chr][ ( int(Obj_start) , int(Obj_End) ) ] = { "type" : Compnt_Type , "desc" : [ CompntId , CompntStart , CompntEnd ,  Orientation ] }
	elif direction == "old" :
		for chr in sorted(agp_db.keys()) :
			for start in sorted(agp_db[chr].keys()) :
				Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = agp_db[chr][start]
				if Compnt_Type == "W" :
					if CompntId not in regions :
						regions[CompntId] = {}
					regions[CompntId][ ( int(CompntStart) , int(CompntEnd) ) ] = { "type" : Compnt_Type , "desc" : [ Obj_Name , Obj_start , Obj_End , Orientation ] }
	else :
		print >> sys.stderr, "[ERROR] Unkown database reference base"
		sys.exit(1)
	return regions


def invert_agp( agp_db , test = False , old_seq_length_db = "" ) :
	# given (agp_db: new_id <- old_id) returns (rev_agp_db: old_id <- new_id)
	# NOTE (1): the code needs that the AGP has no overlapping regions in the old annotation.
	# NOTE (2): Does not handles the breaking of overlapping regions, but can test.
	#  A test of uniqueness may be run, but if it fails rises an error, prints the offending lines and kills the main process (sys.exit(2))
	# NOTE (3): If "old_seq_length_db" is not given, rev_agp_db may not cover the original sequence to the end

	rev_agp_db = {}

	if test :
		print >> sys.stderr, "- Testing AGP ranges uniqueness"
		test_results = test_range_uniqueness(agp_db)
		if not test_results == {} :
			# Test failed
			print >> sys.stdout, "-- [ERROR]" + str(len(test_results.keys())) + " original sequences show ranges used multiple times. See error.overlapping_ranges.info for more info"
			print >> sys.stdout , "------------------------------"
			print >> sys.stdout , "- Quitting"
			print >> sys.stdout , "------------------------------"
			print_overlapping_info(test_results, "error.overlapping_ranges.info")
			sys.exit(2)

	# generate new components for each old CompntId
	component_db = {}

	for seq_id in agp_db :
		try :
			elements = sorted( agp_db[seq_id].keys() )
		except :
			print >> sys.stderr, agp_db[seq_id]
			exit(3)
		for element_start in elements :
			try :
				Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = agp_db[seq_id][element_start]
			except :
				print >> sys.stderr, seq_id
				print >> sys.stderr, element_start
				print >> sys.stderr, agp_db[seq_id][element_start]
				exit(1)
			if Compnt_Type == "W" :
				# Translate only sequence components, not gaps
				if CompntId not in component_db :
					component_db[CompntId] = {}
				component_db[CompntId][int(CompntStart)] = [ CompntId , CompntStart , CompntEnd , Obj_Name , Obj_start , Obj_End ,  Orientation ]

	# regenerate agp structure from elements list and patch up the missing regions
	for seq_id in component_db :
		rev_agp_db[seq_id] = {}
		PartNum = 0
		actual_CompntEnd = 0
		for start in sorted(component_db[seq_id].keys()) :
			CompntId , CompntStart , CompntEnd , Obj_Name , Obj_start , Obj_End ,  Orientation = component_db[seq_id][start]
			if not int(CompntStart) == actual_CompntEnd + 1 :
				# Hit do not cover the region directly after the previous one (or the start)
				# Generate a filler component hitting itself
				PartNum += 1
				fill_start = actual_CompntEnd + 1
				fill_end = int(CompntStart) - 1
				rev_agp_db[seq_id][PartNum] = [ seq_id , str(fill_start) , str(fill_end) , str(PartNum) , "W" , seq_id , str(fill_start) , str(fill_end) ,  "+" ]
				print >> sys.stderr, "[invert_agp:NOTE] Unused region from original sequences: " + seq_id + ":" + str(fill_start) + "-" + str(fill_end)

			# make new component out of the element
			PartNum += 1
			rev_agp_db[seq_id][PartNum] = [ seq_id , CompntStart , CompntEnd , str(PartNum) , "W" , Obj_Name , Obj_start , Obj_End , Orientation ]
			actual_CompntEnd = int(CompntEnd)

		if not old_seq_length_db == "" :
			# Check if the last component covers the end of the sequence
			if not actual_CompntEnd == old_seq_length_db[seq_id] :
				# not covered until the end, add a component to patch the region
				PartNum += 1
				fill_start = actual_CompntEnd + 1
				fill_end = old_seq_length_db[seq_id]
				rev_agp_db[seq_id][PartNum] = [ seq_id , str(fill_start) , str(fill_end) , str(PartNum) , "W" , seq_id , str(fill_start) , str(fill_end) ,  "+" ]
				print >> sys.stderr, "[invert_agp:NOTE] Unused region from original sequences: " + seq_id + ":" + str(fill_start) + "-" + str(fill_end)

	return rev_agp_db


def agp2fasta( agp_db , fasta_in_db , direction = "old_to_new" ) :
	fasta_out_db = {}

	if direction == "old_to_new" :
		# fasta_in_db contains original seq
		# fasta_out_db delivers new seq
		for key in agp_db :
			print >> sys.stderr, "## Creating FASTA sequence for " + key
			fasta_out_db[key] = ""
			for start in sorted(agp_db[key].keys()) :
				if agp_db[key][start][4] == "W" :
					old_seq_id = agp_db[key][start][5]
					old_seq_start = int(agp_db[key][start][6]) - 1
					old_seq_stop = int(agp_db[key][start][7])
					my_seq = fasta_in_db[old_seq_id][old_seq_start:old_seq_stop]
					print >> sys.stderr, "## Adding sequence " + old_seq_id + ":" + str(old_seq_start) + "-" + str(old_seq_stop) + " | Region length: " + str(old_seq_stop - old_seq_start +1) + " | Added bases: " + str(len(my_seq))
					if agp_db[key][start][8] == "-" :
						try :
							new_seq = str(Seq(my_seq).reverse_complement())
						except :
							print >> sys.stdout, str(my_seq)
							sys.exit(1)
					else :
						new_seq = my_seq
				else :
					new_seq = "N" * int(agp_db[key][start][5])
					print >> sys.stderr, "## Adding gap of " + str(len(new_seq)) + "bp"

				fasta_out_db[key] += new_seq

	elif direction == "new_to_old" :
		# fasta_in_db contains new seq
		# fasta_out_db delivers original seq
		# Split new into original
		for key in agp_db :
			myseq = fasta_in_db[key]
			for start in sorted(agp_db[key].keys()) :
				if agp_db[key][start][4] == "W" :
					old_seq_id = agp_db[key][start][5]
					seq_start = int(agp_db[key][start][1]) - 1
					seq_stop = int(agp_db[key][start][2])
					old_seq = myseq[seq_start:seq_stop]
					if agp_db[key][start][8] == "-" :
						old_seq = str(Seq(old_seq).reverse_complement())

					fasta_out_db[old_seq_id] = old_seq

	return fasta_out_db


def translate_from_AGP_whole_genome(agp_dict) :
	translation_db = {}
	# translation_db[old_seq_id][(old_seq_start , old_seq_stop)] = [new_seq_id , coords_offset , direction ]
	# if direction == "+" -> offset = new_start - old_start
	# if direction == "-" -> offset = new_stop + old_start

	for original_seq in agp_dict.keys() :
		for start in agp_dict[original_seq].keys() :
			element = agp_dict[original_seq][start]
			if element[4] == "W" :
				# sequence entry, not gap
				new_seq_id = element[0]
				old_seq_id = element[5]
				old_seq_start = int(element[6])
				old_seq_stop = int(element[7])
				direction = element[8]
				if direction == "+" :
					offset = int(element[1]) - int(element[6])
				else :
					#direction == "-"
					offset = int(element[2]) + int(element[6])
				if old_seq_id not in translation_db:
					translation_db[old_seq_id] = {}
				translation_db[old_seq_id][(old_seq_start , old_seq_stop)] = [new_seq_id , offset , direction ]

	return translation_db


#def translate_from_AGP_single_seq(agp_dict) :
#	translation_db = {}
#	# translation_db[old_seq_id][(old_seq_start , old_seq_stop)] = [new_seq_id , coords_offset , direction ]
#	# if direction == "+" -> offset = new_start - old_start
#	# if direction == "-" -> offset = new_stop + old_start
#
#	for start in agp_dict.keys() :
#		element = agp_dict[start]
#		if element[4] == "W" :
#			# sequence entry, not gap
#			new_seq_id = element[0]
#			old_seq_id = element[5]
#			old_seq_start = int(element[6])
#			old_seq_stop = int(element[7])
#			direction = element[8]
#			if direction == "+" :
#				offset = int(element[1]) - int(element[6])
#			else :
#				#direction == "-"
#				offset = int(element[2]) + int(element[6])
#			if old_seq_id not in translation_db :
#				translation_db[old_seq_id] = {}
#			translation_db[old_seq_id][(old_seq_start , old_seq_stop)] = [new_seq_id , offset , direction ]
#
#	return translation_db


def translate_from_AGP_whole_genome_reverse( agp_dict ) :
	translation_db = {}
	# translation_db[new_seq_id][(new_seq_start , new_seq_stop)] = [old_seq_id , coords_offset , direction ]
	# if direction == "+" -> offset = new_start - old_start
	# if direction == "-" -> offset = new_stop + old_start

	#TODO: handle duplicated ids -> uniquify names

	for original_seq in agp_dict.keys() :
		for start in agp_dict[original_seq].keys() :
			element = agp_dict[original_seq][start]
			if element[4] == "W" :
				# sequence entry, not gap
				new_seq_id = element[0]
				new_seq_start = int(element[1])
				new_seq_stop = int(element[2])
				old_seq_id = element[5]
				direction = element[8]
				if direction == "+" :
					offset = int(element[6]) - int(element[1])
				else :
					#direction == "-"
					offset = int(element[2]) + int(element[6])
				if new_seq_id not in translation_db :
					translation_db[new_seq_id] = {}
				translation_db[new_seq_id][(new_seq_start , new_seq_stop)] = [ old_seq_id , offset , direction ]

	return translation_db


#def translate_from_AGP_single_seq_reverse( agp_dict ) :
#	translation_db = {}
#	# translation_db[new_seq_id][(new_seq_start , new_seq_stop)] = [old_seq_id , coords_offset , direction ]
#	# if direction == "+" -> offset = new_start - old_start
#	# if direction == "-" -> offset = new_stop + old_start
#
#	for start in agp_dict.keys() :
#		element = agp_dict[start]
#		if element[4] == "W" :
#			# sequence entry, not gap
#			new_seq_id = element[0]
#			new_seq_start = int(element[1])
#			new_seq_stop = int(element[2])
#			old_seq_id = element[5]
#			direction = element[8]
#			if direction == "+" :
#				offset = int(element[6]) - int(element[1])
#			else :
#				#direction == "-"
#				offset = int(element[2]) + int(element[6])
#			if new_seq_id not in translation_db :
#				translation_db[new_seq_id] = []
#			translation_db[new_seq_id][(new_seq_start , new_seq_stop)] = [ old_seq_id , offset , direction ]
#
#	return translation_db


def make_agp_from_list( querylist , queryfastalen , gaplen , seqoutname , outfilename) :

	############################################################################################################################################################
	# AGP file format, tab delimited, 1-base coordinates #######################################################################################################
	############################################################################################################################################################
	#	##agp-version	2
	#	# 0                1          2          3        4            5                   6                    7                  8
	#	# Obj_Name	       Obj_start  Obj_End	 PartNum  Compnt_Type  CompntId/GapLength  CompntStart_GapType  CompntEnd_Linkage  Orientation_LinkageEvidence
	#	Super-Scaffold_9   1	      2          1	      N  	       2	               scaffold	            yes	               map
	#	Super-Scaffold_9   3          2581380    2	      W	           SMEL_DRAFT_562	   1	                2581378	           -
	#	Super-Scaffold_9   2581381    2716079    3	      N	           134699	           scaffold	            yes	               map
	############################################################################################################################################################
	#	0 - Obj_Name	- Object Out_sequence name (options.out + "_{1,2}_" + Target_sequence)
	#	1 - Obj_start 	- Object start on Out_sequence (start form 1, = 1 + prev_Obj_End)
	#	2 - Obj_End		- Object end on Out_sequence ( start + CompntEnd - CompntStart | start + GapLength - 1)
	#	3 - PartNum		- Progressive element count (start 1)
	#	4 - Compnt_Type 				- (W -> sequence) 						| (U -> Gap unknown size) [not use (N -> Gap known size) ]
	#	5 - CompntId/GapLength			- query_seq name 						| gap_length
	#	6 - CompntStart_GapType			- query_seq start [ == 1]				| gap_type [== scaffold]
	#	7 - CompntEnd_Linkage			- query_seq end [ == query_seq(len)]	| Linkage [== yes]
	#	8 - Orientation_LinkageEvidence	- strand								| Kind of evidence linking [align_genus]
	############################################################################################################################################################

	### Query list element format
	# Gap:	[61252	,	(0:61252)		,	gap		,	61252	,	0]
	#		[length	,	(T_start:Tstop)	, 	"gap" 	, 	length 	, 	0]
	# Object:	[b40-14.HS_iter4_seq3733|+	,	(61252:6463804)	,	93612:7595148	,	-6402552			,	4526208]
	#			[ID|strand					,	(T_start:Tstop)	,	Q_start:Q_stop	,	-(alignment length)	,	matches]

	gaplen = int(gaplen)
	partNum = 0
	Obj_Name = seqoutname
	Obj_end = 0

	for Id in querylist :
		CompntId_name = Id[:-2]
		Orientation = Id[-1]
		partNum += 1
		# Print gap if not printing the first sequence, print a gap
		if partNum > 1 :
			Obj_start = Obj_end + 1
			Obj_end = Obj_start + gaplen - 1
			print >> outfilename, "\t".join(str(x) for x in [ Obj_Name , str(Obj_start) , str(Obj_end) , str(partNum) , "N" , str(gaplen) , "scaffold" , "yes" , "align_genus" ])
			partNum += 1
		# Print the component
		CompntStart = 1
		CompntEnd = int(queryfastalen[CompntId_name])
		Obj_start = Obj_end + 1
		Obj_end = Obj_start + CompntEnd - CompntStart
		print >> outfilename, "\t".join(str(x) for x in [ Obj_Name , str(Obj_start) , str(Obj_end) , str(partNum) , "W" , CompntId_name , CompntStart , CompntEnd , Orientation])


def get_original_from_agp( agp_filename , fasta_dict ) :

	agp_dict = read_agp(agp_filename)
	sequence_pos_agp , gap_agp = agp_split( agp_dict )
	sequence_pos_dict = agp2range( sequence_pos_agp, "new")
	gap_dict = agp2range( gap_agp , "new")
	original_fasta_dict = agp2fasta( agp_dict , fasta_dict , "new_to_old" )
	original_gap_dict = get_agp_from_gap_list( get_gap_from_fasta_db( original_fasta_dict ) , original_fasta_dict )

	return agp_dict , sequence_pos_dict , gap_dict , original_fasta_dict , original_gap_dict


def agp_split( agp_dict ) :
	agp_seq = {}
	agp_gap = {}

	for seq_id in agp_dict :
		agp_seq[seq_id] = {}
		agp_gap[seq_id] = {}
		for start in sorted(agp_dict[seq_id].keys()) :
			if agp_dict[seq_id][start][4] == "W" :
				agp_seq[seq_id][start] = agp_dict[seq_id][start]
			else :
				agp_gap[seq_id][start] = agp_dict[seq_id][start]

	return agp_seq , agp_gap


def agp_ungapped_to_gap(agp_db, exclusion_db , min_length ) :
	# Input format:		agp_db[seq_id][int(start)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation ]
	#					exclusion_db[chr][ [chr, int(start) , int(end) , gene_id] ... [...]]
	# Output format:	regions_to_gaps_db[chr][(start,stop)] = {"upstream":(up_gap_start,up_gap_stop) , "downstream":(down_gap_start,down_gap_stop) }
	regions_to_gaps_db = {}
	for chr in agp_db :
		regions_to_gaps_db[chr] = {}
		prev_gap = [0,0]
		prev_region = ""
		for start in sorted(agp_db[chr].keys()) :
			#print >> sys.stderr , agp_db[chr][start]
			Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = agp_db[chr][start]
			# first item
			if Compnt_Type == "W" :
				# seq
				seq_region = ( int(Obj_start) , int(Obj_End) )
				if prev_region == "" :
					# first ungapped region, no gap before -> start
					regions_to_gaps_db[chr][seq_region] = {}
					regions_to_gaps_db[chr][seq_region]["upstream"] = prev_gap
					regions_to_gaps_db[chr][seq_region]["downstream"] = ""
					prev_region = seq_region
				else :
					if not prev_gap == "JUMP" :
						# prev_gap == "" |-> first ungapped region
						# prev_gap == (start,stop) |-> new ungapped region
						regions_to_gaps_db[chr][seq_region] = {}
						regions_to_gaps_db[chr][seq_region]["upstream"] = prev_gap
						regions_to_gaps_db[chr][seq_region]["downstream"] = ""
						prev_region = seq_region
					else :
						# prev_gap == "JUMP" |-> the element before was a skipped gap
						# Update the region in the DB extending it tup to the end of this region
						try:
							new_seq_region = ( int(prev_region[0]) , int(seq_region[1]) )
						except :
							print >> sys.stderr , "Error in previos region. Quitting!"
							print >> sys.stderr , prev_region
							sys.exit(1)
						regions_to_gaps_db[chr][new_seq_region] = {}
						regions_to_gaps_db[chr][new_seq_region]["upstream"] = regions_to_gaps_db[chr][prev_region]["upstream"]
						regions_to_gaps_db[chr][new_seq_region]["downstream"] = ""
						del regions_to_gaps_db[chr][prev_region]
						prev_region = new_seq_region
			else :
				# gap
				end = int(Obj_End)
				start = int(Obj_start)
				if exclusion_db == "[]" :
					excluded_regions = []
				else :
					try:
						excluded_regions = exclusion_db[chr]
					except :
						excluded_regions = []

				# test if at the beginnig
				if prev_region == "" :
					prev_gap = [start,end]
				else :
					if (start - end + 1) < int(min_length) or is_region_excluded( start , end , excluded_regions ) :
						# gap not suitable for breaking -> jump it
						prev_gap = "JUMP"
					else :
						regions_to_gaps_db[chr][prev_region]["downstream"] = [start,end]
						prev_gap = [start,end]

		if regions_to_gaps_db[chr][prev_region]["downstream"] == "" :
			# Close last item with end gap as there is no gap from the las ungapped region and the end
			regions_to_gaps_db[chr][prev_region]["downstream"] = [ int(prev_region[1]) +1 , int(prev_region[1]) +1 ]

	return regions_to_gaps_db


def is_region_excluded( start , stop , exclusion_db ) :
	# Input: exclusion_db = [ [chr, int(start) , int(end) , gene_id] ... [...]]
	excluded = False
	if not exclusion_db == [] :
		# If there is any region to exclude, test it
		for element in sorted(exclusion_db) :
			if element[1] <= int(start) <= element[2] or element[1] <= int(stop) <= element[2] :
				# the searched locus start or end falls within an excluded region
				excluded = True
	return excluded


def agp2range( agp_db , mode="new") :
	range_dict = {}
	for seq_id in agp_db :
		range_dict[seq_id] = []

		if mode == "new" :
			range_dict[seq_id] = []
			for start in sorted(agp_db[seq_id].keys()) :
				seq_start = int(agp_db[seq_id][start][1]) - 1
				seq_stop = int(agp_db[seq_id][start][2])
				range_dict[seq_id].append([ seq_id , seq_start , seq_stop ])
		elif mode == "old" :
			for start in sorted(agp_db[seq_id].keys()) :
				if agp_db[seq_id][start][4] == "W" :
					old_seq_id = agp_db[seq_id][start][5]
					old_seq_start = int(agp_db[seq_id][start][6]) - 1
					old_seq_stop = int(agp_db[seq_id][start][7])
					range_dict[seq_id].append([ old_seq_id , old_seq_start , old_seq_stop ])
		else :
			print >> sys.stderr , "[ERROR] Unknown range extraction mode"
			sys.exit(5)

	return range_dict


def test_range_uniqueness( agp_dict ) :
	overlapping = {}

	# create db of all old sequences ranges used
	all = {}
	for seq_id in agp_dict :
		for start in sorted(agp_dict[seq_id].keys()) :
			if agp_dict[seq_id][start][4] == "W" :
				Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = agp_dict[seq_id][start]
				if CompntId not in all :
					all[CompntId] = {}
				all[CompntId][( int(CompntStart) , int(CompntEnd) )] = [ CompntId , CompntStart , CompntEnd ,  Orientation , Obj_Name , Obj_start , Obj_End ]

	# find overlapping ranges
	for CompntId in all:
		if len(all[CompntId].keys()) == 1 :
			# CompntId present in only one range
			# No overlap possible
			# Skip
			continue
		else:
			# Multiple ranges from the same CompntId
			# Check if any are overlapping
			all_ranges = sorted(all[CompntId].keys())
			for range_ids in combinations(range(len(all_ranges)),2) :
				range1 = all_ranges[range_ids[0]]
				range2 = all_ranges[range_ids[2]]
				# Check if any range is overlapping the following ones
				# In case, add it to the overlapping dict
				if range1[0] <= range2[0] < range1[1] :
					# range2 starts within range1 -> Overlap
					# add both overlapping ranges to overlapping database
					if CompntId not in overlapping :
						overlapping[CompntId] = {}
					if range1 not in overlapping[CompntId] :
						overlapping[CompntId][range1] = all[CompntId][range1]
					if range2 not in overlapping[CompntId] :
						overlapping[CompntId][range2] = all[CompntId][range2]

	return overlapping


def print_overlapping_info(overlapping_ranges_dict, out_file_name) :
	out_file = open(out_file_name , "w")
	for id in overlapping_ranges_dict:
		print >> out_file , ">" + str(id)
		for range in sorted(overlapping_ranges_dict[id].keys()) :
			print >> out_file , "\t".join(overlapping_ranges_dict[id][range])
	out_file.close()


def get_broken_agp(sequences, regions_to_break, prefix) :
	# given a fasta sequence list and the pairs of gaps to break, generate new agp and broken sequences db
	# regions_to_break structure:
	# regions_to_break[seq_id][ [ breakpoint_type , start_gap , stop_gap ] , [...] , ... ]
	# breakpoint_type should be of 2 categories only: "gap" or "extremity"
	#	if "sequence" found -> raise error ### [BUG] breakpoint region search should have solved it

	new_sequences = {}
	new_agp_db = {}

	for seq_id in sequences :
		# agp_db content:
		# agp_db[seq_id][int(start)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation ]
		seq_len = int(len(sequences[seq_id]))

		if not seq_id in regions_to_break :
			# seq_id doesn't need to be broken,
			# add it as it is to the agp
			new_id = prefix + seq_id
			new_agp_db[new_id] = {}
			new_agp_db[new_id][1] = [ prefix + seq_id , "1" , str(seq_len) , "1" , "W" , seq_id , "1" , str(seq_len) , "+"]
		else :
			# seq_id must be broken around each gap
			# Check for extremity entry: begin == [ seq_id , 0 , 0 ] or end == [ seq_id , seq_len , seq_len ]
			chunk_num = 0
			# generete list of components
			component = {}
			next_CompntStart = ""
			prev_stop = 0
			for break_region in sorted(regions_to_break[seq_id]) :
				print >> sys.stderr, break_region
				try :
					type_gap , start_gap , stop_gap  = break_region
				except :
					sys.exit(2)
				if type_gap == "extremity" :
					#print >> sys.stderr, "Extremity"
					continue
				else :
					# Intrested in gaps only
					chunk_num += 1
					# Components within a break region are enclosed between
					# 	start: the end of left gap (start_gap_stop_pos)
					#	stop: the beginning of the right gap (stop_gap_start_pos)
					# Components in-between-break
					# 	start: after the end of the actual break region right gap (actual.stop_gap_stop_pos)
					# 	stop: before the beginning of the following break region left gap (next.start_gap_start_pos)
					CompntStart = int(prev_stop) + 1
					CompntEnd = int(start_gap) - 1
					component[int(chunk_num)] = [int(CompntStart) , int(CompntEnd)]
					prev_stop = stop_gap
					#print >> sys.stderr, component[int(PartNum)]

			chunk_num += 1
			CompntStart = int(prev_stop) + 1
			component[chunk_num] = [CompntStart , int(seq_len)]

			# generate agp entries:
			for chunk_num in sorted(component.keys()) :
				new_id = prefix + seq_id + "_" + str(chunk_num)
				start , stop = component[chunk_num]
				new_agp_db[new_id] = {}
				new_agp_db[new_id][1] = [ new_id , "1" , str(int(stop) - int(start) + 1 ) , "1" , "W" , seq_id , str(start) , str(stop) , "+"]

	new_sequences = agp2fasta( new_agp_db , sequences , "old_to_new")

	return new_sequences , new_agp_db


def agp_translate_agp( new_to_common_agp_db , common_to_old_agp_db ) :
	# agp_db[chr][int(start)]=
	#	Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId  , CompntStart , CompntEnd , Orientation
	#	Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , GapLength , GapType     , Linkage   , LinkageEvidence

	# new_to_common_agp_db: new_seq_id <- seq_id
	# common_to_old_agp_db: seq_id <- old_seq_ids
	# output: old_to_new_agp_db : new_seq_id <- old_seq_ids
	old_to_new_agp_db = {}
	for chr_id in sorted(new_to_common_agp_db.keys()) :
		#print >> sys.stderr ,  "- chr_id:" + chr_id
		old_to_new_PartNum = 0
		component_list = {}

		for block_start in sorted(new_to_common_agp_db[chr_id].keys()) :
			new_Obj_Name , new_Obj_start , new_Obj_End , new_PartNum , common_Compnt_Type , common_CompntId  , common_CompntStart , common_CompntEnd , common_Orientation = new_to_common_agp_db[chr_id][block_start]
			if not common_Compnt_Type == "W" :
				# The block is a gap, add it to the list as it is
				old_to_new_PartNum += 1
				component_list[old_to_new_PartNum] = [str(old_to_new_PartNum) , common_Compnt_Type , common_CompntId  , common_CompntStart , common_CompntEnd , common_Orientation ]

			else :
				# The block is a sequence region:
				# 	retrieve in common_to_old_agp_db the legacy components that do compose the component
				#	update the orientation
				#	trim unused regions
				if common_CompntId not in common_to_old_agp_db :
					print >> sys.stdout , "[WARNING] Sequence id missing from common to legacy sequences AGP file, impossible to complete the structure translation"
					print >> sys.stderr , "[WARNING] Sequence id " + common_CompntId + " is missing from common to legacy sequences AGP file, impossible to complete the structure translation"
					return {}
				else :
					legacy_components = common_to_old_agp_db[common_CompntId]
					used_legacy_components = {}

					# search for used parts of components
					cap_left = int(common_CompntStart)
					cap_right = int(common_CompntEnd)
					#print >> sys.stderr ,  "- region:" + common_CompntId + ":" + str(cap_left) + "-" + str(cap_right) + " | " + str(new_to_common_agp_db[chr_id][block_start])

					for common_start in sorted(legacy_components.keys()) :
						common_Obj_Name , common_Obj_start , common_Obj_End , legacy_PartNum , legacy_Compnt_Type , legacy_CompntId  , legacy_CompntStart , legacy_CompntEnd , legacy_Orientation = legacy_components[common_start]
						if ( cap_left <= int(common_Obj_start) <= cap_right ) and ( cap_left <= int(common_Obj_End) <= cap_right ) :
							#print >> sys.stderr ,  "-- legacy_component : " + str(legacy_components[common_start])
							#print >> sys.stderr ,  "--- legacy component entirely contained in the selected region "
							# Entire legacy component should be used
							used_legacy_components[int(legacy_PartNum)] = legacy_Compnt_Type , legacy_CompntId  , legacy_CompntStart , legacy_CompntEnd , legacy_Orientation
						elif ( cap_right < int(common_Obj_start) ) or ( int(common_Obj_End) < cap_left ) :
							# Entire legacy component before or after the selected region
							continue
						else :
							if legacy_Compnt_Type == "W" :
								#print >> sys.stderr ,  "-- legacy_component : " + str(legacy_components[common_start])
								#print >> sys.stderr ,  "--- legacy component should be used only partially as extends beyond the selected region "
								# trim left and/or right to cap_upper and/or cap_lower
								# left_delta and right_delta: If positive -> need to trim, otherwise use as is
								left_delta = cap_left - int(common_Obj_start)
								#print >> sys.stderr ,  "--- left delta:" + str(left_delta)
								if left_delta > 0 :
									# If positive -> need to trim as extends ahead of the selected region, otherwise no need to edit as starts with the selected region
									if legacy_Orientation == "+" :
										try:
											# trim beginning
											legacy_CompntStart = str(int(legacy_CompntStart) + left_delta)
										except :
											print >> sys.stdout ,  "[ERROR] Error translating AGP components: selected region error"
											print >> sys.stderr ,  "[ERROR] Error translating AGP components: faulty line in agp file, column 7 expecting an integer that was not found"
											print >> sys.stderr ,  "line           >> " + "\t".join([str(x) for x in new_to_common_agp_db[chr_id][block_start] ])
											print >> sys.stderr ,  "Legacy element >> " + "\t".join([str(x) for x in legacy_components[common_start] ])
											print >> sys.stderr ,  "left_delta     >> " + str(left_delta)
											sys.exit(1)
									else :
										try:
											# reversed orientation, trim the end
											legacy_CompntEnd = str(int(legacy_CompntEnd) - left_delta)
										except:
											print >> sys.stdout ,  "[ERROR] Error translating AGP components: selected region error"
											print >> sys.stderr ,  "[ERROR] Error translating AGP components: faulty line in agp file, column 8 expecting an integer that was not found"
											print >> sys.stderr ,  "line           >> " + "\t".join([str(x) for x in new_to_common_agp_db[chr_id][block_start] ])
											print >> sys.stderr ,  "Legacy element >> " + "\t".join([str(x) for x in legacy_components[common_start] ])
											print >> sys.stderr ,  "left_delta >> " + str(left_delta)
											sys.exit(1)


								right_delta = int(common_Obj_End) - int(cap_right)
								#print >> sys.stderr ,  "--- right delta:" + str(right_delta)
								if right_delta > 0 :
									# If positive -> need to trim as extends beyond the selected region, otherwise no need to edit as ends within the selected region
									if legacy_Orientation == "+" :
										# Trim the end
										try:
											legacy_CompntEnd = str(int(legacy_CompntEnd) - right_delta)
										except :
											print >> sys.stdout ,  "[ERROR] Error translating AGP components: selected region error"
											print >> sys.stderr ,  "[ERROR] Error translating AGP components: faulty line in agp file, column 8 expecting an integer that was not found"
											print >> sys.stderr ,  "line           >> " + "\t".join([str(x) for x in new_to_common_agp_db[chr_id][block_start] ])
											print >> sys.stderr ,  "Legacy element >> " + "\t".join([str(x) for x in legacy_components[common_start] ])
											print >> sys.stderr ,  "right_delta >> " + str(right_delta)
											sys.exit(1)
									else :
										try:
											# trim beginning
											legacy_CompntStart = str(int(legacy_CompntStart) + right_delta)
										except :
											print >> sys.stdout ,  "[ERROR] Error translating AGP components: selected region error"
											print >> sys.stderr ,  "[ERROR] Error translating AGP components: faulty line in agp file, column 7 expecting an integer that was not found"
											print >> sys.stderr ,  "line >> " + "\t".join([str(x) for x in new_to_common_agp_db[chr_id][block_start] ])
											print >> sys.stderr ,  "element >> " + "\t".join([str(x) for x in legacy_components[common_start] ])
											print >> sys.stderr ,  "right_delta >> " + str(right_delta)
											sys.exit(1)
								used_legacy_components[int(legacy_PartNum)] = legacy_Compnt_Type , legacy_CompntId  , legacy_CompntStart , legacy_CompntEnd , legacy_Orientation
								#print >> sys.stderr ,  "---- New component: " + str(used_legacy_components[int(legacy_PartNum)])
							else :
								left_delta = cap_left - int(common_Obj_start)
								right_delta = int(common_Obj_End) - int(cap_right)
								print >> sys.stdout ,  "[ERROR] Error translating AGP components: selected region error"
								print >> sys.stderr ,  "[ERROR] Error translating AGP components: the selected region includes part of a gap sequence at one extremity. Region needs refinement"
								print >> sys.stderr ,  "AGP line         >> " + "\t".join([str(x) for x in new_to_common_agp_db[chr_id][block_start] ])
								print >> sys.stderr ,  "Block region     >> " + "\t".join([str(x) for x in [ new_to_common_agp_db[chr_id][block_start][5] , int(new_to_common_agp_db[chr_id][block_start][6]) - 1 , new_to_common_agp_db[chr_id][block_start][7] , new_to_common_agp_db[chr_id][block_start][8] ] ] )
								print >> sys.stderr ,  "Gap involved     >> " + "\t".join([str(x) for x in legacy_components[common_start] ])
								gap_length = int(legacy_components[common_start][5])
								if left_delta > 0 :
									print >> sys.stderr ,  "Left trim        >> " + str(gap_length - left_delta) + "bp"
									if new_to_common_agp_db[chr_id][block_start][8] == "+" :
										print >> sys.stderr ,  "New Block region >> " + "\t".join([str(x) for x in [ new_to_common_agp_db[chr_id][block_start][5] , int(new_to_common_agp_db[chr_id][block_start][6]) - 1 + (gap_length - left_delta) , new_to_common_agp_db[chr_id][block_start][7] , new_to_common_agp_db[chr_id][block_start][8] ] ] )
									if new_to_common_agp_db[chr_id][block_start][8] == "-" :
										print >> sys.stderr ,  "New Block region >> " + "\t".join([str(x) for x in [ new_to_common_agp_db[chr_id][block_start][5] , int(new_to_common_agp_db[chr_id][block_start][6]) - 1 , new_to_common_agp_db[chr_id][block_start][7] - (gap_length - right_delta) , new_to_common_agp_db[chr_id][block_start][8] ] ] )
								elif right_delta > 0 :
									print >> sys.stderr ,  "Right trim       >> " + str(gap_length - right_delta) + "bp"
									if new_to_common_agp_db[chr_id][block_start][8] == "-" :
										print >> sys.stderr ,  "New Block region >> " + "\t".join([str(x) for x in [ new_to_common_agp_db[chr_id][block_start][5] , int(new_to_common_agp_db[chr_id][block_start][6]) - 1 + (gap_length - left_delta) , new_to_common_agp_db[chr_id][block_start][7] , new_to_common_agp_db[chr_id][block_start][8] ] ] )
									if new_to_common_agp_db[chr_id][block_start][8] == "+" :
										print >> sys.stderr ,  "New Block region >> " + "\t".join([str(x) for x in [ new_to_common_agp_db[chr_id][block_start][5] , int(new_to_common_agp_db[chr_id][block_start][6]) - 1 , new_to_common_agp_db[chr_id][block_start][7] - (gap_length - right_delta) , new_to_common_agp_db[chr_id][block_start][8] ] ] )

								sys.exit(1)
								#delta = left_delta + right_delta
								#new_length = str(int(legacy_CompntId) - delta)
								#used_legacy_components[int(legacy_PartNum)] = legacy_Compnt_Type , new_length  , legacy_CompntStart , legacy_CompntEnd , legacy_Orientation

					#print >> sys.stderr ,  "-- legacy components to add: " + str(len(used_legacy_components))
					# Add selected regions with the correct orientation
					if common_Orientation == "+" :
						# Orientation has been kept
						for component in sorted(used_legacy_components.keys()) :
							old_to_new_PartNum += 1
							legacy_Compnt_Type , legacy_CompntId  , legacy_CompntStart , legacy_CompntEnd , legacy_Orientation = used_legacy_components[component]
							component_list[old_to_new_PartNum] = [str(old_to_new_PartNum) , legacy_Compnt_Type , legacy_CompntId  , legacy_CompntStart , legacy_CompntEnd , legacy_Orientation ]

					elif common_Orientation == "-" :
						# Orientation has been inverted, invert all sequences
						for component in sorted(used_legacy_components.keys() , reverse=True) :
							old_to_new_PartNum += 1
							legacy_Compnt_Type , legacy_CompntId  , legacy_CompntStart , legacy_CompntEnd , legacy_Orientation = used_legacy_components[component]
							if legacy_Orientation == "+" :
								legacy_Orientation = "-"
							else:
								legacy_Orientation = "+"
							component_list[old_to_new_PartNum] = [str(old_to_new_PartNum) , legacy_Compnt_Type , legacy_CompntId  , legacy_CompntStart , legacy_CompntEnd , legacy_Orientation ]

					else :
						print >> sys.stdout , "[WARNING] Sequence missing orientation information, impossible to complete the structure translation"
						print >> sys.stderr , "[WARNING] Sequence missing orientation information, impossible to complete the structure translation"
						return {}

		# Reconstruct new_agp coordinates on the component_list content
		old_to_new_agp_db[chr_id] = {}
		Obj_Name = chr_id
		Obj_start = 0
		Obj_End	= 0
		for PartNum in sorted(component_list.keys()) :
			Obj_start = Obj_End +1
			old_to_new_PartNum , common_Compnt_Type , common_CompntId  , common_CompntStart , common_CompntEnd , common_Orientation = component_list[PartNum]
			if not common_Compnt_Type == "W" :
				Obj_End = Obj_start + int(common_CompntId) - 1
			else :
				Obj_End = Obj_start + int(common_CompntEnd) - int(common_CompntStart)
			old_to_new_agp_db[chr_id][Obj_start] = [ Obj_Name , str(Obj_start) , str(Obj_End) , old_to_new_PartNum , common_Compnt_Type , common_CompntId  , common_CompntStart , common_CompntEnd , common_Orientation ]

	return old_to_new_agp_db


def agp_translate_agp_translate_agp( new_agp_db , agp_db , legacy_agp_db , legacy_seq_len_db ):
	# new_agp_db: new_seq_id <- seq_id
	# agp_db: seq_id <- legacy_ids
	# legacy_agp_db: new_legacy_ids <- legacy_ids
	# output: new_seq_id <- new_legacy_ids

	# invert legacy_agp_db
	# legacy_agp_db_rev :  legacy_ids <- new_legacy_ids
	legacy_agp_db_rev = invert_agp( legacy_agp_db , False , legacy_seq_len_db )
	# convert: (seq_id <- legacy_ids) <- (legacy_ids <- new_legacy_ids) => (seq_id <- new_legacy_ids)
	intermediate_agp = agp_translate_agp( agp_db , legacy_agp_db_rev )
	# convert: (new_seq_id <- seq_id) <- (seq_id <- new_legacy_ids) => (new_seq_id <- new_legacy_ids)
	old_to_new_agp_db = agp_translate_agp( new_agp_db , intermediate_agp )

	# TODO: Double checking for novel breakpoints

	return old_to_new_agp_db


def bed_to_agp_oneliner( bed_dict , prefix ) :
	agp_db = {}
	# generate AGP from BED file
	# AGP_db format:
	# agp_db[seq_id][int(start)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation ]
	# BED format:
	# [ chrom , chromStart(0->) , chromEnd(1->) , name , score , strand , ... ]
	# ignore file fields 4,5,7->
	# If strand is "." -> "+"
	# Generate a new entry for each bed range
	for seq_id in bed_dict :
		for entry_id in sorted(bed_dict[seq_id].keys()) :
			entry = bed_dict[seq_id][entry_id]
			CompntStart = entry[1] + 1
			CompntEnd = entry[2]
			Obj_Name = prefix + seq_id + "-" + str(CompntStart) + "_" + str(CompntEnd)
			Obj_End = 1 + int(CompntEnd) - int(CompntStart)

			Orientation = "+"
			if ( len(entry) > 5 ) and ( entry[6] == "-" ) :
				Orientation = "-"

			agp_db[Obj_Name][1] = [ Obj_Name , "1" , Obj_End , "1" , "W" , seq_id , CompntStart , CompntEnd ,  Orientation ]

	return agp_db


def bed_to_agp_onefile( bed_dict , gap_size , out_seq_name ) :
	agp_db = {}
	# generate AGP from BED file
	# AGP_db format:
	# agp_db[seq_id][int(start)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation ]
	# BED format:
	# bed_dict[num][ chrom , chromStart(0->) , chromEnd(1->) , name , score , strand , ... ]
	# ignore file fields 3,4,6->
	# If strand is "." -> "+"
	# Generate a new entry for each bed range
	Obj_end = 0
	partNum = 0
	Obj_Name = out_seq_name
	gaplen = int(gap_size)
	agp_db[Obj_Name] = {}
	for entry_id in sorted(bed_dict.keys()) :
		entry = bed_dict[entry_id]
		component_id = entry[0]
		CompntStart = int(entry[1]) + 1
		CompntEnd = int(entry[2])
		try :
			Orientation = entry[5]
		except :
			print >> sys.stderr, "[WARNING] BED line has no strand information, +(plus) will be used: " + "\t".join( [str(x) for x in entry] )
			Orientation = "+"
		if Orientation == "." :
			print >> sys.stderr, "[WARNING] BED line has no strand information, +(plus) will be used: " + "\t".join( [str(x) for x in entry] )
			Orientation = "+"

		# add gap if not printing the first sequence, print a gap
		if partNum > 0 :
			partNum += 1
			Obj_start = Obj_end + 1
			Obj_end = Obj_start + gaplen - 1
			agp_db[Obj_Name][int(Obj_start)] =  [ Obj_Name , Obj_start , Obj_end , partNum , "N" , gaplen , "scaffold" , "yes" , "align_genus" ]

		# add the component
		partNum += 1
		Obj_start = Obj_end + 1
		Obj_end = Obj_start + CompntEnd - CompntStart
		agp_db[Obj_Name][int(Obj_start)] = [ Obj_Name , Obj_start , Obj_end , partNum , "W" , component_id , CompntStart , CompntEnd ,  Orientation ]

	return agp_db


def agp_to_block( agp_db ) :
	block_db = {}

	for chr_id in sorted(agp_db.keys()) :
		block_db[chr_id] = {}
		block_id = 0
		for comp_start in sorted( agp_db[chr_id].keys() ) :
			Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId  , CompntStart , CompntEnd , Orientation = agp_db[chr_id][comp_start]
			if not Compnt_Type == "W" :
				continue
			else :
				block_id += 1
				block_db[chr_id][block_id] = [ CompntId  , int(CompntStart) - 1 , int(CompntEnd) , Orientation ]
	return block_db


def block_to_agp( block_db , gap_size , prefix = "" , substitute_id=False) :
	# Input agp_db format:
	#	agp_db[OutSeqID_1]
	#	#	agp_db[OutSeqID_1][int(start_1)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type [=="W"] , CompntId  , CompntStart , CompntEnd , Orientation     ]
	# 	#	agp_db[OutSeqID_1][int(start_2)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type [!="W"] , GapLength , GapType     , Linkage   , LinkageEvidence ]
	#	# start is 1 based!
	#
	#
	# Output block_db format:
	#	block_db[OutSeqID_1]
	#	#	block_db[OutSeqID_1][1] = [SeqID_1 , Start , Stop , Strand]
	#	#	block_db[OutSeqID_1][2] = [SeqID_2 , Start , Stop , Strand]
	#	#	...
	#	#	block_db[OutSeqID_1][N] = [SeqID_N , Start , Stop , Strand]
	#	block_db[OutSeqID_2]
	#	#	block_db[OutSeqID_2][1] = [SeqID_X , Start , Stop , Strand]
	#	#	block_db[OutSeqID_2][2] = [SeqID_Y , Start , Stop , Strand]
	#	#	...
	#	#	block_db[OutSeqID_2][M] = [SeqID_Z , Start , Stop , Strand]
	#	Start is zero based!

	agp_db = {}
	counter = 0
	for seq_name in sorted(block_db.keys()) :
		counter += 1
		if substitute_id :
			Obj_Name = prefix + str(counter)
		else :
			Obj_Name = prefix + seq_name
		agp_db[Obj_Name] = {}
		blocks = block_db[seq_name]
		gaplen = int(gap_size)
		partNum = 0
		Obj_end = 0

		for block_num in sorted(blocks.keys()) :
			seq_id , Start , Stop , Orientation = blocks[block_num]
			CompntStart = int(Start) + 1
			CompntEnd = int(Stop)
			# add gap if not printing the first sequence, print a gap
			if partNum > 0 :
				partNum += 1
				Obj_start = Obj_end + 1
				Obj_end = Obj_start + gaplen - 1
				agp_db[Obj_Name][int(Obj_start)] = [Obj_Name , Obj_start , Obj_end , partNum , "N" , gaplen , "scaffold" , "yes" , "align_genus" ]

			# add the component
			partNum += 1
			Obj_start = Obj_end + 1
			Obj_end = Obj_start + CompntEnd - CompntStart
			agp_db[Obj_Name][int(Obj_start)] = [ Obj_Name , Obj_start , Obj_end , partNum , "W" , seq_id , CompntStart , CompntEnd ,  Orientation ]

	return agp_db


def get_joined( gap , agp_dict ) :
	found = False
	left_seq_id = ""
	right_seq_id = ""
	for start in sorted(agp_dict.keys()) :
		Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation  = agp_dict[start]
		if Compnt_Type=="W" :
			if not found :
				left_seq_id = CompntId
			else :
				if right_seq_id == "" :
					right_seq_id = CompntId
				break
		else :
			# search for the given gap
			if ( int(Obj_start) == int(gap[1]) ) and ( int(Obj_End) == int(gap[2]) ):
				found = True

	return left_seq_id , right_seq_id


def get_agp_from_gap_list( gap_db , seq_db ) :
	seq_len_db = get_length_from_fasta_db(seq_db)
	agp_db = {}
	# Input format: gap_db[chr][ [chr , gap_start(0 based) , gap_stop(1 based) ] ... [...] ]
	# Output format: agp_db[seq_id][int(start)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation ]
	# 																										 GapLength , GapType    , Linkage   , LinkageEvidence
	for chr in sorted(gap_db.keys()) :
		agp_db[chr] = {}
		prev_gap_stop = 0
		part_num = 0
		for gap in sorted(gap_db[chr]) :
			if gap[1] == 0 :
				# first element is a gap
				# add only the gap region to the AGP
				part_num += 1
				gap_agp_line = [ str(chr) , str(int(gap[1]) +1 ) , str(gap[2]) , str(part_num) , "N" , str(int(gap[2]) - int(gap[1]) ) , "scaffold" , "yes" , "map"]
				agp_db[chr][int(gap[1]) +1] = gap_agp_line
			else:
				# Add a sequence region between the precedent and the actual gap
				part_num += 1
				seq_line = [ str(chr) , str(int(prev_gap_stop) +1 ) , str(gap[1]) , str(part_num) , "W" , str(chr) + "_comp" + str(part_num) , "1" , str( int(gap[1]) - int(prev_gap_stop) ) , "+" ]
				agp_db[chr][int(prev_gap_stop) +1] = seq_line
				part_num += 1
				gap_agp_line = [ str(chr) , str(int(gap[1]) +1 ) , str(gap[2]) , str(part_num) , "N" , str(int(gap[2]) - int(gap[1]) ) , "scaffold" , "yes" , "map"]
				agp_db[chr][int(gap[1]) +1] = gap_agp_line

			prev_gap_stop = str(gap[2])

		if not prev_gap_stop == seq_len_db[chr] :
			# add a sequence region to cover until the end
			part_num += 1
			seq_line = [ str(chr) , str(int(prev_gap_stop) +1 ) , str(seq_len_db[chr]) , str(part_num) , "W" , str(chr) + "_comp" + str(part_num) , "1" , str( seq_len_db[chr] - int(prev_gap_stop) ) , "+" ]
			agp_db[chr][int(prev_gap_stop) +1] = seq_line

	return agp_db


def find_component_in_agp( pos , component_db ) :
	for component in sorted(component_db.keys()) :
		Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = component_db[component]
		if ( Compnt_Type == "W" ) and ( int(Obj_start) <= int(pos) <= int(Obj_End) ) :
			return CompntId


def ranges_to_agp( seq_id , ranges , prefix , name_prefix , gaplen = 10000 ) :
	# Always 3 ranges: Upstream , gap/gap_alt , downstream
	print >> sys.stderr , "##### " + name_prefix
	print >> sys.stderr , "###### Making AGP and FASTA from genomic regions"
	print >> sys.stderr , "####### gaplen = " + str(gaplen)
	#json.dump(ranges , sys.stderr , indent=4)
	Obj_Name = prefix
	CompntId_name = seq_id
	Orientation = "+"
	agp_file = name_prefix + ".agp"
	fasta_file = name_prefix + ".fasta"
	outfilename = open(agp_file , 'w')
	regions_dict = { Obj_Name : {} }
	region_types = [ "UP" , "SEQ" , "DOWN" ]
	partNum = 0
	Obj_start = 1
	Obj_end = 0
	new_sequence = ""
	full_signal = ""
	print >> sys.stderr , "####### AGP result: "
	for region in ranges :
		if not partNum == 0 :
			partNum+=1
			Obj_start = Obj_end + 1
			Obj_end = Obj_start + gaplen - 1
			gap_seq = "N"*gaplen
			gap_signal = "G"*gaplen
			new_sequence += gap_seq
			full_signal += gap_signal
			region_type = "GAP"
			Compnt_Type = "N"
			regions_dict[Obj_Name][partNum] = {
				"coords" : [ Obj_Name , int(Obj_start) , int(Obj_end) ] ,
				"region_type" : region_type ,
				"Compnt_Type" : Compnt_Type ,
				"desc" : gaplen
				}
			print >> outfilename, "\t".join(str(x) for x in [ Obj_Name , str(Obj_start) , str(Obj_end) , str(partNum) , Compnt_Type , str(gaplen) , "scaffold" , "yes" , "align_genus" ])
			print >> sys.stderr, "######## " + "\t".join(str(x) for x in [ Obj_Name , str(Obj_start) , str(Obj_end) , str(partNum) , Compnt_Type , str(gaplen) , "scaffold" , "yes" , "align_genus" ])

		partNum += 1
		region_type = region_types.pop(0)

		try :
			id, start , stop , seq , sig  = region
			CompntStart = int(start)
			CompntEnd = int(stop)
			Obj_start = Obj_end + 1
			Obj_end = Obj_start + CompntEnd - CompntStart
			regions_dict[Obj_Name][partNum] = {
				"coords" : [ Obj_Name , int(Obj_start) , int(Obj_end) ] ,
				"region_type" : region_type ,
				"Compnt_Type" : "W" ,
				"desc" : [ id , CompntStart , CompntEnd ,  Orientation ]
				}
			print >> outfilename, "\t".join(str(x) for x in [ Obj_Name , str(Obj_start) , str(Obj_end) , str(partNum) , "W" , id , CompntStart , CompntEnd , Orientation])
			print >> sys.stderr, "######## " + "\t".join(str(x) for x in [ Obj_Name , str(Obj_start) , str(Obj_end) , str(partNum) , "W" , id , CompntStart , CompntEnd , Orientation])
		except :
			#print >> sys.stderr , region
			substitute_gap = 10
			try:
				id, start , stop = region
			except :
				print >> sys.stderr , region
				sys.exit(1)
			CompntStart = int(start)
			CompntEnd = int(stop)
			Obj_start = Obj_end + 1
			Obj_end = Obj_start + substitute_gap - 1
			seq = "N"*substitute_gap
			sig = "G"*substitute_gap
			regions_dict[Obj_Name][partNum] = {
				"coords" : [ Obj_Name , int(Obj_start) , int(Obj_end) ] ,
				"region_type" : region_type ,
				"Compnt_Type" : "gap_substitute" ,
				"desc" : [ id , CompntStart , CompntEnd ,  Orientation ]
				}
			print >> outfilename, "\t".join(str(x) for x in [ Obj_Name , str(Obj_start) , str(Obj_end) , str(partNum) , "N" , str(substitute_gap) , "scaffold_" + id + ":" + str(start) + "-" + str(stop) , "yes" , "align_genus" ])
			print >> sys.stderr, "######## " + "\t".join(str(x) for x in [ Obj_Name , str(Obj_start) , str(Obj_end) , str(partNum) , "N" , str(substitute_gap) , "scaffold_" + id + ":" + str(start) + "-" + str(stop) , "yes" , "align_genus" ])
		new_sequence += seq
		full_signal += sig


	outfilename.close()
	outfasta = open(fasta_file , 'w')
	print >> outfasta, ">" + prefix
	print >> outfasta, new_sequence
	outfasta.close()

	print >> sys.stderr , "####### FASTA result: " + prefix + " | length: " + str(len(new_sequence))

	# regions_dict[Obj_Id]
	#	 regions_dict[Obj_Id][partNum]
	#	 	{
	#	 		"coords" : [ Obj_Id , int(Obj_start) , int(Obj_End) ] ,
	#			"region_type" : [ UP | SEQ | DOWN | GAP ]
	#	 		"Compnt_Type" : [ "N" , "W" , "gap_substitute" ]
	#	 		"desc" : [ CompntId_name , CompntStart , CompntEnd ,  Orientation ] / int(gap_len)
	#	 	}

	return agp_file , regions_dict , new_sequence , fasta_file , full_signal


def rename_agp_sequences( agp_db , prefix ) :
	new_agp_db = {}
	numeric_id = 0
	for chr in sorted( agp_db.keys() ) :
		numeric_id += 1
		new_name = prefix + "_" + str(numeric_id)
		new_agp_db[new_name] = {}
		for pos in sorted(agp_db[chr].keys()) :
			new_agp_db[new_name][int(pos)] = agp_db[chr][pos]
			new_agp_db[new_name][int(pos)][0] = new_name
	return new_agp_db


def add_block_to_agp(agp_db , chr , block , gap_size , spacer , refining_status ) :
	new_agp_dict = agp_db[chr]
	# AGP db format: agp_db[seq_id][int(start)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation ]
	# 																										 GapLength , GapType    , Linkage   , LinkageEvidence
	if new_agp_dict == {} :
		PartNum = 0
		Obj_End = 0
		Obj_Name = chr
	else :
		last_element = new_agp_dict[ max( new_agp_dict.keys() )  ]
		Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , prevd , prevstart , prevend ,  prevOrientation = last_element

	if "region_trimmed" in block :
		CompntId , blockStart , blockEnd , Orientation = block["region_trimmed"]
	elif "region_corrected" in block :
		CompntId , blockStart , blockEnd , Orientation = block["region_corrected"]
	else :
		CompntId , blockStart , blockEnd , Orientation = block["region_given"]

	CompntStart = int(blockStart) + 1
	CompntEnd = int(blockEnd)
	# Make new agp component from block
	PartNum = int(PartNum) + 1
	Obj_start = int(Obj_End) + 1
	Obj_End = int(Obj_start) + int(CompntEnd) - int(CompntStart)
	new_agp_dict[int(Obj_start)] = [ Obj_Name , Obj_start , Obj_End , PartNum , "W" , CompntId , CompntStart , CompntEnd , Orientation ]
	print >> sys.stderr , "#### Region saved: " + str([ CompntId , blockStart , blockEnd , Orientation ])
	if refining_status == "disjointed" :
		if int(gap_size) > 0 :
			# create a gap after the sequence
			PartNum = int(PartNum) + 1
			Obj_start = int(Obj_End) + 1
			Obj_End = Obj_start + int(gap_size) - 1
			new_agp_dict[int(Obj_start)] = [ Obj_Name , Obj_start , Obj_End , PartNum , "N" , int(gap_size) , "scaffold" , "yes" , "align_genus" ]
	elif refining_status == "overlapping" :
		if int(spacer) > 0 :
			# Add a spacer after the sequence
			PartNum = int(PartNum) + 1
			Obj_start = int(Obj_End) + 1
			Obj_End = Obj_start + int(spacer) - 1
			new_agp_dict[int(Obj_start)] = [ Obj_Name , Obj_start , Obj_End , PartNum , "N" , int(spacer) , "contig" , "yes" , "align_genus" ]
	elif not refining_status == "last" :
		print >> sys.stdout , "[ERROR] Unknown junction type"
		print >> sys.stderr , "[ERROR] Unknown junction type"
		print >> sys.stderr , "\"" + refining_status + "\" reported while acceptable values are [disjointed|overlapping|last] "
		sys.exit(1)

	return new_agp_dict


def clean_self_agp(agp_db) :
	new_agp = {}
	for seq_id in agp_db :
		for start in sorted(agp_db[seq_id].keys()) :
			component = agp_db[seq_id][start]
			Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = component
			if not Obj_Name == CompntId :
				if seq_id not in new_agp :
					new_agp[seq_id] = {}
				new_agp[seq_id][start] = component
			else :
				print >> sys.stderr , "[WARNING] AGP component for " + Obj_Name + " in the region " + str(Obj_start) + "-" + str(Obj_End) + " refers to itself from " + str(CompntStart) + "-" + str(CompntEnd) + ". Component removed"
	return new_agp

