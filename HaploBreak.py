#!/usr/bin/env python


import argparse
from lib_files.HaploFunct import *
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *
import sys

gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


def main() :

	###### Options and help ######

	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fasta", dest="fasta",
					help="Genome sequences in FASTA format. Use unmasked or soft masked sequences only. [REQUIRED]", metavar="genome.fasta")
	parser.add_argument("-b", "--breakpoints", dest="breaks",
					help="List of breakpoints to be applied [REQUIRED]. Tab separated format, 1) sequence id; 2)comma separated values of regions to break: start1-stop1[,start2-stop2, ... ,startN-stopN]", metavar="breakpoints.tsv")
	parser.add_argument("-o", "--out", dest="out", default="out" ,
					help="Output file names [Default: out]", metavar="outname")
	parser.add_argument("-p", "--prefix", dest="prefix", default="" ,
					help="New sequence names prefix to be added before the old name. Leave unset to keep the same names [Default: NONE]", metavar="new")
	parser.add_argument("-g", "--gff3", dest="gff3",
					help="Genome annotation in GFF3 format.", metavar="annotation.gff3")
	parser.add_argument("-a", "--agp", dest="agp",
					help="AGP file defining the coordinates of legacy sequences composing the input assembly", metavar="previous_to_actual.genome.agp")
	parser.add_argument("--bed", dest="bed",
					help="Region annotation in BED format. Accepted BED3 and BED6 formats, additional columns will be considered as annotation and reported unedited", metavar="annotation.bed")

	parser.add_argument("-d", "--distance", dest="distance", default="25000" ,
					help="Maximum distance (in bp) allowed to search precise breakpoint into junctions and trigger in FASTA gap search [Default: 25,000]", metavar="25000")
	parser.add_argument("-m", "--minimum", dest="min", default="200" ,
					help="Minimum gap size (in bp) on fasta for break point selection [Default: 200]", metavar="200")
	parser.add_argument("--allow_annotation_breaking", dest="break_annot", default=False, action="store_true",
					help="Allow breaking annotated genes if GFF3 is reported")

	scriptDirectory = os.path.dirname(os.path.realpath(__file__)) + "/support_scripts"
	print >> sys.stdout, "Running HaploBreaker tool from HaploSync version " + get_version()
	print >> sys.stdout, "To reproduce this run use the following command: " + " ".join( pipes.quote(x) for x in sys.argv)
	print >> sys.stdout, "----"
	# Sanity Check

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)

	options = parser.parse_args()

	if not options.fasta :
		# noinspection PyCompatibility
		print >> sys.stderr , "[ERROR] Genome FASTA file missing"
		parser.print_help()
		sys.exit(1)

	if not options.breaks :
		print >> sys.stderr , "[ERROR] Breakpoints file missing"
		parser.print_help()
		sys.exit(1)

	if not options.gff3 :
		print >> sys.stderr , "[NOTE] Annotation file missing: Coordinates conversion will not take place (you must run it downstream) and breakpoints selection will not trigger gene loci disruption protection"

	if options.agp :
		print >> sys.stderr , "[NOTE] AGP file given: breaking on path gap will be prioritized. Any breaking of original sequences will be notified. Output will be delivered also for original sequences"
	else :
		print >> sys.stderr, "[WARNING] No AGP file was given, \"N\" encoded gaps in the sequence will be used for searching safe breakpoints only"
	max_distance = int(options.distance)

	##### Reads inputs
	#### Read sequence
	print >> sys.stdout, "- Importing genome sequences"
	fasta_in = read_fasta(options.fasta)
	sequences = {}
	for id in fasta_in :
		sequences[id] = remove_islets_fasta(fasta_in[id].upper())

	print >> sys.stdout, "-- Calculating reference sequences lengths"
	sequences_len = get_length_from_fasta_db(sequences)

	print >> sys.stdout, "--- Number of reference sequences: " + str(len(sequences_len.keys()))

	#### Read breaking points
	print >> sys.stdout, "- Importing breaking points"
	breakpoint_db = {}
	for line in open( options.breaks , "r") :
		try :
			seq_id , breaks = line.rstrip().rstrip(",").split("\t")
		except :
			continue
		if seq_id not in sequences :
			print >> sys.stderr, "[[WARNING]] Breakpoints reported for sequence missing from FASTA file: [" + seq_id + "]. Data discarded"
			continue
		if breaks == "" :
			# No breakpoints
			continue
		else :
			# Sanity check
			breaks_list = breaks.split(",")
			for break_couple in breaks_list :
				try :
					start , stop = break_couple.rstrip().lstrip().split("-")
				except :
					print >> sys.stderr, "[ERROR] Error for breakpoint  coordinates in " + seq_id + " : " + str(break_couple.rstrip()) + " found, not a valid pair"
					sys.exit(1)
				if start.upper() == "BEGIN" or start.upper() == "START" :
					start = 0
				else :
					try :
						start = int(start)
					except :
						print >> sys.stderr, "[ERROR] Non integer coordinates for breakpoint in " + seq_id + " : " + str(start) +" found"
						sys.exit(1)
				if stop.upper() == "END" or stop.upper() == "STOP" :
					stop = sequences_len(seq_id)
				else :
					try:
						stop = int(stop)
					except :
						print >> sys.stderr, "[ERROR] Non integer coordinates for breakpoint in " + seq_id + " : " + str(stop) +" found"
						sys.exit(1)
				if seq_id not in breakpoint_db :
					breakpoint_db[seq_id] = []
				breakpoint_db[seq_id].append( [ start , stop ] )

	print >> sys.stdout, "-- Breaking points imported:"
	for seq_id in sorted(sequences.keys()) :
		break_line = ""
		if seq_id in breakpoint_db :
			for element in sorted(breakpoint_db[seq_id]) :
				break_line += "(" + str(element[0]) + ":" + str(element[1]) + "); "
		if not break_line == "" :
			print >> sys.stdout, "--- Sequence " + str(seq_id) + ": " + break_line

	#### - IF GIVEN - Read annotation
	if options.gff3:
		print >> sys.stdout, "- Importing GFF3"
		annotation_gff3 , mRNA_db = read_gff3(options.gff3)
		gene_db = feature_ranges(annotation_gff3 , "gene")
	else :
		annotation_gff3 = ""
		gene_db = {}
		for seq_id in sorted(sequences.keys()) :
			gene_db[seq_id] = []

	#### List gaps in FASTA sequences
	print >> sys.stdout, "- Extracting position of gaps in sequences"

	if options.gff3 and not options.break_annot :
		#gap_list = get_gap_from_fasta_db(sequences , {})
		#print >> sys.stderr , "gap list generated without gene masking"
		#print >> sys.stderr , gap_list
		gap_list = get_gap_from_fasta_db( sequences , gene_db )
		#print >> sys.stderr , "gap list generated with gene masking"
		#print >> sys.stderr , gap_list
	else :
		gap_list = get_gap_from_fasta_db(sequences , {})
		print >> sys.stderr , "gap list generated"

	gaps_file_name = options.out + ".breakable_gaps.list.txt"
	gaps_file = open(gaps_file_name , 'w')
	for seq in gap_list :
		for element in sorted(gap_list[seq]) :
			print >> gaps_file , "\t".join([ str(x) for x in element])
	gaps_file.close()

	gapped_agp_db = get_agp_from_gap_list(gap_list , sequences )
	ungapped_db = agp_ungapped_to_gap(gapped_agp_db , gene_db , int(options.min) )

	#### - IF GIVEN - Read AGP
	if options.agp :
		print >> sys.stdout, "- Importing AGP and original sequence information"
		agp_db , original_sequence_position , junctions , sequences_legacy , original_gap_list = get_original_from_agp( options.agp , sequences )
		sequences_legacy_length_db = get_length_from_fasta_db(sequences_legacy)
		#### Test uniqueness of old ids
		print >> sys.stdout, "-- Testing uniqueness AGP ids"
		overlapping_original_ranges = test_range_uniqueness(agp_db)
		if overlapping_original_ranges == {} :
			print >> sys.stdout, "--- None of the original ids was used more than once"
			translation_db = translate_from_AGP_whole_genome_reverse(agp_db)
			if options.gff3:
				print >> sys.stdout, "--- Converting annotation"
				legacy_annotation_gff3 = translate_gff3( annotation_gff3 , translation_db , options.out + ".broken_loci.given_to_legacy_transfer.txt" )
				write_gff3(legacy_annotation_gff3 , options.out + ".legacy.annotation.gff3" , get_length_from_fasta_db( sequences_legacy ) )
			else :
				legacy_annotation_gff3 = {}
				for seq_id in sequences_legacy :
					legacy_annotation_gff3[seq_id] = []
			legacy_gap_list = get_gap_from_fasta_db( sequences_legacy , legacy_annotation_gff3 )
			legacy_gaps_file_name = options.out + ".legacy_breakable_gaps.list.txt"
			legacy_gaps_file = open(legacy_gaps_file_name , 'w')
			for seq in sorted(legacy_gap_list.keys()) :
				for element in sorted(legacy_gap_list[seq]) :
					print >> legacy_gaps_file , "\t".join([ str(x) for x in element])
			legacy_gaps_file.close()

		else :
			print >> sys.stdout, "--- " + str(len(overlapping_original_ranges.keys())) + " original ids were used multiple times"
			# TODO: Uniquify and run on uniquified
			print >> sys.stdout , "--- Please check and correct"
			print >> sys.stdout , "------------------------------"
			print >> sys.stdout , "- Quitting"
			print >> sys.stdout , "------------------------------"
			print_overlapping_info(overlapping_original_ranges, ".error.overlapping_ranges.info")
			sys.exit(2)
			#print >> sys.stdout, "---- Generating intermediate unique files"
			# TODO: uniquify_agp in AGP_lib.py
			# TODO: Save uniquified files
			# TODO: updated variables in script

	##### Find gap to break
	print >> sys.stdout, "- Searching for best breaking position"

	regions_to_break = {}
	regions_to_break_old = {}
	blacklist = {}
	#print sys.stderr , breakpoint_db.keys()
	for seq_id in sorted(breakpoint_db.keys()) :
		blacklist[seq_id] = []
		print >> sys.stdout, "-- Analysing sequence: " + str(seq_id)
		print >> sys.stderr, "## Analysing sequence: " + str(seq_id)
		for break_couple in sorted(breakpoint_db[seq_id]) :
			start , stop = break_couple
			print >> sys.stdout, "--- Breakpoints: " + str(start) + ":" + str(stop)
			print >> sys.stderr, "### Breakpoints: " + str(start) + ":" + str(stop)
			if options.agp :
				#### AGP IS GIVEN - search in AGP junctions for breaking point
				# If the region is between two junction -> blacklist the legacy sequence in-between
				print >> sys.stdout, "---- Searching in given AGP for original sequence junction"
				print >> sys.stderr, "#### Searching in given AGP for original sequence junction"

				start_gap = find_nearest_gap(start , junctions[seq_id] , sequences_len[seq_id] , "both" , max_distance )
				stop_gap = find_nearest_gap(stop , junctions[seq_id] , sequences_len[seq_id] , "both" , max_distance )
				# Test the junctions
				#print >> sys.stderr, start_gap
				#print >> sys.stderr, stop_gap

				if ( not start_gap[0] == "sequence" ) and (not stop_gap[0] == "sequence" ) :
					# The breakpoints hit a junction or the end of the sequence
					# Check it is not the same
					print >> sys.stdout, "----- Breakpoints lead to known junctions"
					print >> sys.stderr, "##### Breakpoints lead to known junctions on " + seq_id + ": [" + str(start_gap[1]) + ":" + str(start_gap[2]) + "] -> [" + str(stop_gap[1]) + ":" + str(stop_gap[2]) + "]"

					if ( start_gap == stop_gap ) or ( start_gap[0] == "extremity" and stop_gap[0] == "extremity" ) :
						# Same junction/extremity or from one extremity to the other:

						if ( start_gap == stop_gap ) :
							print >> sys.stdout, "------ Right and left breakpoints lead to the same junction."
						else :
							print >> sys.stdout, "------ Breakpoints lead to sequence extremities"

						# assign the junction to the nearest, break the sequence for the other
						start_sq_distance = min( (int(start)- start_gap[1])**2 , (int(start) - start_gap[2])**2 )
						stop_sq_distance = min( (int(stop)- start_gap[1])**2 , (int(stop) - start_gap[2])**2 )

						if start_sq_distance <= stop_sq_distance :
							print >> sys.stdout, "------ Searching a gap in sequence compatible with right breakpoint."
							legacy_stop_gap_seq_id, legacy_stop = translate_coords( stop , translation_db[seq_id] )
							print >> sys.stdout, "------ Breakpoint position in legacy sequence: " + seq_id + ":" + str(stop) + "->" + str(legacy_stop_gap_seq_id) + ":" + str(legacy_stop)
							print >> sys.stderr, "###### Breakpoint position in legacy sequence: " + seq_id + ":" + str(stop) + "->" + str(legacy_stop_gap_seq_id) + ":" + str(legacy_stop)
							legacy_stop_gap = find_nearest_gap( legacy_stop , legacy_gap_list[legacy_stop_gap_seq_id] , sequences_legacy_length_db[legacy_stop_gap_seq_id] , "both" , sequences_legacy_length_db[legacy_stop_gap_seq_id] )

						else :
							print >> sys.stdout, "------ Searching a gap in sequence compatible with left breakpoint."
							legacy_start_gap_seq_id , legacy_start = translate_coords( start , translation_db[seq_id] )
							print >> sys.stdout, "------ Breakpoint position in legacy sequence: " + seq_id + ":" + str(start) + "->" + str(legacy_start_gap_seq_id) + ":" + str(legacy_start)
							print >> sys.stderr, "###### Breakpoint position in legacy sequence: " + seq_id + ":" + str(start) + "->" + str(legacy_start_gap_seq_id) + ":" + str(legacy_start)
							legacy_start_gap = find_nearest_gap( legacy_start , legacy_gap_list[legacy_start_gap_seq_id] , sequences_legacy_length_db[legacy_start_gap_seq_id] , "both" , sequences_legacy_length_db[legacy_start_gap_seq_id] )

					else :
						# Different junctions/extremity
						# Blacklist the sequences in between
						print >> sys.stdout, "------ Junctions to break found: [" + str(start_gap[1]) + ":" + str(start_gap[2]) + "] <-> [" + str(stop_gap[1]) + ":" + str(stop_gap[2]) + "]"
						print >> sys.stderr, "###### Junctions to break found: [" + str(start_gap[1]) + ":" + str(start_gap[2]) + "] <-> [" + str(stop_gap[1]) + ":" + str(stop_gap[2]) + "]"
						for component_start in agp_db[seq_id] :
							if (start_gap[1] <= int(component_start) <= stop_gap[2]) and (agp_db[seq_id][component_start][4] == "W") :
								# if the component in the AGP is a sequence and falls withing the breakpoints -> blacklist
								blacklist[seq_id].append(agp_db[seq_id][component_start][5])
						print >> sys.stderr, "###### Blacklisted for " + str(seq_id) + " : " + ",".join([str(x) for x in blacklist[seq_id] ])
						print >> sys.stdout, "------ Blacklisted " + str(len(blacklist[seq_id])) + " sequences"

				else :
					# At least one breakpoint falls inside a sequence.
					# Break the original in the nearest gap
					legacy_start_gap = False
					legacy_stop_gap = False

					if start_gap[0] == "sequence" :
						# No junction is compatible with the breakpoint
						# Search a breakpoint in sequence gaps
						print >> sys.stdout, "----- No junction compatible with left breakpoint. Searching suitable gap in sequence"
						print >> sys.stderr, "##### No junction compatible with left breakpoint. Searching suitable gap in sequence"
						legacy_start_gap_seq_id , legacy_start = translate_coords( start , translation_db[seq_id] )
						print >> sys.stdout, "------ Breakpoint position in legacy sequence: " + seq_id + ":" + str(start) + "->" + str(legacy_start_gap_seq_id) + ":" + str(legacy_start)
						print >> sys.stderr, "###### Breakpoint position in legacy sequence: " + seq_id + ":" + str(start) + "->" + str(legacy_start_gap_seq_id) + ":" + str(legacy_start)
						legacy_start_gap = find_nearest_gap( legacy_start , legacy_gap_list[legacy_start_gap_seq_id] , sequences_legacy_length_db[legacy_start_gap_seq_id] , "both" , sequences_legacy_length_db[legacy_start_gap_seq_id] )
					else :
						# junction is compatible with the breakpoint
						print >> sys.stdout, "----- Junction compatible with left breakpoint found: " + seq_id + ":" + str(start_gap[1]) + "-" + str(start_gap[2])
						print >> sys.stderr, "##### Junction compatible with left breakpoint found: " + seq_id + ":" + str(start_gap[1]) + "-" + str(start_gap[2])

					if stop_gap[0] == "sequence" :
						# No junction is compatible with the breakpoint
						# Search a breakpoint in sequence gaps
						print >> sys.stdout, "----- No junction compatible with right breakpoint. Searching suitable gap in sequence"
						print >> sys.stderr, "##### No junction compatible with right breakpoint. Searching suitable gap in sequence"
						legacy_stop_gap_seq_id, legacy_stop = translate_coords( stop , translation_db[seq_id] )
						print >> sys.stdout, "------ Breakpoint position in legacy sequence: " + seq_id + ":" + str(stop) + "->" + str(legacy_stop_gap_seq_id) + ":" + str(legacy_stop)
						print >> sys.stderr, "###### Breakpoint position in legacy sequence: " + seq_id + ":" + str(stop) + "->" + str(legacy_stop_gap_seq_id) + ":" + str(legacy_stop)
						legacy_stop_gap = find_nearest_gap( legacy_stop , legacy_gap_list[legacy_stop_gap_seq_id] , sequences_legacy_length_db[legacy_stop_gap_seq_id] , "both" , sequences_legacy_length_db[legacy_stop_gap_seq_id] )
					else :
						# junction is compatible with the breakpoint
						print >> sys.stdout, "----- Junction compatible with right breakpoint found: " + seq_id + ":" + str(stop_gap[1]) + "-" + str(stop_gap[2])
						print >> sys.stderr, "##### Junction compatible with right breakpoint found: " + seq_id + ":" + str(stop_gap[1]) + "-" + str(stop_gap[2])

					if legacy_start_gap and legacy_stop_gap and (legacy_start_gap == legacy_stop_gap) :
						print >> sys.stdout, "----- [WARNING] : The same gap is the nearest to both breakpoints and the same gap. Manual inspection suggested"
						print >> sys.stderr, "[WARNING] : Breakpoints - " + str(start) + ":" + str(stop) + " leading to the same gap"
						print >> sys.stderr, "[WARNING] : Nearest gap - " + str(legacy_stop_gap[1]) + ":" + str(legacy_stop_gap[2])
						print >> sys.stderr, "[WARNING] : Junctions in AGP:"
						#for junction in junctions[seq_id] :
						#	print >> sys.stderr, junction
						#print >> sys.stderr, "[WARNING] : Usable gaps in sequence:"
						#for gap in gap_list[seq_id] :
						#	print >> sys.stderr, gap
						print >> sys.stdout, "------ Forcing research to upstream of left and downstream of right"
						print >> sys.stderr, "###### Forcing research to upstream of left and downstream of right"
						nearest = nearest_to_gap( legacy_start , legacy_stop , legacy_start_gap )
						if legacy_start < legacy_stop :
							if nearest == legacy_stop :
								legacy_start_gap = find_nearest_gap( legacy_start , legacy_gap_list[legacy_start_gap_seq_id] , sequences_legacy_length_db[legacy_start_gap_seq_id] , "upstream" , sequences_legacy_length_db[legacy_start_gap_seq_id] )
								print >> sys.stdout, "------- Forced search: upstream gap of left breakpoint:" + str(legacy_start_gap)
								print >> sys.stdout, "------- Forced search: distance of left breakpoint from upstream gap:" + str(int(legacy_start) - int(legacy_start_gap[2]))
								print >> sys.stderr, "####### Forced search: upstream gap of left breakpoint:" + str(legacy_start_gap)
								print >> sys.stderr, "####### Forced search: distance of left breakpoint from upstream gap:" + str(int(legacy_start) - int(legacy_start_gap[2]))
							else :
								legacy_stop_gap = find_nearest_gap( legacy_stop , legacy_gap_list[legacy_stop_gap_seq_id] , sequences_legacy_length_db[legacy_stop_gap_seq_id] , "downstream" , sequences_legacy_length_db[legacy_stop_gap_seq_id] )
								print >> sys.stdout, "------- Forced search: downstream gap of left breakpoint:" + str(legacy_stop_gap)
								print >> sys.stdout, "------- Forced search: distance of right breakpoint from downstream gap" + str(int(legacy_stop_gap[1]) - int(legacy_stop))
								print >> sys.stderr, "####### Forced search: downstream gap of left breakpoint:" + str(legacy_stop_gap)
								print >> sys.stderr, "####### Forced search: distance of right breakpoint from downstream gap" + str(int(legacy_stop_gap[1]) - int(legacy_stop))
						else:
							# Legacy sequence is in reverse direction
							if nearest == legacy_stop :
								legacy_start_gap = find_nearest_gap( legacy_start , legacy_gap_list[legacy_start_gap_seq_id] , sequences_legacy_length_db[legacy_start_gap_seq_id] , "downstream" , sequences_legacy_length_db[legacy_start_gap_seq_id] )
								print >> sys.stdout, "------- Forced search: upstream gap of left breakpoint:" + str(legacy_start_gap)
								print >> sys.stdout, "------- Forced search: distance of left breakpoint from upstream gap:" + str(int(legacy_start) - int(legacy_start_gap[2]))
								print >> sys.stderr, "####### Forced search: upstream gap of left breakpoint:" + str(legacy_start_gap)
								print >> sys.stderr, "####### Forced search: distance of left breakpoint from upstream gap:" + str(int(legacy_start) - int(legacy_start_gap[2]))
							else :
								legacy_stop_gap = find_nearest_gap( legacy_stop , legacy_gap_list[legacy_stop_gap_seq_id] , sequences_legacy_length_db[legacy_stop_gap_seq_id] , "upstream" , sequences_legacy_length_db[legacy_stop_gap_seq_id] )
								print >> sys.stdout, "------- Forced search: downstream gap of left breakpoint:" + str(legacy_stop_gap)
								print >> sys.stdout, "------- Forced search: distance of right breakpoint from downstream gap" + str(int(legacy_stop_gap[1]) - int(legacy_stop))
								print >> sys.stderr, "####### Forced search: downstream gap of left breakpoint:" + str(legacy_stop_gap)
								print >> sys.stderr, "####### Forced search: distance of right breakpoint from downstream gap" + str(int(legacy_stop_gap[1]) - int(legacy_stop))

					if legacy_start_gap :
						if legacy_start_gap_seq_id not in regions_to_break_old :
							regions_to_break_old[legacy_start_gap_seq_id] = []
						regions_to_break_old[legacy_start_gap_seq_id].append( legacy_start_gap )
					if legacy_stop_gap :
						if legacy_stop_gap_seq_id not in regions_to_break_old :
							regions_to_break_old[legacy_stop_gap_seq_id] = []
						regions_to_break_old[legacy_stop_gap_seq_id].append( legacy_stop_gap )

				if seq_id not in regions_to_break :
					regions_to_break[seq_id] = []
				regions_to_break[seq_id].append(start_gap)
				regions_to_break[seq_id].append(stop_gap)

			else :
				#### AGP IS NOT GIVEN
				## Only in-sequence-gap search
				print >> sys.stdout, "---- No AGP given, searching in the sequence encoded gaps"
				stop_gap_new = find_nearest_gap( stop ,  gap_list[seq_id] , sequences_len[seq_id] , "both" , sequences_len[seq_id] )
				start_gap_new = find_nearest_gap( start , gap_list[seq_id] , sequences_len[seq_id] , "both" , sequences_len[seq_id] )
				if seq_id not in regions_to_break :
					regions_to_break[seq_id] = []
				regions_to_break[seq_id].append( start_gap_new )
				regions_to_break[seq_id].append( stop_gap_new )

	broken_gaps_file_name = options.out + ".broken_gaps.given_sequences.list.txt"
	broken_gaps_file = open(broken_gaps_file_name , 'w')
	print >> sys.stderr, "# Final list of gap break"
	print >> sys.stderr, "## Given sequences"
	for seq_id in sorted(regions_to_break.keys()) :
		# uniquify the list of gap to break
		old_list = regions_to_break[seq_id]
		new_list = list(set(tuple(i) for i in old_list))
		regions_to_break[seq_id] = new_list
		print >> broken_gaps_file, '>' + seq_id
		for element in regions_to_break[seq_id] :
			print >> broken_gaps_file, "\t".join([ str(x) for x in element ])
	#print >> sys.stderr, regions_to_break
	broken_gaps_file.close()

	if not regions_to_break_old == {} :
		print >> sys.stderr, "## Legacy sequences"
		broken_gaps_legacy_file_name = options.out + ".broken_gaps.legacy_sequences.list.txt"
		broken_gaps_legacy_file = open(broken_gaps_legacy_file_name , 'w')
		for legacy_stop_gap_seq_id in sorted(regions_to_break_old.keys()) :
			# uniquify the list of gap to break
			old_list = regions_to_break_old[legacy_stop_gap_seq_id]
			new_list = list(set(tuple(i) for i in old_list))
			regions_to_break_old[legacy_stop_gap_seq_id] = new_list
			print >> broken_gaps_legacy_file, ">" + legacy_stop_gap_seq_id
			for element in regions_to_break_old[legacy_stop_gap_seq_id] :
				print >> broken_gaps_legacy_file, "\t".join([ str(x) for x in element ])
		#print >> sys.stdout, regions_to_break_old
		broken_gaps_legacy_file.close()


	##### Generate new AGP using breakpoints
	print >> sys.stdout, "- Updating sequences with novel breaking points"
	print >> sys.stderr, "# Updating sequences with novel breaking points"

	if options.agp :
		print >> sys.stdout, "-- Splitting legacy sequences by gap breaking"
		new_legacy_sequences , new_legacy_agp_db = get_broken_agp( sequences_legacy, regions_to_break_old , options.prefix)
		new_legacy_agp_file = write_agp( new_legacy_agp_db , options.out + ".new_legacy.agp" )
		#print >> sys.stderr, new_legacy_agp_db
		#legacy_to_new_legacy_agp_db = invert_agp( new_legacy_agp_db )
		#legacy_to_new_legacy_agp_file = write_agp( legacy_to_new_legacy_agp_db , options.out + ".legacy_to_new_legacy.agp" )
		#print >> sys.stdout, "-- Updating sequences to novel legacy sequences relationship"
		#legacy_to_new_agp_db = agp_translate_agp( agp_db , legacy_to_new_legacy_agp_db )
		#legacy_to_new_agp_file = write_agp( legacy_to_new_agp_db , options.out + ".legacy_to_new.agp" )
	else :
		print >> sys.stdout, "-- Splitting sequences by gap breaking"
		#print >> sys.stderr, regions_to_break
		new_sequences , new_agp_db = get_broken_agp(sequences, regions_to_break , options.prefix)
		new_agp_file = write_agp( new_agp_db , options.out + ".new.agp" )

	#### Print new sequences
	print >> sys.stdout, "-- Printing new sequences"
	if options.agp :
		write_fasta_from_db( new_legacy_sequences , options.out + ".new_legacy.fasta" )
	else :
		write_fasta_from_db( new_sequences , options.out + ".new.fasta" )

	#### Print blacklist
	print >> sys.stdout, "-- Printing file of blacklisted sequences"
	blacklist_file_name = options.out + ".blacklist.txt"
	blacklist_file = open(blacklist_file_name , 'w')
	for seq_id in sorted(blacklist.keys()) :
		info = ",".join(str(x) for x in blacklist[seq_id])
		if not info.rstrip().lstrip() == "" :
			blacklist_line = seq_id + "\t" + info
			print >> blacklist_file , blacklist_line
	blacklist_file.close()

	#### update and print new annotation files (if set)
	if options.gff3 :
		print >> sys.stdout, "-- Updating annotation cooridnates"
		print >> sys.stdout, "--- Converting coordinates to updated sequences"
		
		if options.agp :
			print >> sys.stdout, "--- Converting coordinates to updated legacy sequences"
			new_legacy_translate_db = translate_from_AGP_whole_genome(new_legacy_agp_db)
			#print >> sys.stderr, "##"
			#print >> sys.stderr, new_legacy_agp_db
			#print >> sys.stderr, "##"
			#print >> sys.stderr, legacy_annotation_gff3
			new_legacy_annotation_gff3 = translate_gff3( legacy_annotation_gff3 , new_legacy_translate_db , options.out + ".broken_loci.annotation.legacy_updated.txt" )
			print >> sys.stdout, "---- Generating legacy GFF3 file for updated legacy sequences"
			write_gff3(new_legacy_annotation_gff3 , options.out + ".new_legacy.annotation.gff3" , get_length_from_fasta_db( new_legacy_sequences ) )
		else:
			new_translate_db = translate_from_AGP_whole_genome(new_agp_db)
			new_annotation_gff3 = translate_gff3( annotation_gff3 , new_translate_db , options.out + ".broken_loci.annotation.txt" )
			print >> sys.stdout, "---- Generating legacy GFF3 file for updated sequences"
			write_gff3(new_annotation_gff3 , options.out + ".new.annotation.gff3" , get_length_from_fasta_db( new_sequences ) )

	# Convert BED file if given
	if options.bed :
		bed_regions = read_bed_sorted_list( options.bed )
		print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] = Translating coordinates of features in bed the file'
		print >> sys.stderr , '# Translating coordinates of features in bed the file'
		if options.agp:
			new_legacy_bed = translate_bed_sorted_list( bed_regions ,  new_legacy_agp_db )
			out_bed_file_name = options.out + ".new_legacy.bed"
			out_bed_file = open(out_bed_file_name , 'w')
			for line in sorted(new_legacy_bed) :
				print >> out_bed_file, "\t".join([str(x) for x in line])
			out_bed_file.close()
		else:
			new_bed = translate_bed_sorted_list( bed_regions ,  new_agp_db )
			out_bed_file_name = options.out + ".new.bed"
			out_bed_file = open(out_bed_file_name , 'w')
			for line in sorted(new_bed) :
				print >> out_bed_file, "\t".join([str(x) for x in line])
			out_bed_file.close()

	##### Finished

	print >> sys.stdout , "------------------------------"
	print >> sys.stdout , "- Done"
	print >> sys.stdout , "------------------------------"
	print >> sys.stderr , "##############################"
	print >> sys.stderr , "# Done"
	print >> sys.stderr , "##############################"


if __name__ == '__main__':
	main()
