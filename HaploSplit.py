#!/usr/bin/env python

import argparse
import itertools
from lib_files.HaploFunct import *
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *
from lib_files.GFF_lib import *

gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


def main() :

	###### Options and help ######

	parser = argparse.ArgumentParser()

	#### Input sequences
	parser.add_argument("-i", "--input" , dest="query",
					help="FASTA file of the input sequences", metavar="query.fasta [Required]")
	parser.add_argument("--GFF3", dest="gff3",
					help="Annotation in GFF3 format", metavar="genes.gff3")
	parser.add_argument("-a" , "--input_agp" , dest="input_agp",
					help="AGP file of query sequences composition" , metavar="query.agp" )
	parser.add_argument("--input_groups" , dest="input_groups",
					help="Tab separated values file of group association of input sequences. Two columns reporting the sequences ID and the association groups id. Sequence can be associated to multiple groups with multiple rows. If not provided, [-e|--exclusion] and [-k|--known] will be used for grouping. To force the usage of alignment hits pass as value \"0\" (zero)"  , metavar="input_groups.tsv")
	parser.add_argument("--legacy_groups" , dest="legacy_groups",
					help="Tab separated values file of group association of components present the input sequences. Two columns reporting the component ID and the association group it belongs to. Sequence can be associated to multiple groups with multiple rows.", metavar="legacy_groups.tsv")

	#### Markers
	parser.add_argument("-m", "--markers", dest="marker_map",
					help="Table of markers order when sorted genomic position. Tab separated file with 3 columns: 1) Chromosome ID, 2) Marker sorting position, 3) Marker ID. If a guide genome is used for collinearity search (-g/--guide), chromosomes ID must match the the files", metavar="markers_list")
	parser.add_argument("-n" , "--map", dest="markers_hits",
					help="Map of markers position on input sequences. BED-like tab separated file with 3 columns 1) Sequence ID, 2) Match start (0-based), 3) Match start (1-based), 4) Marker ID", metavar="map.tsv")
	parser.add_argument("-f" , "--filter", dest="filter_hits", default=False , action="store_true",
					help="If set, remove from the analysis the sequences with intra-sequence markers duplication. [Default: filter markers only, keep all sequences]")
	parser.add_argument("--extended", dest="extended_region", default=False , action="store_true",
					help="Extend the association of a sequence to all markers belonging to the same chromosome regardless of sequence orientation. [Default: only the markers in the expected orientation are considered for tiling analysis]")


	#### Guide genome mapping
	parser.add_argument("-g", "--guide", dest="reference",
					help="FASTA file with reference sequence", metavar="guide.fasta")
	parser.add_argument("-l" , "--local_alignment", dest="hits", default="local.paf",
					help="PAF mapping file with local alignment hits of query sequences on guide sequences. If --align is set, the file will be created or overwritten, otherwise it will be used as input", metavar="local.paf")

	# Output generation
	parser.add_argument("-o", "--out", dest="out", default="out",
					help="Output files name prefix [default: out]", metavar="NAME")
	parser.add_argument("-p", "--prefix", dest="prefix", default="NEW",
					help="Prefix for output sequences IDs [default: NEW]", metavar="PREFIX")
	parser.add_argument("--agp", dest="agp", default=True, action="store_false",
					help="Do not export structures in AGP format")
	parser.add_argument("--fasta", dest="sequence", default=True, action="store_false",
					help="Do not export sequences in FASTA format")
	parser.add_argument("--gapsize", dest="gap", default="1000",
					help="Minimum gap size placeholder for AGP and FASTA (in bp) [default: 1,000bp]" , metavar="N")
	parser.add_argument("--concatenate", dest="conc", default="",
					help="Set the gap size used to concatenate unplaced sequences (in bp) [default: do not concatenated]", metavar="N")

	# Hit tiling control
	parser.add_argument("--distance1", dest="distance1", default="2000000",
					help="Set the maximum distance (bp) between sequence projections on the guide genome to be considered as adjacent of the same tiling path alignment for Hap1 [default: 2,000,000bp]", metavar="N")
	parser.add_argument("--distance2", dest="distance2", default="4000000",
					help="Set the maximum distance (bp) between sequence projections on the guide genome to be considered as adjacent of the same tiling path alignment for Hap2 [default: 4,000,000bp]",metavar="N" )

	# Mapping control
	parser.add_argument("--align" , dest="map", action="store_true",
					help="Do run mapping [overwrite input file.paf if existing]")
	parser.add_argument("-t", "--tool", dest="mapping", default=" --cs -x asm20 -r 1000",
					help="Mapping command for minimap2", metavar="\" --cs -x asm20 -r 1000\"")
	parser.add_argument("-c", "--cores", dest="cores", default=4,
					help="Cores used in mapping process [default: 4]", metavar="N")
	parser.add_argument("--hitgap" , dest="hitgap", default="100000",
					help="Allowed  maximum distance (in bp) between local alignment hits to be merged as part of the same global alignment. [default: 100,000bp]", metavar="N")

	# Tiling path control
	parser.add_argument("-e", "--exclusion", dest="exclusion",
					help="Tab separated file of sequences pairs to prevent being placed in the same haplotype", metavar="exclusion.tsv")
	parser.add_argument("-k", "--known", dest="known",
					help="Tab separated file of sequences known to be in the same haplotype. Two columns reporting the sequence ID and the associated group ID", metavar="known.tsv")
	parser.add_argument("--allow_rearrangements", dest="rearrangements", default=False, action="store_true",
					#help="Restrict constrains on sequence to a the given chromosome, allowing rearrangements on different chromosomes [default: Do not allow]")
					help=argparse.SUPPRESS )
	parser.add_argument("--alternative_groups", dest="alternative",
					help="Tab separated file with two columns reporting ids (comma separated lists) of sequences that should be used in alternative haplotypes.", metavar="alternative_groups.tsv")
	parser.add_argument("-1", "--path1", dest="use1",
					#help="Use ONLY the following list of (stranded) contigs to create the first haplotype. Tab separated file with Target_ID followed by comma separated list of contigs (contig1|+,contig2|-,[...],contigN|+)", metavar="1st.txt")
					help=argparse.SUPPRESS )
	parser.add_argument("-2", "--path2", dest="use2",
					#help="Use ONLY the following list of (stranded) contigs to create the second haplotype. Tab separated file with Target_ID followed by comma separated list of contigs (contig1|+,contig2|-,[...],contigN|+)", metavar="2nd.txt")
					help=argparse.SUPPRESS )
	parser.add_argument("--required_as_path", dest="forced_as_path" , default=False, action="store_true",
					help="Manage lists in required sequences form \"--R1\" and \"--R2\" as tiling paths")
	parser.add_argument("--R1", dest="Require1",
					help="Require to use the following list of (stranded) contigs in the first haplotype. Tab separated file with Target_ID followed by comma separated list of contigs (contig1|+,contig2|-,[...],contigN|+)", metavar="1st.txt")
	parser.add_argument("--F1" , "--forcemarkers1" , dest="force_direction1", default=False, action="store_true",
					help="Force the use of markers according to chromosome and direction reported in --R1")
	parser.add_argument("--min1", dest="minR1", default="0" ,
					help="Minimum length of sequence allowed to be a requirement for the first haplotype (default: 0)", metavar="N")
	parser.add_argument("--R2", dest="Require2",
					help="Require to use the following list of (stranded) contigs in the second haplotype. Tab separated file with Target_ID followed by comma separated list of contigs (contig1|+,contig2|-,[...],contigN|+)", metavar="2nd.txt")
	parser.add_argument("--F2" , "--forcemarkers2" , dest="force_direction2", default=False, action="store_true",
					help="Force the use of markers according to chromosome and direction reported in --R2")
	parser.add_argument("--min2", dest="minR2", default="0" ,
					help="Minimum length of sequence allowed to be a requirement for the second haplotype (default: 0)", metavar="N")
	parser.add_argument("--conflict_resolution", dest="conflict_resolution" , default=2,
					help=argparse.SUPPRESS )
	parser.add_argument("--B1", "--blacklist1", dest="Blacklist1",
					help="Blacklisted (stranded) contigs NOT to be used in the first haplotype. Tab separated file with Target_ID followed by comma separated list of contigs (contig1,contig2,[...],contigN)", metavar="2nd.txt")
	parser.add_argument("--B2", "--blacklist2",  dest="Blacklist2",
					help="Blacklisted (stranded) contigs NOT to be used in the second haplotype. Tab separated file with Target_ID followed by comma separated list of contigs (contig1,contig2,[...],contigN)", metavar="2nd.txt")

	# processing steps control
	parser.add_argument("--N2", dest="No2", action="store_true",
					help="Don't run the search for the 2nd path")
	parser.add_argument("--reuse_intermediate" , dest="reuse_intermediate", default=False, action="store_true",
					help="Reuse the sequences and the intermediate results of a previous analysis and rerun just the plot" )
	parser.add_argument("-v", "--dry", dest="dry", action="store_true",
					help="Dry run: check uniqueness of marker in input sequences and quit")

	# QC control
	parser.add_argument("--skip_chimeric_qc", dest="skip_qc", default=False , action="store_true",
					help="Do not perform QC report on potentially chimeric input sequences")
	parser.add_argument("--avoid_rejected_qc" , dest="avoidrejectedqc", default=False, action="store_true",
					help="Avoid running the quality control of the sequences assigned to a chromosome but not placed in the reconstructed pseudomolecule [Default: run]")
	parser.add_argument("--only_markers" , dest="only_markers", default=False, action="store_true",
					help="Limit the quality control of the rejected sequences only to those with markers [Default: Do all]")
	parser.add_argument("--haplodup" , dest="haplodup", default=False, action="store_true",
					help="Run HaploDup on assembly results. If --GFF3 is set, GMAP will be used to generate interactive plots for deduplication analysis, otherwise only dotplots will be delivered")
	parser.add_argument("--debug" , dest="debug", default=False, action="store_true",
					help=argparse.SUPPRESS )
	parser.add_argument("--disable_marker_ploidy_check", dest="disable_marker_ploidy_check" , default=False, action="store_true",
					help=argparse.SUPPRESS )

	scriptDirectory = os.path.dirname(os.path.realpath(__file__)) + "/support_scripts"

	print >> sys.stdout, "Running HaploSplit tool from HaploSync version " + get_version()
	print >> sys.stdout, "To reproduce this run use the following command: " + " ".join( pipes.quote(x) for x in sys.argv)
	print >> sys.stdout, "----"

	# Sanity Check

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)

	options = parser.parse_args()

	if not options.reference and not options.marker_map and not options.markers_hits:
		print >> sys.stderr , "[ERROR] Information for phasing pseudomolecules was given. -g/--guide and/or -m/--markers with -n/--map must be given"
		parser.print_help()
		sys.exit(1)

	if (not options.marker_map and options.markers_hits) or (options.marker_map and not options.markers_hits):
		print >> sys.stderr , "[ERROR] Marker information is partial. -m/--markers and -n/--map must be given concurrently"
		parser.print_help()
		sys.exit(1)

	if not options.query :
		print >> sys.stderr , "[ERROR] Query genome FASTA file missing"
		parser.print_help()
		sys.exit(1)

	###### Main ######

	paths = set_paths(os.path.join(sys.path[0], 'HaploSync.conf.toml'))
	nucmer_path = paths["nucmer"]
	showcoords_path = paths["show-coords"]

	ref_to_hap1 = {}
	hap1_to_ref = {}
	hap1_ids = []
	ref_ids = []
	all_seq_length = {}
	if not options.No2 :
		hap2_ids = []
		ref_to_hap2 = {}
		hap2_to_ref = {}
		hap2_to_hap1 = {}
		hap1_to_hap2 = {}

	### Read query and compute lengths
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Reading query sequences"
	query_fasta_db = read_fasta(options.query)
	query_len = {}
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Calculating query sequences lengths"
	for query_name in query_fasta_db :
		query_len[query_name] = int(len(query_fasta_db[query_name]))
		all_seq_length[query_name] = int(len(query_fasta_db[query_name]))
		all_seq_length[query_name + "|+"] = int(len(query_fasta_db[query_name]))
		all_seq_length[query_name + "|-"] = int(len(query_fasta_db[query_name]))
		all_seq_length[query_name + "|."] = int(len(query_fasta_db[query_name]))
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Number of query sequences: " + str(len(query_len.keys()))

	### Read input AGP if given
	if options.input_agp :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Reading query structure from AGP file"
		query_agp_db = read_agp(options.input_agp)
	else :
		query_agp_db = {}

	if options.gff3 :
		### Read GFF3 input
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Reading annotation from GFF3"
		annotation_gff3, mRNA_db = read_gff3(options.gff3)
	else :
		annotation_gff3 = {}
		mRNA_db = {}

	reference_sequence_list = []

	### Read marker map and markers hits
	marker_hits_by_seq = {}
	marker_hits_by_id = {}
	unknown_markers = []

	if options.marker_map and options.markers_hits :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Reading markers map"
		marker_map_by_seq = {}
		marker_map_by_id = {}
		# Read map
		for line in open( options.marker_map ) :
			if line == "" :
				continue
			if line[0] == "#" :
				continue
			try :
				seq_id , pos , marker_id = line.rstrip().split("\t")
			except :
				print >> sys.stdout, "Error in markers map file " + options.marker_map + ", line do not match file format: " + line.rstrip()
				print >> sys.stdout, "Error in markers map file " + options.marker_map + ", line do not match file format: " + line.rstrip()
				sys.exit(1)
			else :
				if seq_id not in marker_map_by_seq :
					reference_sequence_list.append(seq_id)
					marker_map_by_seq[seq_id] = []
				marker_map_by_seq[seq_id].append([ int(pos) , marker_id ] )
				marker_map_by_id[marker_id] = [ seq_id , int(pos) , marker_id ]

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Reading markers hits"

		# Read marker position
		for line in open( options.markers_hits ) :
			if line == "" :
				continue
			if line[0] == "#" :
				continue
			try :
				seq_id , pos_1 , pos_2 , marker_id = line.rstrip().split("\t")
			except :
				print >> sys.stdout, "Error in markers hits file " + options.markers_hits + ", line do not match file format."
				print >> sys.stderr, "Error in markers hits file " + options.markers_hits + ", line do not match file format."
				print >> sys.stderr, "Line format Expected: Seq_ID\tStart\tStop\tMarker_ID"
				print >> sys.stderr, "Line format found: " + line.rstrip()
				sys.exit(1)
			else :
				if marker_id in marker_map_by_id :
					if seq_id not in marker_hits_by_seq :
						marker_hits_by_seq[seq_id] = []
					marker_pos = marker_map_by_id[marker_id][1]
					marker_chr = marker_map_by_id[marker_id][0]
					start = min( int(pos_1) , int(pos_2) )
					stop = max( int(pos_1) , int(pos_2) )
					marker_hits_by_seq[seq_id].append([ int(start) , int(stop) , marker_id , marker_chr , int(marker_pos) ] )
					if marker_id not in marker_hits_by_id :
						marker_hits_by_id[marker_id] = []
					marker_hits_by_id[marker_id].append( [ seq_id , int(start) , int(stop) ] )
				else :
					unknown_markers.append(marker_id)
		if not unknown_markers == [] :
			unknown_markers_file_name = options.out + ".unknown_markers.txt"
			unknown_markers_file = open(unknown_markers_file_name ,'w')
			unknown_markers = list(set(unknown_markers))
			for unknown_marker_id in sorted(unknown_markers) :
				print >> unknown_markers_file , unknown_marker_id
			unknown_markers_file.close()
			print >> sys.stderr, "[WARNING] Markers alignment on the query sequences reports " + str(len(unknown_markers)) + " unknown marker ids from the genetic map. Unknown markers will be exclude, their ids are reported in " + unknown_markers_file_name + " file."

	### Read reference, if given, and compute lengths
	if options.reference:
		add_seq_id=False
		if reference_sequence_list == [] :
			add_seq_id=True
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Reading reference sequences"
		reference = read_fasta(options.reference)
		reference_len = {}
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Calculating reference sequences lengths"
		# check if map sequences id are all represented -> issue warning
		for ref_name in reference_sequence_list :
			if ref_name not in reference :
				print >> sys.stderr , "[WARNING] Reference " + ref_name + " is present in the markers map but missing from guide genome fasta"
		for ref_name in reference :
			if add_seq_id :
				reference_sequence_list.append(ref_name)
			else :
				if ref_name not in reference_sequence_list :
					# Sequence name in fasta file do not match any map -> issue warning
					print >> sys.stderr , "[WARNING] Guide sequence " + ref_name + " is present in the guide genome but absent from markers map. Sequence will be excluded from the analysis"
			reference_len[ref_name] = int(len(reference[ref_name]))
			all_seq_length[ref_name] = int(len(reference[ref_name]))
			all_seq_length[ref_name + "|+"] = int(len(reference[ref_name]))
			all_seq_length[ref_name + "|-"] = int(len(reference[ref_name]))
			all_seq_length[ref_name + "|."] = int(len(reference[ref_name]))
			ref_ids.append(ref_name)
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Number of reference sequences: " + str(len(reference_len.keys()))

	### For each chromosome, Create graph
	best_1_paths = {}
	best_1_paths_edges = {}
	best_2_paths = {}
	best_2_paths_edges ={}
	all_used = []
	used_by_chr_hap1 = {}
	used_by_chr_hap2 = {}
	unused_by_chr = {}


	# Read mutually exclusive pairs
	# 	To avoid tiling path with mutually exclusive sequences:
	#   1) remove direct edges between all unwanted pairs
	#   2) Iterate:
	#		a) get the longest tiling path
	#		b) find connected pairs (check if both are in tiling path)
	#		c) remove the shortest sequence fo all wrong from the network
	unwanted_pairs = {}
	if options.exclusion :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Reading mutually exclusive sequences pairs"
		for line in open(options.exclusion) :
			if ( line == "" ) or ( line[0] == "#") :
				continue
			else:
				try :
					seq_1_id , seq_2_id = line.rstrip().split("\t")
				except :
					print >> sys.stdout , "[ERROR] Line in exclude pair file (" + options.exclusion + ") does not contain two ids"
					print >> sys.stderr , "[ERROR] Line in exclude pair file (" + options.exclusion + ") does not contain two ids"
					print >> sys.stderr , "[ERROR] Line expected: Seq_1_id\tSeq_2_id"
					print >> sys.stderr , "[ERROR] Line found: " + line.rstrip()
					sys.exit(1)
				if (seq_1_id + "|+") not in unwanted_pairs :
					unwanted_pairs[seq_1_id + "|+"] = []
					unwanted_pairs[seq_1_id + "|-"] = []
					unwanted_pairs[seq_1_id + "|."] = []
				unwanted_pairs[seq_1_id + "|+"] += three_orientation_list(seq_2_id , False)
				unwanted_pairs[seq_1_id + "|-"] += three_orientation_list(seq_2_id , False)
				unwanted_pairs[seq_1_id + "|."] += three_orientation_list(seq_2_id , False)
				if (seq_2_id + "|+") not in unwanted_pairs :
					unwanted_pairs[seq_2_id + "|+"] = []
					unwanted_pairs[seq_2_id + "|-"] = []
					unwanted_pairs[seq_2_id + "|."] = []
				unwanted_pairs[seq_2_id + "|+"] += three_orientation_list(seq_1_id , False)
				unwanted_pairs[seq_2_id + "|-"] += three_orientation_list(seq_1_id , False)
				unwanted_pairs[seq_2_id + "|."] += three_orientation_list(seq_1_id , False)
	#unwanted_pairs_file = open(options.out + ".unwanted_pairs.txt", 'w')
	#json.dump(unwanted_pairs , unwanted_pairs_file , indent=4)
	#unwanted_pairs_file.close()

	# Known relationships
	known_groups = {}
	known_groups_by_seqid = {}
	if options.known :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Reading known groups of sequences in the same haplotype"
		for line in open(options.known) :
			if ( line == "" ) or ( line[0] == "#") :
				continue
			else:
				try :
					seq_id , group_id = line.rstrip().split("\t")
				except :
					print >> sys.stdout , "[ERROR] Line in exclude pair file (" + options.known + ") does not contain two ids"
					print >> sys.stderr , "[ERROR] Line in exclude pair file (" + options.known + ") does not contain two ids"
					print >> sys.stderr , "[ERROR] Line expected: Seq_id\tGroup_id"
					print >> sys.stderr , "[ERROR] Line found: " + line.rstrip()
					sys.exit(1)
				if group_id not in known_groups :
					known_groups[group_id] = []
				known_groups[group_id].append(seq_id)
		#print >> sys.stderr , known_groups.keys()
		for group_id in known_groups.keys() :
			for seq_id in known_groups[group_id] :
				if seq_id not in known_groups_by_seqid :
					known_groups_by_seqid[seq_id] = []
					known_groups_by_seqid[seq_id + "|+" ] = []
					known_groups_by_seqid[seq_id + "|-" ] = []
					known_groups_by_seqid[seq_id + "|." ] = []
				for seq_2_id in known_groups[group_id] :
					known_groups_by_seqid[seq_id].append(seq_2_id)
					known_groups_by_seqid[seq_id] += three_orientation_list(seq_2_id , False)
					known_groups_by_seqid[seq_id + "|+" ].append(seq_2_id)
					known_groups_by_seqid[seq_id + "|+" ] += three_orientation_list(seq_2_id , False)
					known_groups_by_seqid[seq_id + "|-" ].append(seq_2_id)
					known_groups_by_seqid[seq_id + "|-" ] += three_orientation_list(seq_2_id , False)
					known_groups_by_seqid[seq_id + "|." ].append(seq_2_id)
					known_groups_by_seqid[seq_id + "|." ] += three_orientation_list(seq_2_id , False)
		#json.dump(known_groups_by_seqid , open("known_groups.by_seq_id.txt", 'w') , indent=4)
		#json.dump(known_groups , open("known_groups.by_cluster.txt", 'w') , indent=4)

	alternative_sequences = {}
	if options.alternative :
		for line in open( options.alternative ) :
			if line[0] == "#" or line.rstrip() == "" :
				continue
			else:
				group_1 , group_2 = line.rstrip().split("\t")
				#print >> sys.stderr , "group_1 and group_2"
				group_1_list = group_1.split(",")
				group_2_list = group_2.split(",")
				#print >> sys.stderr , group_1_list
				#print >> sys.stderr , group_2_list
				for seq_id in group_1_list :
					#print >> sys.stderr , seq_id
					if seq_id not in alternative_sequences :
						alternative_sequences[seq_id] = []
					alternative_sequences[seq_id] += group_2_list
					alternative_sequences[seq_id] = list(set(alternative_sequences[seq_id]))

				for seq_id in group_2_list :
					if seq_id not in alternative_sequences :
						alternative_sequences[seq_id] = []
					alternative_sequences[seq_id] += group_1_list
					alternative_sequences[seq_id] = list(set(alternative_sequences[seq_id]))

	# Force path if necessary
	fixed_path_1 = defaultdict(list)
	fixed_path_2 = defaultdict(list)

	if options.use1 or options.use2 :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Reading given paths"

	if options.use1 :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Sequences used in 1st path"
		for line in open(options.use1) :
			try :
				Target_sequence , loaded_path = line.rstrip().split("\t")
				prompted_list = loaded_path.split(",")
				best_1_paths_edges[Target_sequence] = []
				best_1_paths[Target_sequence] = []
				fixed_path_1[Target_sequence] = []
				for id in prompted_list :
					best_1_paths[Target_sequence].append([id,0,0,0,0])
					fixed_path_1[Target_sequence].append(id)
					best_1_paths_edges[Target_sequence].append(id)
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Seq " + str(Target_sequence) + " : " + ",".join(best_1_paths_edges[Target_sequence])
			except:
				pass

	if options.use2 :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Sequences used in 2nd path"
		for line in open(options.use2) :
			try :
				Target_sequence , loaded_path = line.rstrip().split("\t")
				prompted_list = loaded_path.split(",")
				best_2_paths_edges[Target_sequence] = []
				best_2_paths[Target_sequence] = []
				fixed_path_1[Target_sequence] = []
				for id in prompted_list :
					best_2_paths[Target_sequence].append([id,0,0,0,0])
					fixed_path_2[Target_sequence].append(id)
					best_2_paths_edges[Target_sequence].append(id)
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Seq " + str(Target_sequence) + " : " + ",".join(best_2_paths_edges[Target_sequence])
			except :
				pass


	# Create required lists
	forced_list_1 = defaultdict(list)
	forced_list_2 = defaultdict(list)
	discarded=[]
	forced_list_1_id_oriented=[]
	forced_list_1_id=[]
	forced_list_2_id_oriented=[]
	forced_list_2_id=[]
	# Create Blacklists
	blacklist_1 = defaultdict(list)
	blacklist_2 = defaultdict(list)

	if options.Require1 or options.Require2 :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Reading required sequence lists"

	if options.Require1 :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Sequence required in 1st path"
		for line in open(options.Require1,"r") :
			try :
				Target_sequence , loaded_path = line.rstrip().split("\t")
				prompted_list = loaded_path.split(",")
				forced_list_1[Target_sequence] = []
				for id in prompted_list :
					id_name = id[:-2]
					id_length = query_len[id_name]
					if id_length >= int(options.minR1) :
						forced_list_1[Target_sequence].append(id)
						forced_list_1_id_oriented.append(id)
						forced_list_1_id.append(id_name)
					else :
						discarded.append(id_name)
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Seq " + str(Target_sequence) + ": " + ",".join(forced_list_1[Target_sequence])
			except :
				pass


	if options.Require2 :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Sequence required in 2nd path"
		for line in open(options.Require2,"r") :
			try :
				Target_sequence , loaded_path = line.rstrip().split("\t")
				prompted_list = loaded_path.split(",")
				forced_list_2[Target_sequence] = []
				for id in prompted_list :
					id_name = id[:-2]
					id_length = query_len[id_name]
					if id_length >= int(options.minR2) :
						forced_list_2[Target_sequence].append(id)
						forced_list_2_id_oriented.append(id)
						forced_list_2_id.append(id_name)
					else :
						discarded.append(id_name)
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Seq " + str(Target_sequence) + ": " + ",".join(forced_list_2[Target_sequence])
			except :
				pass
		doubled = []
		for seq_id in forced_list_2_id :
			if seq_id in forced_list_1_id :
				doubled.append(seq_id)
		if not doubled == [] :
			print >> sys.stdout, '[ERROR] Input error. Same sequence(s) required for both haplotypes'
			print >> sys.stderr, '[ERROR] Input error. Same sequence(s) required for both haplotypes'
			print >> sys.stderr, '[ERROR] duplicated sequences IDs:'
			print >> sys.stderr, ", ".join([str(x) for x in doubled])
			sys.exit(1)

	if not discarded == [] :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Discarded required sequences because of short length (hap1 <" + str(options.minR1) +" or hap2 <" + str(options.minR2) + ")"
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === " + ",".join(discarded)

	# Generate blacklist
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Blacklists"

	for ref_name in reference_sequence_list :
		blacklist_1[ref_name]=[]
		blacklist_2[ref_name]=[]

	if options.Blacklist1 :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Rejected for 1st path"
		for line in open(options.Blacklist1,"r") :
			try :
				Tid , Qids = line.rstrip().split("\t")
				if Tid in blacklist_1 :
					for name in Qids.split(",") :
						blacklist_1[Tid] += three_orientation_list(name)
					print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Seq " + str(Tid) + " : " + ",".join(blacklist_1[Tid])
			except:
				pass

	for key_in in sorted(blacklist_1.keys()) :
		# Blacklist sequences required somewhere else
		Required_somewhere_else = []

		for key_out in sorted(blacklist_1.keys()) :
			if not key_in == key_out :
				Required_somewhere_else += forced_list_1[key_out]
				Required_somewhere_else += fixed_path_1[key_out]
			Required_somewhere_else += forced_list_2[key_out]
			Required_somewhere_else += fixed_path_2[key_out]

		for name in list(set(Required_somewhere_else)) :
			blacklist_1[key_in] += three_orientation_list(name)

		print >> sys.stderr , "### To reject for Path1@" + str(key_in) + " : " + ",".join(blacklist_1[key_in])


	if options.Blacklist2 :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Rejected for 2nd path"
		for line in open(options.Blacklist2,"r") :
			try :
				Tid , Qids = line.rstrip().split("\t")
				if Tid in blacklist_2 :
					for name in Qids.split(",") :
						blacklist_2[Tid] += three_orientation_list(name)
					print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Seq " + str(Tid) + " : " + ",".join(blacklist_2[Tid])
			except:
				pass


	for key_in in sorted(blacklist_2.keys()) :
		# Blacklist sequences required somewhere else
		Required_somewhere_else = []

		for key_out in sorted(blacklist_2.keys()) :
			if not key_in == key_out :
				Required_somewhere_else += forced_list_2[key_out]
				Required_somewhere_else += fixed_path_2[key_out]
			Required_somewhere_else += forced_list_1[key_out]
			Required_somewhere_else += fixed_path_1[key_out]

		for name in list(set(Required_somewhere_else)) :
			blacklist_2[key_in] += three_orientation_list(name)

		print >> sys.stderr , "### To reject for Path2@" + str(key_in) + " : " + ",".join(blacklist_2[key_in])

	# If markers a given, run consistency check, duplication control and reporting
	if options.marker_map :
		if not options.markers_hits :
			print >> sys.stdout , "[ERROR] Marker map present but no marker hits given"
			print >> sys.stderr , "[ERROR] Marker map present but no marker hits given"
			sys.exit(1)
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Testing markers uniqueness"

		if not options.disable_marker_ploidy_check :
			if not options.No2 :
				ploidy = 2
			else :
				ploidy = 1
		else :
			ploidy = 100

		sys.stdout, '[' + str(datetime.datetime.now()) + "] == Checking copy number"
		multi_copy_markers , unique_marker_hits_by_id = check_marker_copies( marker_hits_by_id , ploidy )
		#print >> sys.stderr , "multi_copy_markers len:" + str(len(multi_copy_markers))
		# multi_copy_markers: list of markers IDs in more copies than ploidy
		# unique_marker_hits_by_id:	dict of unique markers
		sys.stdout, '[' + str(datetime.datetime.now()) + "] == Checking duplications within sequence"

		chimera_list , markers_itradup , unique_distinct_marker_hits_by_id = check_in_sequence_duplications(marker_hits_by_seq , unique_marker_hits_by_id , multi_copy_markers , query_fasta_db , annotation_gff3 , query_agp_db , options.out , int(options.cores) , paths , options.skip_qc , ploidy )
		#json.dump(unique_distinct_marker_hits_by_id, sys.stderr , indent=4)
		print >> sys.stderr , "chimera_list len:" + str(len(chimera_list))
		print >> sys.stderr , "markers_itradup len:" + str(len(markers_itradup))
		#sys.exit(1)
		# chimera_list = list of chimeric sequences
		# markers_itradup = list of duplicated markers in chimeric sequences
		# unique_distinct_marker_hits_by_id = db of markers in #ploidy copies in #ploidy distinct sequences
		#  In dry mode, stop after marker QC
		if options.dry :
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Dry run completed"
			print >> sys.stderr, "# Dry run completed"
			print >> sys.stdout , "------------------------------"
			print >> sys.stdout,  "- Done"
			print >> sys.stdout , "------------------------------"
			print >> sys.stderr , "##############################"
			print >> sys.stderr , "# Done"
			print >> sys.stderr , "##############################"
			sys.exit(0)

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Cleaning markers from duplications"
		print >> sys.stderr , "# Cleaning markers from duplications"
		clean_marker_set_by_seq = clean_markers( unique_distinct_marker_hits_by_id , chimera_list , marker_map_by_seq , marker_map_by_id , options.out , query_fasta_db , {"show-coords" : paths["show-coords"] , "nucmer" : paths["nucmer"] } , int(options.cores) , options.filter_hits , options.extended_region , forced_list_1 , forced_list_2 , options.force_direction1 , options.force_direction2)
		#json.dump(clean_marker_set_by_seq, open( "clean_marker_set_by_seq.txt",'w'), indent=4)
		# Format:
		# 	clean_marker_set_by_seq[seq_id] = {
		#		clean_marker_set_by_seq[seq_id]["id"] : seq_id
		# 		clean_marker_set_by_seq[seq_id]["chr"] : chr_id
		# 		clean_marker_set_by_seq[seq_id]["markers"] : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ]
		# 		clean_marker_set_by_seq[seq_id]["markers"]["+"|"-"] with orientaion=="." : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ]
		# 		clean_marker_set_by_seq[seq_id]["range"] : [marker_pos_min , marker_pos_max] ,
		# 		clean_marker_set_by_seq[seq_id]["orientation"] : ["+" or "-" or "."] }

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Building tiling paths on markers"
		print >> sys.stderr , "# Building tiling paths on markers"
		# 1 - Split clean_marker_set_by_seq by chromosome
		clean_marker_set_by_chr = {}
		for seq_id in clean_marker_set_by_seq.keys() :
			chr_id = clean_marker_set_by_seq[seq_id]["chr"][0]
			if chr_id not in clean_marker_set_by_chr:
				clean_marker_set_by_chr[chr_id] = []
			clean_marker_set_by_chr[chr_id].append(clean_marker_set_by_seq[seq_id])
		#json.dump(clean_marker_set_by_seq, open("clean_marker_set_by_seq.txt" , 'w') , indent=4)

		unwanted_sequences_alternative_to_forced = {}
		for Target_sequence in forced_list_1.keys() :
			for seq_id in forced_list_1[Target_sequence] :
				if seq_id in alternative_sequences :
					alternative_to_seq_id = alternative_sequences[seq_id]
					for Target_sequence_2 in forced_list_1.keys() :
						if Target_sequence_2 not in unwanted_sequences_alternative_to_forced :
							unwanted_sequences_alternative_to_forced[Target_sequence_2] = []
						unwanted_sequences_alternative_to_forced[Target_sequence_2]+=alternative_to_seq_id

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Hap1 "
		print >> sys.stderr , "# Hap1 "

		for chr_id in sorted(clean_marker_set_by_chr.keys()) :
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Processing " + str(chr_id)
			print >> sys.stderr , "## Processing " + str(chr_id)
			# Add start and stop to markers list
			#print >> sys.stderr , "clean_marker_set_by_chr:"
			#print >> sys.stderr, clean_marker_set_by_chr[chr_id]
			stop_pos = max([ int(x[0]) for x in marker_map_by_seq[chr_id] ]) + 1
			oriented_marker_set = add_orientation_to_ID(clean_marker_set_by_chr[chr_id])
			#print >> sys.stderr , "oriented_marker_set:"
			#print >> sys.stderr, oriented_marker_set
			## QC of concordance of the forced list orientation with marker_set
			marker_set , validated_forced_list_1  , validated_forced_list_2 , validated_blacklist_1 , validated_blacklist_2 = validate_marker_set( oriented_marker_set , forced_list_1[chr_id] , forced_list_2[chr_id] , blacklist_1[chr_id] , blacklist_2[chr_id] , options.conflict_resolution )
			# validated_forced_list_X do not contain any sequence ID that has no marker
			# sequences previously missing orientation info are directed according to the force_list, if present
			#print >> sys.stderr , "marker_set:"
			#print >> sys.stderr, marker_set
			#json.dump(marker_set, open("validated_marker_set."+ chr_id +".txt" , 'w') , indent=4)

			if not options.use1 :
				preferred_db = {}
				# No given path for first tiling, search for first best tiling with blacklist and forced_list
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Searching best tiling path"
				print >> sys.stderr , "### Searching 1st best tiling path"
				# Check forced list with unwanted_pairs (both in -> rise error)
				unwanted_to_remove_list = []
				for constrained_seq_id in (unwanted_pairs.keys()) :
					if constrained_seq_id in validated_forced_list_1 :
						# sequence in an unwanted pair is in the forced set, check the mate
						mate_ids = unwanted_pairs[constrained_seq_id]
						if not is_list_overlapping_list( mate_ids , validated_forced_list_1 ):
							unwanted_to_remove_list += mate_ids
						else :
							# Both sequences requested to be in the same haplotype -> error!
							print >> sys.stdout , "[ERROR] Conflicting sequences constrains: Requested sequences pair present in exclusion information file"
							print >> sys.stderr , "[ERROR] Conflicting sequences constrains: Requested sequences pair present in exclusion information file"
							print >> sys.stderr , "[ERROR] Sequence pair present in both requested tiling path " + options.Require1 + " and excluded pairs from the same haplotype in file " + options.exclusion
							print >> sys.stderr , "[ERROR] Requested sequence: " + constrained_seq_id
							sys.exit(1)
				if chr_id in unwanted_sequences_alternative_to_forced :
					unwanted_to_remove_list += unwanted_sequences_alternative_to_forced[chr_id]
				best_found = False
				unwanted_to_remove_list = list(set(unwanted_to_remove_list))
				unwanted_to_reuse_list = []

				new_forced = validated_forced_list_1[:]
				while not best_found :
					#print >> sys.stderr , "unwanted_to_remove_list: "
					#print >> sys.stderr, unwanted_to_remove_list
					marker_graph_1 = remove_sequence_from_graph( unwanted_to_remove_list , chr_id , stop_pos , marker_set , new_forced , validated_blacklist_1 )
					#print >> sys.stderr, marker_graph_1
					best_1_nodes = nx.dag_longest_path(marker_graph_1, weight='marker_num')
					#print >> sys.stderr , best_1_nodes
					best_1_edges , best_1_edges_names = get_marker_subgraph_from_path( marker_graph_1 , best_1_nodes )
					#print >> sys.stderr , best_1_edges_names
					# Run unwanted pair search, eventually update graph and rerun search
					#print >> sys.stderr, "searching unwanted"
					unwanted_to_remove = search_unwanted( best_1_edges_names , unwanted_pairs , all_seq_length , validated_forced_list_1)

					if unwanted_to_remove == "" :
						all_preferred_used = True
						to_reintegrate = []
						preferred_to_remove = []
						# Check if the best path contains all preferred sequences >> otherwise remove them and try improving integrating the incompatible sequences

						for seq_id in preferred_db.keys() :
							if not seq_id in best_1_edges_names :
								print >> sys.stderr, "##### Found equal optimal path without " + seq_id + ", trying to improve by reintegrating the incompatible sequences"
								#print >> sys.stderr, "###### To reintegrate: " + str(preferred_db[seq_id])
								all_preferred_used = False
								unwanted_to_reuse_list += three_orientation_list(seq_id , True)
								preferred_to_remove.append(seq_id)
								to_reintegrate += preferred_db[seq_id]

						#print >> sys.stderr, "to_reintegrate: " + str(to_reintegrate)
						if all_preferred_used :
							# Cannot not be improved as best uses all preferred >> best found
							best_found = True
						else :
							# Pathway do not uses all preferred, it is equal to one using preferred but can be ameliorated if can make use of the sequences rejected by the prefereed
							# Unused preferred >> unwanted_to_remove_list (will not reduce the performance not using that sequence
							# Incompatible with unused preferred >> reintegrate to try to improve the pathway
							new_forced = best_1_edges_names
							to_reintegrate = list(set(to_reintegrate))
							for reintegrate_id in to_reintegrate :
								if (reintegrate_id in unwanted_to_remove_list) and (not reintegrate_id in unwanted_to_reuse_list) :
									unwanted_to_remove_list.remove(reintegrate_id)
									print >> sys.stderr, "###### Reintegrating: " + str(reintegrate_id)

							for seq_id in preferred_to_remove :
								del preferred_db[seq_id]
								unwanted_to_remove_list += three_orientation_list(seq_id , True)

							unwanted_to_remove_list = list(set(unwanted_to_remove_list))

					else :
						# unwanted_to_remove["keep"] = keep_seq_id|+
						# unwanted_to_remove["remove"] = [ reject_seq_id|+ , reject_seq_id|- , reject_seq_id|. ]
						print >> sys.stderr, "#### Found pair of sequences in the path that should not be placed in the same haplotype: " + str(unwanted_to_remove["keep"]) + " (to keep) | " + str(unwanted_to_remove["remove"])  + " (to remove)"

						# Check if the one we are going to remove was a preferred before, eventually reconsider the removal decision
						to_reintegrate = []
						if unwanted_to_remove["remove"][0] in preferred_db :
							to_reintegrate += preferred_db[unwanted_to_remove["remove"][0]]
							del preferred_db[unwanted_to_remove["remove"][0]]
						if unwanted_to_remove["remove"][1] in preferred_db :
							to_reintegrate += preferred_db[unwanted_to_remove["remove"][1]]
							del preferred_db[unwanted_to_remove["remove"][1]]
						if unwanted_to_remove["remove"][2] in preferred_db :
							to_reintegrate += preferred_db[unwanted_to_remove["remove"][2]]
							del preferred_db[unwanted_to_remove["remove"][2]]
							# The sequence to be removed was once preferred to other(s) sequence(s) in the path.
							# Reintegrate the deleted one(s), if they are compatible

						# Block - Start - Note - This should solve the issue of sequences IDs stuck in forced bucket even after incompatibility resolution
						updated_forced = []
						for element in new_forced :
							if not element in unwanted_to_remove["remove"] :
								updated_forced.append(element)
						# Block - End

						new_forced = updated_forced

						if to_reintegrate == [] :
							unwanted_to_remove_list += unwanted_to_remove["remove"]
						else :
							unwanted_to_remove_list = list(set(unwanted_to_remove_list))
							# Reintegrate the old sequences by removing them from the unwanted list
							for reintegrate_id in to_reintegrate :
								if (reintegrate_id in unwanted_to_remove_list) and (reintegrate_id not in unwanted_to_reuse_list) :
									unwanted_to_remove_list.remove(reintegrate_id)
									new_forced.append(reintegrate_id)
									print >> sys.stderr, "###### Reintegrating: " + reintegrate_id
							unwanted_to_remove_list += unwanted_to_remove["remove"]
						unwanted_to_remove_list = list(set(unwanted_to_remove_list))

						if unwanted_to_remove["keep"] not in preferred_db :
							preferred_db[unwanted_to_remove["keep"]] = []
						preferred_db[unwanted_to_remove["keep"]] += unwanted_to_remove["remove"]

				best_1_paths[chr_id] = best_1_edges
				used , best_1_paths_edges[chr_id] = make_list_from_marker_path( best_1_edges )
				print >> sys.stderr , "#### Used sequence (in order): " + ",".join([str(x) for x in best_1_paths_edges[chr_id]])
				all_used += used
				used_by_chr_hap1[chr_id] = used

		# Find all sequences paired to the sequences used in the best_1_paths
		paired_to_used = []
		if not options.rearrangements :
			print >> sys.stderr , "## Setting constrains based on sequences used for Hap 1"
			for Target_sequence in best_1_paths_edges.keys() :
				for seq_id in best_1_paths_edges[Target_sequence] :
					if seq_id in known_groups_by_seqid :
						paired_to_used += known_groups_by_seqid[seq_id]
						#print >> sys.stderr ,  "#### To discard: " + str(known_groups_by_seqid[seq_id])
			paired_to_used = list(set(paired_to_used))

		# Blacklist all sequences alternate to best_1_paths_edges[Target_sequence]
		unwanted_sequences_alternative_to_used = {}
		for Target_sequence in best_1_paths_edges.keys() :
			unwanted_sequences_alternative_to_used[Target_sequence] = []
			for seq_id in best_1_paths_edges[Target_sequence] :
				if seq_id in alternative_sequences :
					alternative_to_seq_id = alternative_sequences[seq_id]
					for Target_sequence_2 in best_1_paths_edges.keys() :
						if not Target_sequence_2 == Target_sequence :
							if Target_sequence_2 not in unwanted_sequences_alternative_to_used :
								unwanted_sequences_alternative_to_used[Target_sequence_2] = []
							unwanted_sequences_alternative_to_used[Target_sequence_2]+=alternative_to_seq_id

		unwanted_sequences_alternative_to_forced = {}
		for Target_sequence in forced_list_2.keys() :
			for seq_id in forced_list_2[Target_sequence] :
				if seq_id in alternative_sequences :
					alternative_to_seq_id = alternative_sequences[seq_id]
					for Target_sequence_2 in forced_list_2.keys() :
						if Target_sequence_2 not in unwanted_sequences_alternative_to_forced :
							unwanted_sequences_alternative_to_forced[Target_sequence_2] = []
						unwanted_sequences_alternative_to_forced[Target_sequence_2]+=alternative_to_seq_id

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Hap2 "
		print >> sys.stderr , "# Hap2 "
		for chr_id in sorted(clean_marker_set_by_chr.keys()) :
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Processing " + str(chr_id)
			print >> sys.stderr , "## Processing " + str(chr_id)
			stop_pos = max([ int(x[0]) for x in marker_map_by_seq[chr_id] ]) + 1
			oriented_marker_set = add_orientation_to_ID(clean_marker_set_by_chr[chr_id])
			marker_set , validated_forced_list_1  , validated_forced_list_2 , validated_blacklist_1 , validated_blacklist_2 = validate_marker_set( oriented_marker_set , forced_list_1[chr_id] , forced_list_2[chr_id] , blacklist_1[chr_id] , blacklist_2[chr_id] , options.conflict_resolution )

			if not options.use2 :
				if not options.No2 :
					preferred_db = {}
					# No given path for second tiling, search for second best tiling with blacklist and forced_list
					print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Searching best tiling path"
					print >> sys.stderr , "### Searching 2nd best tiling path"
					unwanted_to_remove_list = paired_to_used
					unwanted_to_remove_list += unwanted_sequences_alternative_to_used[chr_id]	
					if chr_id in unwanted_sequences_alternative_to_forced :
						unwanted_to_remove_list += unwanted_sequences_alternative_to_forced[chr_id]
					for seq_id in best_1_paths_edges[chr_id] :
						unwanted_to_remove_list += three_orientation_list(seq_id , True)
						if seq_id in known_groups_by_seqid :
							print >> sys.stderr ,  "#### To discard because of " + seq_id +" in hap1: " + str(known_groups_by_seqid[seq_id])
							unwanted_to_remove_list += known_groups_by_seqid[seq_id]

					# Check that no forced sequence is in unwanted_to_remove_list -> rise error and quit (Should never happen)
					for forced_seq_id in validated_forced_list_2 :
						if forced_seq_id in unwanted_to_remove_list :
							print >> sys.stdout , "[ERROR] Conflicting sequences constrains: Requested sequences should belong to Haplotype 1"
							print >> sys.stderr , "[ERROR] Requested sequence " + forced_seq_id + " is in conflict with Haplotype 1 reconstruction. From sequences known to be in the same haplotype reported in file " + options.known + ", the sequence should be along with others in Haplotype 1"
							sys.exit(1)
					# Check forced list with unwanted_pairs (both in -> rise error)
					for constrained_seq_id in (unwanted_pairs.keys()) :
						if constrained_seq_id in validated_forced_list_2 :
							# sequence in an unwanted pair is in the forced set, check the mate(s)
							mate_id = unwanted_pairs[constrained_seq_id]
							if not is_list_overlapping_list(mate_id,  validated_forced_list_2 ):
								unwanted_to_remove_list += mate_id
							else :
								# Both sequences requested to be in the same haplotype -> error!
								print >> sys.stdout , "[ERROR] Conflicting sequences constrains: Requested sequences pair present in exclusion information file"
								print >> sys.stderr , "[ERROR] Requested sequence " + constrained_seq_id + " present in both requested tiling path " + options.Require2 + " and excluded pairs from the same haplotype in file " + options.exclusion
								sys.exit(1)

					unwanted_to_remove_list = list(set(unwanted_to_remove_list))
					unwanted_to_reuse_list = []
					best_found = False
					new_forced = validated_forced_list_2[:]
					while not best_found :
						#print >> sys.stderr , "unwanted_to_remove_list: "
						#print >> sys.stderr, unwanted_to_remove_list
						#print >> sys.stderr , "validated_blacklist_2: "
						#print >> sys.stderr, validated_blacklist_2
						marker_graph_2 = remove_sequence_from_graph( unwanted_to_remove_list , chr_id , stop_pos, marker_set , new_forced , validated_blacklist_2 )
						best_2_nodes = nx.dag_longest_path(marker_graph_2, weight='marker_num')
						#print >> sys.stderr , best_2_nodes
						best_2_edges , best_2_edges_names = get_marker_subgraph_from_path( marker_graph_2 , best_2_nodes )
						# Run unwanted pair search, eventually update graph and rerun search
						#print >> sys.stderr , best_2_edges
						#print >> sys.stderr , best_2_edges_names
						unwanted_to_remove = search_unwanted( best_2_edges_names , unwanted_pairs , all_seq_length , validated_forced_list_2)

						if unwanted_to_remove == "" :
							all_preferred_used = True
							to_reintegrate = []
							preferred_to_remove = []
							# Check if the best path contains all preferred sequences >> otherwise remove them and try improving integrating the incompatible sequences

							for seq_id in preferred_db.keys() :
								if not seq_id in best_2_edges_names :
									print >> sys.stderr, "##### Found equal optimal path without " + seq_id + ", trying to improve by reintegrating the incompatible sequences"
									#print >> sys.stderr, "###### To reintegrate: " + str(preferred_db[seq_id])
									all_preferred_used = False
									unwanted_to_reuse_list += three_orientation_list(seq_id , True)
									preferred_to_remove.append(seq_id)
									to_reintegrate += preferred_db[seq_id]

							#print >> sys.stderr, "to_reintegrate: " + str(to_reintegrate)
							if all_preferred_used :
								# Cannot not be improved as best uses all preferred >> best found
								best_found = True
							else :
								# Pathway do not uses all preferred, it is equal to one using preferred.
								# Unused preferred >> unwanted_to_remove_list (will not reduce the performance not using that sequence
								# Incompatible with unused preferred >> reintegrate to try to improve the pathway
								new_forced = best_2_edges_names
								to_reintegrate = list(set(to_reintegrate))
								for reintegrate_id in to_reintegrate :
									if (reintegrate_id in unwanted_to_remove_list) and (not reintegrate_id in unwanted_to_reuse_list) :
										unwanted_to_remove_list.remove(reintegrate_id)
										print >> sys.stderr, "###### Reintegrating: " + str(reintegrate_id)

								for seq_id in preferred_to_remove :
									del preferred_db[seq_id]
									unwanted_to_remove_list += three_orientation_list(seq_id , True)

								unwanted_to_remove_list = list(set(unwanted_to_remove_list))

						else :
							# unwanted_to_remove["keep"] = keep_seq_id|+
							# unwanted_to_remove["remove"] = [ reject_seq_id|+ , reject_seq_id|- , reject_seq_id|. ]
							print >> sys.stderr, "#### Found pair of sequences in the path that should not be placed in the same haplotype: " + str(unwanted_to_remove["keep"]) + " (to keep) | " + str(unwanted_to_remove["remove"])  + " (to remove)"

							# Check if the one we are going to remove was a preferred before, eventually reconsider the removal decision
							to_reintegrate = []
							if unwanted_to_remove["remove"][0] in preferred_db :
								to_reintegrate += preferred_db[unwanted_to_remove["remove"][0]]
								del preferred_db[unwanted_to_remove["remove"][0]]
							if unwanted_to_remove["remove"][1] in preferred_db :
								to_reintegrate += preferred_db[unwanted_to_remove["remove"][1]]
								del preferred_db[unwanted_to_remove["remove"][1]]
							if unwanted_to_remove["remove"][2] in preferred_db :
								to_reintegrate += preferred_db[unwanted_to_remove["remove"][2]]
								del preferred_db[unwanted_to_remove["remove"][2]]
								# The sequence to be removed was once preferred to other(s) sequence(s) in the path.
								# Reintegrate the deleted one(s), if they are compatible

							if to_reintegrate == [] :
								unwanted_to_remove_list += unwanted_to_remove["remove"]
							else :
								unwanted_to_remove_list = list(set(unwanted_to_remove_list))
								# Reintegrate the old sequences by removing them from the unwanted list
								for reintegrate_id in to_reintegrate :
									if (reintegrate_id in unwanted_to_remove_list) and (reintegrate_id not in unwanted_to_reuse_list) :
										unwanted_to_remove_list.remove(reintegrate_id)
										print >> sys.stderr, "###### Reintegrating: " + reintegrate_id
								unwanted_to_remove_list += unwanted_to_remove["remove"]
							unwanted_to_remove_list = list(set(unwanted_to_remove_list))

							if unwanted_to_remove["keep"] not in preferred_db :
								preferred_db[unwanted_to_remove["keep"]] = []
							preferred_db[unwanted_to_remove["keep"]] += unwanted_to_remove["remove"]

					best_2_paths[chr_id] = best_2_edges
					used , best_2_paths_edges[chr_id] = make_list_from_marker_path( best_2_edges )
					print >> sys.stderr , "#### Used sequence (in order): " + ",".join([str(x) for x in best_2_paths_edges[chr_id]])
					all_used += used
					used_by_chr_hap2[chr_id] = used

					unused_by_chr[chr_id] = []
					for element in marker_set :
						seq_id = element["id"]
						if seq_id in all_used :
							continue
						else :
							unused_by_chr[chr_id].append(seq_id)

		# Return variables
		#	all_used
		#	best_1_paths
		#	best_1_paths_edges
		#	best_2_paths
		#	best_2_paths_edges

		# 	3 - If no collinearity analysis is required, prepare to export the paths warning for missing directional information
		if not options.reference :
			# Scan the edges to find the ones with missing order, convert them to "+" and issue a warning
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Filling missing orientation info for sequence in tiling paths"
			print >> sys.stderr , "## Filling missing orientation info for sequence in tiling paths"
			best_1_paths , best_1_paths_edges = fill_orientation( best_1_paths , best_1_paths_edges , options.out + ".missing_orientation.hap1.txt")
			best_2_paths , best_2_paths_edges = fill_orientation( best_2_paths , best_2_paths_edges , options.out + ".missing_orientation.hap2.txt")
	else :
		clean_marker_set_by_seq = {}


	# On reference analysis is dependent on markers tiling paths.
	if options.reference and not options.marker_map and not options.forced_as_path:
		#TODO: correct files output folder
		# Only guide genome is given, no markers. Run mapping and ex novo tiling paths identification

		# update blacklists using known and unwanted pairs
		forced_relationships_file = open( ".forced_relationships.txt" , 'w')
		for Target_sequence in sorted(forced_list_1.keys()) :
			for seq_id in forced_list_1[Target_sequence] :
				#print >> sys.stderr , "- " + Target_sequence + " - "  +seq_id
				if seq_id in unwanted_pairs :
					print >> forced_relationships_file , "-- " + Target_sequence + " - Hap1 - " + seq_id + " unwanted_pairs: " + str( unwanted_pairs[seq_id] )
					mate_ids = unwanted_pairs[seq_id]
					blacklist_1[Target_sequence] += mate_ids
					#print >> sys.stderr , "--- blacklist_1[" + Target_sequence + "]: "
					#print >> sys.stderr , blacklist_1[Target_sequence]

				if seq_id in known_groups_by_seqid :
					known_ids = known_groups_by_seqid[seq_id]
					print >> forced_relationships_file , "-- " + Target_sequence + " - Hap1 - " + seq_id + " known_ids: " + str(known_groups_by_seqid[seq_id] )
					for Target_sequence_2 in sorted(forced_list_1.keys()) :
						if not Target_sequence_2 == Target_sequence :
							blacklist_1[Target_sequence_2] += known_ids
							#print >> sys.stderr , "--- blacklist_1[" + Target_sequence_2 + "] :"
							#print >> sys.stderr , blacklist_1[Target_sequence_2]
					for Target_sequence_3 in sorted(forced_list_2.keys()) :
						blacklist_2[Target_sequence_3] += known_ids
						#print >> sys.stderr , "--- blacklist_2[" + Target_sequence_3 + "] :"
						#print >> sys.stderr , blacklist_2[Target_sequence_3]

				if seq_id in alternative_sequences :
					alternative_to_seq_id = alternative_sequences[seq_id]
					print >> forced_relationships_file , "-- " + Target_sequence + " - Hap1 - " + seq_id + " alternative_to_seq_id: " + str(alternative_to_seq_id)
					for Target_sequence_2 in forced_list_1.keys() :
						blacklist_1[Target_sequence_2] += alternative_to_seq_id
					for Target_sequence_3 in sorted(forced_list_2.keys()) :
						if not Target_sequence_3 == Target_sequence :
							blacklist_2[Target_sequence_3] += alternative_to_seq_id

		for Target_sequence in sorted(forced_list_2.keys()) :
			for seq_id in forced_list_2[Target_sequence] :
				#print >> sys.stderr , "- " + Target_sequence + " - seq_id:" + seq_id
				if seq_id in unwanted_pairs :
					print >> forced_relationships_file , "-- " + Target_sequence + " - Hap2 - " + seq_id + " unwanted_pairs: " + str( unwanted_pairs[seq_id] )
					mate_ids = unwanted_pairs[seq_id]
					blacklist_2[Target_sequence] += mate_ids
					#print >> sys.stderr , "--- blacklist_2[" + Target_sequence + "]: "
					#print >> sys.stderr , blacklist_2[Target_sequence]

				if seq_id in known_groups_by_seqid :
					known_ids = known_groups_by_seqid[seq_id]
					print >> forced_relationships_file , "-- " + Target_sequence + " - Hap2 - " + seq_id + " known_ids: " + str(known_groups_by_seqid[seq_id])
					for Target_sequence_2 in sorted(forced_list_2.keys()) :
						if not Target_sequence_2 == Target_sequence :
							blacklist_2[Target_sequence_2] += known_ids
							#print >> sys.stderr , "--- blacklist_2[" + Target_sequence_2 + "] :"
							#print >> sys.stderr , blacklist_2[Target_sequence_2]
					for Target_sequence_3 in sorted(forced_list_1.keys()) :
						blacklist_1[Target_sequence_3] += known_ids
						#print >> sys.stderr , "--- blacklist_1[" + Target_sequence_3 + "] :"
						#print >> sys.stderr , blacklist_1[Target_sequence_3]

				if seq_id in alternative_sequences :
					alternative_to_seq_id = alternative_sequences[seq_id]
					print >> forced_relationships_file , "-- " + Target_sequence + " - Hap2 - " + seq_id + " alternative_to_seq_id: " + str(alternative_to_seq_id)
					for Target_sequence_2 in forced_list_2.keys() :
						blacklist_2[Target_sequence_2] += alternative_to_seq_id
					for Target_sequence_3 in sorted(forced_list_2.keys()) :
						if not Target_sequence_3 == Target_sequence :
							blacklist_1[Target_sequence_3] += alternative_to_seq_id


		forced_relationships_file.close()

		for Target_sequence in blacklist_1.keys() :
			blacklist_1[Target_sequence] = list(set(blacklist_1[Target_sequence] ))
		for Target_sequence in blacklist_2.keys() :
			blacklist_2[Target_sequence] = list(set(blacklist_2[Target_sequence] ))

		# Checkpoint: Forced sequences present also in their own blacklist
		forcing_errors = []
		for Target_sequence in sorted(forced_list_1.keys()) :
			for seq_id in forced_list_1[Target_sequence] :
				if seq_id in blacklist_1[Target_sequence] :
					forcing_errors.append([Target_sequence , "hap1" , seq_id])
		for Target_sequence in sorted(forced_list_2.keys()) :
			for seq_id in forced_list_2[Target_sequence] :
				if seq_id in blacklist_2[Target_sequence] :
					forcing_errors.append([Target_sequence , "hap2" , seq_id])

		if not forcing_errors == [] :
			forcing_errors_file = open("forced_and_rejected_sequences.txt" , 'w')
			for element in sorted(forcing_errors) :
				print >> forcing_errors_file , "\t".join(element)
			forcing_errors_file.close()
			print >> sys.stdout, "[ERROR] Sequence(s) contemporary present in requested sequence list and the blacklist of a pseusomolecule"
			print >> sys.stderr, "[ERROR] Sequence(s) contemporary present in requested sequence list and the blacklist of a pseusomolecule. See forced_and_rejected_sequences.txt for more details"
			sys.exit(1)

		if options.map :
			doubleQuery = {}
			doubleQuery_len = {}
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Mapping query sequences on reference"

			query_for = "tmp." + os.path.basename(options.query) + ".for"

			qf_file = open(query_for , 'w' )

			for query_seq in query_fasta_db :
				print >> qf_file , ">" + query_seq + "|+"
				print >> qf_file , str(query_fasta_db[query_seq]).upper()
				doubleQuery[query_seq + "|+"]=str(query_fasta_db[query_seq]).upper()
				doubleQuery_len[query_seq + "|+"]=int(len(str(query_fasta_db[query_seq])))
			qf_file.close()

			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Mapping forward query sequences"
			map_for_file = open("tmp." + os.path.basename(options.hits) + ".for" , "w")

			minimap_path = paths["minimap2"]

			if minimap_path == "" :
				minimap2_search=subprocess.Popen( "which minimap2" , shell=True, stdout=subprocess.PIPE )
				command_line , error = minimap2_search.communicate()
				command_line = command_line.rstrip()
				if command_line == "" :
					print >> sys.stderr , '[ERROR] Minimap expected to be in $PATH, not found'
					exit(1)
			else :
				command_line = minimap_path + "/minimap2"


			if not os.path.exists(command_line) :
				print >> sys.stderr , "[ERROR] wrong or no path to minimap2"
				sys.exit(1)

			minimap2_command = command_line + " " + options.mapping + " "

			mappingCommand = minimap2_command + " --for-only -t " + str(options.cores) + " " + options.reference + " " + query_for + " | awk \'$5==\"+\"\'"
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Command line: " + mappingCommand + " > tmp." + options.hits + ".for"
			mapProcess = subprocess.Popen(mappingCommand, shell=True, stdout=map_for_file)
			output, error = mapProcess.communicate()

			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Reversing query sequences"
			query_rev = "tmp." + os.path.basename(options.query) + ".rev"

			qr_file = open(query_rev , 'w' )

			for query_seq in query_fasta_db :
				print >> qr_file , ">" + query_seq + "|-"
				print >> qr_file , str(Seq(query_fasta_db[query_seq]).reverse_complement()).upper()
				doubleQuery[query_seq+ "|-"]=str(Seq(query_fasta_db[query_seq]).reverse_complement()).upper()
				doubleQuery_len[query_seq+ "|-"]=int(len(str(query_fasta_db[query_seq])))
			qr_file.close()

			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Mapping reverse query sequences"
			map_rev_file = open("tmp." + os.path.basename(options.hits) + ".rev" , "w")
			mappingCommand = minimap2_command +  " --for-only -t " + str(options.cores) + " " + options.reference + " " + query_rev + " | awk \'$5==\"+\"\'"
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Command line: " + mappingCommand + " > tmp." + options.hits + ".rev"
			mapProcess = subprocess.Popen(mappingCommand, shell=True, stdout=map_rev_file)
			output, error = mapProcess.communicate()

			map_multi_file = open("tmp." + os.path.basename(options.hits) + ".multimapping" , "w")
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Merge forward and reverse alignments"
			mapProcess = subprocess.Popen("cat tmp." + os.path.basename(options.hits) + ".for tmp." + os.path.basename(options.hits) + ".rev", shell=True, stdout=map_multi_file)
			output, error = mapProcess.communicate()

			map_multi = "tmp." + os.path.basename(options.hits) + ".multimapping"
			mapProcess = subprocess.Popen("cp " + map_multi + " " + options.hits , shell=True)
			output, error = mapProcess.communicate()
		else :
			map_multi = "tmp." + os.path.basename(options.hits) + ".multimapping"
			mapProcess = subprocess.Popen("cp " + options.hits + " " + map_multi , shell=True)
			output, error = mapProcess.communicate()

			doubleQuery = {}
			doubleQuery_len = {}
			for query_seq in query_fasta_db :
				doubleQuery[query_seq+ "|+"]=str(query_fasta_db[query_seq]).upper()
				doubleQuery_len[query_seq+ "|+"]=int(len(doubleQuery[query_seq+ "|+"]))
				doubleQuery[query_seq+ "|-"]= str(Seq(query_fasta_db[query_seq]).reverse_complement()).upper()
				doubleQuery_len[query_seq+ "|-"]=int(len(doubleQuery[query_seq+ "|-"]))

		#### Merge and unique ####
		unique_hits = hit_mu( "tmp." + os.path.basename(options.hits) + ".multimapping" , "paf" , int(options.hitgap) , reference_len , doubleQuery_len)

		#print >> sys.stdout, unique_hits[unique_hits.keys()[0]]

		###### PAF format ######
		# https://github.com/lh3/miniasm/blob/master/PAF.md
		#
		# Tab separated format
		#
		#   0 - Query sequence name
		#   1 - Query sequence length
		#   2 - Query start (0-based)
		#   3 - Query end (0-based)
		#   4 -	Relative strand: "+" or "-"
		#   5 - Target sequence name
		#   6 - Target sequence length
		#   7 - Target start on original strand (0-based)
		#   8 - Target end on original strand (0-based)
		#   9 - Number of residue matches
		##  10 - Alignment block length
		##  11 - Mapping quality (0-255; 255 for missing)
		##
		##  [12:$] - Optional columns with flags
		##
		#########################

		###### Hit Format ######
		#
		# [ Qid , Tstart , Tstop , Qstart , Qstop , matches , hitLen ]
		#
		#	0 - Qid
		#	1 - Tstart
		#	2 - Tstop
		#	3 - Qstart
		#	4 - Qstop
		#	5 - matches
		#	6 - hitLen
		#
		########################

		### Read hits

		hits = {}
		for id in unique_hits.keys() :
			# unique_hits is in PAF format
			#   0 - Query sequence name
			#   1 - Query sequence length
			#   2 - Query start (0-based)
			#   3 - Query end (0-based)
			#   4 -	Relative strand: "+" or "-"
			#   5 - Target sequence name
			#   6 - Target sequence length
			#   7 - Target start on original strand (0-based)
			#   8 - Target end on original strand (0-based)
			#   9 - Number of residue matches
			#  10 - Alignment block length
			#  11 - Mapping quality (0-255; 255 for missing)
			# ['b40-14.HS_iter4_seq3248|+', 29417, 0, 29417, '+', 'chr13', 29075116, 23432640, 23461883, 27812, 29244, 13]
			Qid , Qlen , Qstart , Qstop , strand , Tid , Tlen , Tstart , Tstop , matches , hitLen , quality = unique_hits[id]

			if Tid not in hits :
				hits[Tid]=[]

			hits[Tid].append([ Qid , Tstart , Tstop , Qstart , Qstop , matches , hitLen ])

		# Generate graphs and extract paths
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Generating graphs and best paths"

		for Target_sequence in sorted(hits.keys()) :
			hits_1 = hits[Target_sequence]
			hits_1.append(["ChrStart" , 0 , 0 , 0 , 0 , 0 , 0 ])
			hits_1.append(["ChrStop" , reference_len[Target_sequence] , reference_len[Target_sequence] , 0 , 0 , 0 , 0 ])

			if not options.use1 :
				preferred_db = {}
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Graphing: " + Target_sequence
				print >> sys.stderr, Target_sequence
				map_total= len(hits_1)

				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Mapping query sequences: " + str(map_total)
				distance1 = int(options.distance1)
				unwanted_to_remove_list = blacklist_1[Target_sequence]

				best_found = False
				unwanted_to_remove_list = list(set(unwanted_to_remove_list))
				unwanted_to_reuse_list = []
				if Target_sequence in forced_list_1 :
					new_forced = forced_list_1[Target_sequence][:]
				else :
					new_forced = []
				while not best_found :
					if not new_forced == [] :
						hit_graph = make_forced_graph(hits_1, int(options.distance1) , new_forced , unwanted_to_remove_list  , options.Require1 + ".err.txt" )
					else :
						hit_graph = make_graph(hits_1, int(options.distance1) , unwanted_to_remove_list )
					print >> sys.stderr, hit_graph.edges.data()
					while not nx.has_path(hit_graph , source=0 , target=reference_len[Target_sequence] ) :
						distance1 += 100000
						print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Increased allowed gap to " + str(distance1)
						print >> sys.stderr, "#### Rebuilding increased allowed gap to " + str(distance1)
						if options.Require2 :
							hit_graph = make_forced_graph(hits_1, distance1 , forced_list_1[Target_sequence] , unwanted_to_remove_list  , options.Require1 + ".err.txt" )
						else :
							hit_graph = make_graph(hits_1, distance1 , unwanted_to_remove_list )
					# Select first
					print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Searching 1st path maximizing coverage: " + Target_sequence + " (0:" + str(reference_len[Target_sequence]) + ")"
					print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Graph order: " + str(hit_graph.order()) + "; Graph size: " + str(hit_graph.size())
					best_1_nodes = nx.dag_longest_path(hit_graph, weight='align')
					print >> sys.stderr , best_1_nodes
					best_1_edges , best_1_edges_names = get_subgraph_from_path( hit_graph , best_1_nodes )
					print >> sys.stderr , best_1_edges

					unwanted_to_remove = search_unwanted( best_1_edges_names , unwanted_pairs , all_seq_length , forced_list_1)

					if unwanted_to_remove == "" :
						all_preferred_used = True
						to_reintegrate = []
						preferred_to_remove = []
						# Check if the best path contains all preferred sequences >> otherwise remove them and try improving integrating the incompatible sequences

						for seq_id in preferred_db.keys() :
							if not seq_id in best_1_edges_names :
								print >> sys.stderr, "##### Found equal optimal path without " + seq_id + ", trying to improve by reintegrating the incompatible sequences"
								#print >> sys.stderr, "###### To reintegrate: " + str(preferred_db[seq_id])
								all_preferred_used = False
								unwanted_to_reuse_list += three_orientation_list(seq_id , True)
								preferred_to_remove.append(seq_id)
								to_reintegrate += preferred_db[seq_id]

						#print >> sys.stderr, "to_reintegrate: " + str(to_reintegrate)
						if all_preferred_used :
							# Cannot not be improved as best uses all preferred >> best found
							best_found = True
						else :
							# Pathway do not uses all preferred, it is equal to one using preferred.
							# Unused preferred >> unwanted_to_remove_list (will not reduce the performance not using that sequence
							# Incompatible with unused preferred >> reintegrate to try to improve the pathway
							new_forced = best_1_edges_names
							to_reintegrate = list(set(to_reintegrate))
							for reintegrate_id in to_reintegrate :
								if (reintegrate_id in unwanted_to_remove_list) and (not reintegrate_id in unwanted_to_reuse_list) :
									unwanted_to_remove_list.remove(reintegrate_id)
									print >> sys.stderr, "###### Reintegrating: " + str(reintegrate_id)

							for seq_id in preferred_to_remove :
								del preferred_db[seq_id]
								unwanted_to_remove_list += three_orientation_list(seq_id , True)

							unwanted_to_remove_list = list(set(unwanted_to_remove_list))

					else :
						# unwanted_to_remove["keep"] = keep_seq_id|+
						# unwanted_to_remove["remove"] = [ reject_seq_id|+ , reject_seq_id|- , reject_seq_id|. ]
						print >> sys.stderr, "#### Found pair of sequences in the path that should not be placed in the same haplotype: " + str(unwanted_to_remove["keep"]) + " (to keep) | " + str(unwanted_to_remove["remove"])  + " (to remove)"

						# Check if the one we are going to remove was a preferred before, eventually reconsider the removal decision
						to_reintegrate = []
						if unwanted_to_remove["remove"][0] in preferred_db :
							to_reintegrate += preferred_db[unwanted_to_remove["remove"][0]]
							del preferred_db[unwanted_to_remove["remove"][0]]
						if unwanted_to_remove["remove"][1] in preferred_db :
							to_reintegrate += preferred_db[unwanted_to_remove["remove"][1]]
							del preferred_db[unwanted_to_remove["remove"][1]]
						if unwanted_to_remove["remove"][2] in preferred_db :
							to_reintegrate += preferred_db[unwanted_to_remove["remove"][2]]
							del preferred_db[unwanted_to_remove["remove"][2]]
							# The sequence to be removed was once preferred to other(s) sequence(s) in the path.
							# Reintegrate the deleted one(s), if they are compatible

						if to_reintegrate == [] :
							unwanted_to_remove_list += unwanted_to_remove["remove"]
						else :
							unwanted_to_remove_list = list(set(unwanted_to_remove_list))
							# Reintegrate the old sequences by removing them from the unwanted list
							for reintegrate_id in to_reintegrate :
								if (reintegrate_id in unwanted_to_remove_list) and (reintegrate_id not in unwanted_to_reuse_list) :
									unwanted_to_remove_list.remove(reintegrate_id)
									print >> sys.stderr, "###### Reintegrating: " + reintegrate_id
							unwanted_to_remove_list += unwanted_to_remove["remove"]
						unwanted_to_remove_list = list(set(unwanted_to_remove_list))

						if unwanted_to_remove["keep"] not in preferred_db :
							preferred_db[unwanted_to_remove["keep"]] = []
						preferred_db[unwanted_to_remove["keep"]] += unwanted_to_remove["remove"]


				print >> sys.stderr ,"### Best tiling @ 1st iter (Haplotype 1): [" + "] -> [".join( ",".join(str(r) for r in x) for x in best_1_edges) + "]"
				best_1_paths[Target_sequence] = best_1_edges

			used , best_1_paths_edges[Target_sequence] = make_list_from_path( best_1_paths[Target_sequence] )
			all_used += used
			used_by_chr_hap1[Target_sequence] = used

			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Used for " + Target_sequence + ": " + ",".join(best_1_paths_edges[Target_sequence])
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Used so far: " + ",".join(all_used)

		# Find all sequences paired to the sequences used in the best_1_paths
		paired_to_used = []
		if not options.rearrangements :
			print >> sys.stderr , "## Setting constrains based on sequences used for Hap 1"
			for Target_sequence in best_1_paths_edges.keys() :
				for seq_id in best_1_paths_edges[Target_sequence] :
					if seq_id in known_groups_by_seqid :
						paired_to_used += known_groups_by_seqid[seq_id]
						#print >> sys.stderr ,  "#### To discard: " + str(known_groups_by_seqid[seq_id])
			paired_to_used = list(set(paired_to_used))

		# Blacklist all sequences alternate to best_1_paths_edges[Target_sequence]
		unwanted_sequences_alternative_to_used = {}
		for Target_sequence in best_1_paths_edges.keys() :
			unwanted_sequences_alternative_to_used[Target_sequence] = []
			for seq_id in best_1_paths_edges[Target_sequence] :
				if seq_id in alternative_sequences :
					alternative_to_seq_id = alternative_sequences[seq_id]
					for Target_sequence_2 in best_1_paths_edges.keys() :
						if not Target_sequence_2 == Target_sequence :
							if Target_sequence_2 not in unwanted_sequences_alternative_to_used :
								unwanted_sequences_alternative_to_used[Target_sequence_2] = []
							unwanted_sequences_alternative_to_used[Target_sequence_2]+=alternative_to_seq_id

		### Remove nodes form best path and search for second best tiling
		for Target_sequence in sorted(hits.keys()) :
			hits_2 = [x for x in hits[Target_sequence] if x[0] not in all_used ]
			unmap_2 = len(hits_2)

			if not options.use2 :
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Graphing: " + Target_sequence
				print >> sys.stderr, Target_sequence
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Updating graph removing best path for " + Target_sequence
				if not options.No2 :
					preferred_db = {}
					### Make new graph
					print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Recreating the graph"
					distance2 = int(options.distance2)
					unwanted_to_remove_list = blacklist_2[Target_sequence]
					unwanted_to_remove_list += paired_to_used
					unwanted_to_remove_list += unwanted_sequences_alternative_to_used[Target_sequence]
					unwanted_to_reuse_list = []
					# Blacklist the sequences that are known to be in hap1 along to the one sued
					for forced_seq_id in forced_list_2[Target_sequence] :
						if forced_seq_id in unwanted_to_remove_list :
							print >> sys.stdout , "[ERROR] Conflicting sequences constrains: Requested sequences also in black list"
							print >> sys.stderr , "[ERROR] Requested sequence " + forced_seq_id + " is also in black list"
							sys.exit(1)

					for forced_seq_id in forced_list_2[Target_sequence] :
						if forced_seq_id in unwanted_to_remove_list :
							print >> sys.stdout , "[ERROR] Conflicting sequences constrains: Requested sequences should belong to Haplotype 1"
							print >> sys.stderr , "[ERROR] Requested sequence " + forced_seq_id + " is in conflict with Haplotype 1 reconstruction. From sequences known to be in the same haplotype reported in file " + options.known + ", the sequence should be along with others in Haplotype 1"
							sys.exit(1)

					unwanted_to_remove_list = list(set(unwanted_to_remove_list))
					if Target_sequence in forced_list_2 :
						new_forced = forced_list_2[Target_sequence][:]
					else :
						new_forced = []
					best_found = False
					while not best_found :
						if not new_forced == [] :
							hit_graph_2 = make_forced_graph(hits_2, int(options.distance2) , new_forced , unwanted_to_remove_list , options.Require2 + ".err.txt" )
						else :
							hit_graph_2 = make_graph(hits_2, int(options.distance2) , unwanted_to_remove_list )

						# Check if graph is connected, otherwise increase allowed distance between nodes
						while not nx.has_path(hit_graph_2 , source=0 , target=reference_len[Target_sequence] ) :
							distance2 += 1000000
							print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Increased allowed gap to " + str(distance2)
							if options.Require2 :
								hit_graph_2 = make_forced_graph(hits_2, distance2 , forced_list_2[Target_sequence], unwanted_to_remove_list , options.Require2 + ".err.txt" )
							else :
								hit_graph_2 = make_graph(hits_2, distance2 , unwanted_to_remove_list )

						### Extract best 2nd
						print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Graph order: " + str(hit_graph_2.order()) + "; Graph size: " + str(hit_graph_2.size())
						print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Searching 2nd path maximizing coverage: " + Target_sequence + " (0:" + str(reference_len[Target_sequence]) + ")"
						best_2_nodes = nx.dag_longest_path(hit_graph_2, weight='align')
						#print >> sys.stderr , best_2_nodes
						best_2_edges , best_2_edges_names = get_subgraph_from_path( hit_graph_2 , best_2_nodes )
						#print >> sys.stderr , best_2_edges
						#print >> sys.stderr , best_2_edges_names

						unwanted_to_remove = search_unwanted( best_2_edges_names , unwanted_pairs , all_seq_length , forced_list_2)

						if unwanted_to_remove == "" :
							all_preferred_used = True
							to_reintegrate = []
							preferred_to_remove = []
							# Check if the best path contains all preferred sequences >> otherwise remove them and try improving integrating the incompatible sequences

							for seq_id in preferred_db.keys() :
								if not seq_id in best_2_edges_names :
									print >> sys.stderr, "##### Found equal optimal path without " + seq_id + ", trying to improve by reintegrating the incompatible sequences"
									#print >> sys.stderr, "###### To reintegrate: " + str(preferred_db[seq_id])
									all_preferred_used = False
									unwanted_to_reuse_list += three_orientation_list(seq_id , True)
									preferred_to_remove.append(seq_id)
									to_reintegrate += preferred_db[seq_id]

							#print >> sys.stderr, "to_reintegrate: " + str(to_reintegrate)
							if all_preferred_used :
								# Cannot not be improved as best uses all preferred >> best found
								best_found = True
							else :
								# Pathway do not uses all preferred, it is equal to one using preffered.
								# Unused preferred >> unwanted_to_remove_list (will not reduce the performance not using that sequence
								# Incompatible with unused preferred >> reintegrate to try to improve the pathway
								new_forced = best_2_edges_names
								to_reintegrate = list(set(to_reintegrate))
								for reintegrate_id in to_reintegrate :
									if (reintegrate_id in unwanted_to_remove_list) and (not reintegrate_id in unwanted_to_reuse_list) :
										unwanted_to_remove_list.remove(reintegrate_id)
										print >> sys.stderr, "###### Reintegrating: " + str(reintegrate_id)

								for seq_id in preferred_to_remove :
									del preferred_db[seq_id]
									unwanted_to_remove_list += three_orientation_list(seq_id , True)

								unwanted_to_remove_list = list(set(unwanted_to_remove_list))

						else :
							# unwanted_to_remove["keep"] = keep_seq_id|+
							# unwanted_to_remove["remove"] = [ reject_seq_id|+ , reject_seq_id|- , reject_seq_id|. ]
							print >> sys.stderr, "#### Found pair of sequences in the path that should not be placed in the same haplotype: " + str(unwanted_to_remove["keep"]) + " (to keep) | " + str(unwanted_to_remove["remove"])  + " (to remove)"

							# Check if the one we are going to remove was a preferred before, eventually reconsider the removal decision
							to_reintegrate = []
							if unwanted_to_remove["remove"][0] in preferred_db :
								to_reintegrate += preferred_db[unwanted_to_remove["remove"][0]]
								del preferred_db[unwanted_to_remove["remove"][0]]
							if unwanted_to_remove["remove"][1] in preferred_db :
								to_reintegrate += preferred_db[unwanted_to_remove["remove"][1]]
								del preferred_db[unwanted_to_remove["remove"][1]]
							if unwanted_to_remove["remove"][2] in preferred_db :
								to_reintegrate += preferred_db[unwanted_to_remove["remove"][2]]
								del preferred_db[unwanted_to_remove["remove"][2]]
								# The sequence to be removed was once preferred to other(s) sequence(s) in the path.
								# Reintegrate the deleted one(s), if they are compatible

							if to_reintegrate == [] :
								unwanted_to_remove_list += unwanted_to_remove["remove"]
							else :
								unwanted_to_remove_list = list(set(unwanted_to_remove_list))
								# Reintegrate the old sequences by removing them from the unwanted list
								for reintegrate_id in to_reintegrate :
									if (reintegrate_id in unwanted_to_remove_list) and (reintegrate_id not in unwanted_to_reuse_list) :
										unwanted_to_remove_list.remove(reintegrate_id)
										print >> sys.stderr, "###### Reintegrating: " + reintegrate_id
								unwanted_to_remove_list += unwanted_to_remove["remove"]
							unwanted_to_remove_list = list(set(unwanted_to_remove_list))

							if unwanted_to_remove["keep"] not in preferred_db :
								preferred_db[unwanted_to_remove["keep"]] = []
							preferred_db[unwanted_to_remove["keep"]] += unwanted_to_remove["remove"]

					print >> sys.stderr ,"### BEST Best tiling @ 2nd iter (Haplotype 2): [" + "] -> [".join( ",".join(str(r) for r in x) for x in best_2_edges) + "]"
					best_2_paths[Target_sequence] = best_2_edges

				else :
					best_2_edges_names = ["ChrStart" , "ChrStop"]


				used , best_2_paths_edges[Target_sequence] = make_list_from_path( best_2_paths[Target_sequence] )
				all_used += used
				used_by_chr_hap2[Target_sequence] = used

				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ===== Used for " + Target_sequence + ": " + ",".join(best_2_paths_edges[Target_sequence])
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ===== Used so far: " + ",".join(all_used)


				#unused_by_chr[Target_sequence] = []
				#for element in hits_2 :
				#	seq_id = element[0]
				#	if seq_id in all_used :
				#		continue
				#	else :
				#		unused_by_chr[Target_sequence].append(seq_id)

	elif (options.reference and options.marker_map) or ( options.reference and options.forced_as_path ) :
		# TODO: Check logging
		tmp_dir=options.out + ".tmp"
		mkdir(tmp_dir)
		orienation_db = {}
		fasta_files = {}
		fasta_len_dict = {}
		print >> sys.stderr, "## Writing sequences"
		for ref_id in reference.keys() :
			ref_fasta = reference[ref_id]
			ref_file_name = tmp_dir + "/" + ref_id + ".fasta"
			fasta_files[ref_id] = ref_file_name
			ref_file = open(ref_file_name , 'w')
			print >> ref_file, ">" + ref_id
			print >> ref_file, ref_fasta
			ref_file.close()
			fasta_len_dict[ref_id] = len(ref_fasta)

		map_used = []

		if options.marker_map :
			# Best tiling path with marker paths were already generated
			# Map unplaced everywhere and placed on their own chromosome
			# best_1_paths , best_1_paths_edges , best_2_paths , best_2_paths_edges are present but do contain sequences with no orientation (ID|.)

			intermediate_path_1_file = open(options.out + ".intermediate_hap1.list" , 'w')
			for chr_id in sorted(best_1_paths_edges.keys()) :
				print >> intermediate_path_1_file , chr_id + "\t" + ",".join(best_1_paths_edges[chr_id])
			intermediate_path_1_file.close()
			intermediate_path_2_file = open(options.out + ".intermediate_hap2.list" , 'w')
			for chr_id in sorted(best_2_paths_edges.keys()) :
				print >> intermediate_path_2_file , chr_id + "\t" + ",".join(best_2_paths_edges[chr_id])
			intermediate_path_2_file.close()

			## 1: find orientation for the sequences that do miss it
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Updating sequence orientation using guide genome "
			print >> sys.stderr, "# Updating sequence orientation using guide genome "
			undirected_sequences = get_no_orientation( best_1_paths_edges ,  best_2_paths_edges )
			#json.dump(undirected_sequences , sys.stdout)
			#print >> sys.stdout, ""
			# undirected_sequences[seq_id] = chr
			## Do map
			### Split sequences in different files
			for query_id in undirected_sequences.keys():
				query_fasta = query_fasta_db[query_id]
				query_file_name = tmp_dir + "/" + query_id + ".fasta"
				fasta_files[query_id] = query_file_name
				query_file = open(query_file_name , 'w')
				print >> query_file, ">" + query_id + "|+"
				print >> query_file, query_fasta.upper()
				print >> query_file , ">" + query_id + "|-"
				print >> query_file , str(Seq(query_fasta).reverse_complement()).upper()
				fasta_len_dict[query_id + "|+"] = len(query_fasta)
				fasta_len_dict[query_id + "|-"] = len(query_fasta)
				query_file.close()
			### Map each undirected sequence on tis own chromosome
			print >> sys.stderr, "## Mapping unoriented sequences on guide chromosomes"
			for query_id in undirected_sequences.keys():
				target_id = undirected_sequences[query_id]
				map_prefix = target_id + ".on." + query_id
				coords_file = tmp_dir + "/" + map_prefix + ".coords"
				coords_file = map_nucmer( fasta_files[target_id] , fasta_files[query_id] ,  int(options.cores) ,  coords_file , nucmer_path , showcoords_path , " --forward " , " -l -r -T -H ")
				map_results = read_nucmer_coords( coords_file )
				#### Find best alignment
				best_alignment = hits_best_tiling_path(map_results, fasta_len_dict)
				## Return dict of new orientations
				orienation_db[query_id] = best_orientation(best_alignment)
			# Update best_1_paths , best_1_paths_edges , best_2_paths , best_2_paths_edges
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Updating sequence orientation"
			print >> sys.stderr, "## Updating sequence orientation"
			best_1_paths , best_1_paths_edges = fill_orientation( best_1_paths , best_1_paths_edges , options.out + ".missing_orientation.hap1.txt" , orienation_db)
			best_2_paths , best_2_paths_edges = fill_orientation( best_2_paths , best_2_paths_edges , options.out + ".missing_orientation.hap2.txt" , orienation_db)

		else :
			if options.forced_as_path :
				for query_id in query_fasta_db.keys():
					fasta_len_dict[query_id + "|+"] = len(query_fasta_db[query_id])
					fasta_len_dict[query_id + "|-"] = len(query_fasta_db[query_id])

				broken_data = False
				best_1_paths_edges = {}
				best_2_paths_edges = {}
				all_used = []
				if options.Require1 :
					for line in open(options.Require1,"r") :
						try :
							Target_sequence , loaded_path = line.rstrip().split("\t")
							prompted_list = loaded_path.split(",")
							best_1_paths_edges[Target_sequence] = prompted_list
							for seq_id in prompted_list :
								if seq_id in all_used :
									print >> sys.stdout, '[ERROR] ' + seq_id + " required more than once"
									print >> sys.stderr, '[ERROR] ' + seq_id + " required more than once"
									broken_data = True
							all_used += prompted_list
						except :
							pass
				if options.Require2 :
					for line in open(options.Require2,"r") :
						try :
							Target_sequence , loaded_path = line.rstrip().split("\t")
							prompted_list = loaded_path.split(",")
							best_2_paths_edges[Target_sequence] = prompted_list
							for seq_id in prompted_list :
								if seq_id in all_used :
									print >> sys.stdout, '[ERROR] ' + seq_id + " required more than once"
									print >> sys.stderr, '[ERROR] ' + seq_id + " required more than once"
									broken_data = True
							all_used += prompted_list
						except :
							pass
				if broken_data :
					sys.exit(1)

			else :
				print >> sys.stdout , "[ERROR] Unknown option combination"
				print >> sys.stderr , "[ERROR] Unknown option combination"
				sys.exit(1)

		# Blacklist sequences known to be in the same pseudomolecule of the ones forced/in path
		for Target_sequence in best_1_paths_edges.keys() :
			for seq_id in best_1_paths_edges[Target_sequence] :
				if seq_id in known_groups_by_seqid:
					assigned_chr_ids = known_groups_by_seqid[seq_id]

					for Target_sequence_1 in blacklist_1.keys() :
						if not Target_sequence_1 == Target_sequence :
							blacklist_1[Target_sequence_1] += assigned_chr_ids
					for Target_sequence_2 in best_2_paths_edges.keys() :
						blacklist_2[Target_sequence_2] += assigned_chr_ids

		for Target_sequence in best_2_paths_edges.keys() :
			for seq_id in best_2_paths_edges[Target_sequence] :
				if seq_id in known_groups_by_seqid:
					assigned_chr_ids = known_groups_by_seqid[seq_id]

					for Target_sequence_1 in blacklist_1.keys() :
						blacklist_1[Target_sequence_1] += assigned_chr_ids
					for Target_sequence_2 in best_2_paths_edges.keys() :
						if not Target_sequence_2 == Target_sequence :
							blacklist_2[Target_sequence_2] += assigned_chr_ids

		#print >> sys.stderr , "- Updated Blacklists"
		for chr_id in blacklist_1.keys() :
			blacklist_1[chr_id] = list(set(blacklist_1[chr_id]))
			#print >> sys.stderr , "-- hap1 - " + chr_id + " - " + str(blacklist_1[chr_id])
		for chr_id in blacklist_2.keys() :
			blacklist_2[chr_id] = list(set(blacklist_2[chr_id]))
			#print >> sys.stderr , "-- hap2 - " + chr_id + " - " + str(blacklist_1[chr_id])

		# 2: generate intermediate sequences for tiling paths
		path_1_prefix = "intermediate1"
		path_2_prefix = "intermediate2"
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Generating intermediate pseudomolecules"
		print >> sys.stderr, "# Generating intermediate pseudomolecules"
		path1_db = generate_fasta_from_path(best_1_paths_edges , tmp_dir , path_1_prefix , query_fasta_db )
		path2_db = generate_fasta_from_path(best_2_paths_edges , tmp_dir , path_2_prefix , query_fasta_db )

		# pathX_db[chr] =
		# 	pathX_db[chr]["id"]
		# 	pathX_db[chr]["fasta_file"]
		# 	pathX_db[chr]["fasta_len"]
		# 	pathX_db[chr]["structure"] = [ ... , [concat_start , concat_stop , component_seq_id ] , ...  ]
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Mapping intermediate pseudomolecules to guide genome"
		print >> sys.stderr, "# Mapping intermediate pseudomolecules to guide genome"
		# 3: map intermediate tiling paths
		for chr_id in sorted(path1_db.keys()) :
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap1 - " + chr_id
			print >> sys.stderr, "## Hap1 - " + chr_id
			target_id = chr_id
			query_id = path1_db[chr_id]["id"]
			fasta_len_dict[query_id] = path1_db[chr_id]["fasta_len"]
			all_seq_length[query_id] = path1_db[chr_id]["fasta_len"]
			target_fasta = tmp_dir + "/" + target_id + ".fasta"
			query_fasta = path1_db[chr_id]["fasta_file"]
			coords_file = tmp_dir + "/" + query_id + ".on." + target_id + ".coords"
			if not options.reuse_intermediate :
				coords_file = map_nucmer( target_fasta , query_fasta ,  int(options.cores) ,  coords_file , nucmer_path , showcoords_path , " --forward " , " -l -r -T -H ")
			path1_db[chr_id]["coords_file"] = coords_file
			#### Find best alignment
			intermediate_hap1_hits_tiling_file = tmp_dir + "/" + query_id + ".on." + target_id + ".tiling_hits.tsv"
			if not options.reuse_intermediate :
				intermediate_hap1_hits_tiling = used_hits_best_tiling_path(coords_file, fasta_len_dict)[( target_id , query_id )]
			if not options.reuse_intermediate :
				write_table(intermediate_hap1_hits_tiling , intermediate_hap1_hits_tiling_file )
			path1_db[chr_id]["best_alignment"] = intermediate_hap1_hits_tiling_file
			## path1_db[chr] =
			# 	path1_db[chr]["id"]
			# 	path1_db[chr]["fasta_file"]
			# 	path1_db[chr]["fasta_len"]
			# 	path1_db[chr]["structure"] = [ ... , [concat_start , concat_stop , component_seq_id ] , ...  ]
			#	path1_db[chr_id]["coords_file"]
			#	path1_db[chr_id]["best_alignment"] = file >> [ ... , [ Tid , int(Tstart) , int(Tstop) , Qid, int(Qstart) , int(Qstop) , int(align_length) ,  int(match_length) ] , ... ]
		for chr_id in sorted(path2_db.keys()) :
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap2 - " + chr_id
			print >> sys.stderr, "## Hap2 - " + chr_id
			target_id = chr_id
			query_id = path2_db[chr_id]["id"]
			fasta_len_dict[query_id] = path2_db[chr_id]["fasta_len"]
			target_fasta = tmp_dir + "/" + target_id + ".fasta"
			query_fasta = path2_db[chr_id]["fasta_file"]
			all_seq_length[query_id] = path2_db[chr_id]["fasta_len"]
			coords_file = tmp_dir + "/" + query_id + ".on." + target_id + ".coords"
			if not options.reuse_intermediate :
				coords_file = map_nucmer( target_fasta , query_fasta ,  int(options.cores) ,  coords_file , nucmer_path , showcoords_path , " --forward " , " -l -r -T -H ")
			path2_db[chr_id]["coords_file"] = coords_file
			#### Find best alignment
			intermediate_hap2_hits_tiling_file = tmp_dir + "/" + query_id + ".on." + target_id + ".tiling_hits.tsv"
			if not options.reuse_intermediate :
				intermediate_hap2_hits_tiling = used_hits_best_tiling_path(coords_file, fasta_len_dict)[( target_id , query_id )]
			if not options.reuse_intermediate :
				write_table(intermediate_hap2_hits_tiling , intermediate_hap2_hits_tiling_file )
			path2_db[chr_id]["best_alignment"] = intermediate_hap2_hits_tiling_file
			## path2_db[chr] =
			# 	path2_db[chr]["id"]
			# 	path2_db[chr]["fasta_file"]
			# 	path2_db[chr]["fasta_len"]
			# 	path2_db[chr]["structure"] = [ ... , [concat_start , concat_stop , component_seq_id ] , ...  ]
			#	path2_db[chr_id]["coords_file"]
			#	path2_db[chr_id]["best_alignment"] = file >> [ ... , [ Tid , int(Tstart) , int(Tstop) , Qid, int(Qstart) , int(Qstop) , int(align_length) ,  int(match_length) ] , ... ]

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Translating alignments coordinates on input sequence"
		print >> sys.stderr, "# Translating alignments coordinates on input sequence"
		# 4: recover alignment information by sequence and list of sequences with/without maping info
		print >> sys.stderr, "## Hap1"
		path1_component_alignment , path1_mapped_sequences , path1_unmapped_sequences = get_component_alignment( path1_db )
		#json.dump( path1_component_alignment , open("path1_component_alignment.json" ,'w') , indent = 2)
		#json.dump( path1_mapped_sequences , open("path1_mapped_sequences.json" ,'w') , indent = 2)
		#json.dump( path1_unmapped_sequences , open("path1_unmapped_sequences.json" ,'w') , indent = 2)
		print >> sys.stderr, "## Hap2"
		path2_component_alignment , path2_mapped_sequences , path2_unmapped_sequences = get_component_alignment( path2_db )
		#json.dump( path2_component_alignment , open("path2_component_alignment.json" ,'w') , indent = 2)
		#json.dump( path2_mapped_sequences , open("path2_mapped_sequences.json" ,'w') , indent = 2)
		#json.dump( path2_unmapped_sequences , open("path2_unmapped_sequences.json" ,'w') , indent = 2)

		all_mapped = []
		for chr_id in path1_mapped_sequences.keys() :
			all_mapped += path1_mapped_sequences[chr_id]
		for chr_id in path2_mapped_sequences.keys() :
			all_mapped += path2_mapped_sequences[chr_id]

		# 5: map unplaced sequences
		## extract all unplaced
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Mapping unplaced sequences"
		print >> sys.stderr, "# Mapping unplaced sequences"
		unplaced_fasta = tmp_dir + "/intermediate_unplaced.fasta"
		unplaced_db = {}
		for seq_id in sorted(query_fasta_db.keys()) :
			if (seq_id + "|+" in  all_used) or (seq_id + "|-" in all_used) :
				continue
			else :
				unplaced_db[seq_id + "|+"] = query_fasta_db[seq_id].upper()
				unplaced_db[seq_id + "|-"] = str(Seq(query_fasta_db[seq_id]).reverse_complement()).upper()
			fasta_len_dict[query_id + "|+"] = len(query_fasta_db[seq_id])
			fasta_len_dict[query_id + "|-"] = len(query_fasta_db[seq_id])
		write_fasta_from_db( unplaced_db , unplaced_fasta )

		## Map unplaced
		coords_file = tmp_dir + "/unplaced.on.guide.coords"
		if not options.reuse_intermediate :
			coords_file = map_nucmer( options.reference , unplaced_fasta ,  int(options.cores) ,  coords_file , nucmer_path , showcoords_path , " --forward " , " -l -r -T -H ")
			## extract best alignment region
		unique_hits = hit_mu(coords_file , "coords" , int(options.hitgap) , all_seq_length , all_seq_length)

		unplaced_unique_hits = {}
		for id in unique_hits.keys() :
			# unique_hits is in PAF format
			#   0 - Query sequence name
			#   1 - Query sequence length
			#   2 - Query start (0-based)
			#   3 - Query end (0-based)
			#   4 -	Relative strand: "+" or "-"
			#   5 - Target sequence name
			#   6 - Target sequence length
			#   7 - Target start on original strand (0-based)
			#   8 - Target end on original strand (0-based)
			#   9 - Number of residue matches
			#  10 - Alignment block length
			#  11 - Mapping quality (0-255; 255 for missing)
			# ['b40-14.HS_iter4_seq3248|+', 29417, 0, 29417, '+', 'chr13', 29075116, 23432640, 23461883, 27812, 29244, 13]
			Qid , Qlen , Qstart , Qstop , strand , Tid , Tlen , Tstart , Tstop , matches , hitLen , quality = unique_hits[id]
			if Tid not in unplaced_unique_hits :
				unplaced_unique_hits[Tid]=[]
			unplaced_unique_hits[Tid].append([ Qid , Tstart , Tstop , Qstart , Qstop , matches , hitLen ])

		# unplaced_unique_hits[Tid] =
		# 	[ 	... ,
		# 		[ Qid , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ] ,
	 	#		... ]

		# 6: Generate network of first haplotype
		#	Respect forced_list_1[chr_id] , blacklist_1[chr_id] , forced_list_2[chr_id] , blacklist_2[chr_id]
		## Add sequences placed with markers to forced_list_1

		all_mapped += get_mapped_ids(unplaced_unique_hits)

		# Clean up forced lists of unmapped sequences
		clean_forced_list_1 = remove_missing_from_list(forced_list_1, all_mapped , "required sequences for haplotype 1")
		clean_forced_list_2 = remove_missing_from_list(forced_list_2, all_mapped , "required sequences for haplotype 2")

		#json.dump(clean_forced_list_1 , open("clean_forced_list_1.json",'w') , indent= 2 )
		#json.dump(clean_forced_list_2 , open("clean_forced_list_2.json",'w') , indent= 2 )

		merged_forced_list_1 = merge_sorted_lists( path1_mapped_sequences , clean_forced_list_1 )
		merged_forced_list_2 = merge_sorted_lists( path2_mapped_sequences , clean_forced_list_2 )

		#json.dump(merged_forced_list_1 , open("merged_forced_list_1.json",'w') , indent= 2 )
		#json.dump(merged_forced_list_2 , open("merged_forced_list_2.json",'w') , indent= 2 )

		new_forced_list_1 = remove_missing_from_list(merged_forced_list_1, all_mapped , "required sequences for haplotype 1")
		new_forced_list_2 = remove_missing_from_list(merged_forced_list_2, all_mapped , "required sequences for haplotype 2")

		#json.dump(new_forced_list_1 , open("new_forced_list_1.json",'w') , indent= 2 )
		#json.dump(new_forced_list_2 , open("new_forced_list_2.json",'w') , indent= 2 )

		unwanted_sequences_alternative_to_forced = {}
		for Target_sequence in new_forced_list_1.keys() :
			for seq_id in new_forced_list_1[Target_sequence] :
				if seq_id in alternative_sequences :
					alternative_to_seq_id = alternative_sequences[seq_id]
					for Target_sequence_2 in new_forced_list_1.keys() :
						if Target_sequence_2 not in unwanted_sequences_alternative_to_forced :
							unwanted_sequences_alternative_to_forced[Target_sequence_2] = []
						unwanted_sequences_alternative_to_forced[Target_sequence_2]+=alternative_to_seq_id

		for Target_sequence in new_forced_list_2.keys() :
			for seq_id in new_forced_list_2[Target_sequence] :
				if seq_id in alternative_sequences :
					alternative_to_seq_id = alternative_sequences[seq_id]
					for Target_sequence_1 in new_forced_list_2.keys() :
						if Target_sequence_1 not in unwanted_sequences_alternative_to_forced :
							unwanted_sequences_alternative_to_forced[Target_sequence_1] = []
						unwanted_sequences_alternative_to_forced[Target_sequence_1]+=alternative_to_seq_id

		## Add unwanted mates of forced_list_1 to blacklist_1
		for chr_id in sorted(new_forced_list_1.keys()) :
			#print >> sys.stderr, "- chr_id: " + chr_id
			for forced_seq_id in new_forced_list_1[chr_id] :
				if forced_seq_id in unwanted_pairs :
					blacklist_1[chr_id] += unwanted_pairs[forced_seq_id]
					#print >> sys.stderr, "-- forced seq_id: " + forced_seq_id + " - unwanted ids: " + str(unwanted_pairs[forced_seq_id])
				for chr_id_2 in sorted(new_forced_list_1.keys()) :
					if chr_id_2 == chr_id :
						continue
					else :
						#print >> sys.stderr, "--- hap1 - " + chr_id_2
						blacklist_1[chr_id_2] += three_orientation_list( forced_seq_id , True )
				for chr_id_3 in sorted(new_forced_list_2.keys()) :
					#print >> sys.stderr, "--- hap2 - " + chr_id_3
					blacklist_2[chr_id_3] += three_orientation_list( forced_seq_id , True )

		for chr_id in sorted(new_forced_list_2.keys()) :
			for forced_seq_id in new_forced_list_2[chr_id] :
				if forced_seq_id in unwanted_pairs :
					blacklist_2[chr_id] += unwanted_pairs[forced_seq_id]
					#print >> sys.stderr, "-- forced seq_id: " + forced_seq_id + " - unwanted ids: " + str(unwanted_pairs[forced_seq_id])
				for chr_id_2 in sorted(new_forced_list_2.keys()) :
					if chr_id_2 == chr_id :
						continue
					else :
						#print >> sys.stderr, "--- hap2 - " + chr_id_2
						blacklist_2[chr_id_2] += three_orientation_list( forced_seq_id , True )
				for chr_id_3 in sorted(new_forced_list_1.keys()) :
					#print >> sys.stderr, "--- hap1 - " + chr_id_3
					blacklist_1[chr_id_3] += three_orientation_list( forced_seq_id , True )

		#print >> sys.stderr , "- Updated Blacklists"
		for chr_id in blacklist_1.keys() :
			blacklist_1[chr_id] = list(set(blacklist_1[chr_id]))
			#print >> sys.stderr , "-- hap1 - " + chr_id + " - " + str(blacklist_1[chr_id])
		for chr_id in blacklist_2.keys() :
			blacklist_2[chr_id] = list(set(blacklist_2[chr_id]))
			#print >> sys.stderr , "-- hap2 - " + chr_id + " - " + str(blacklist_1[chr_id])

		## Split by chromosome, filter unusable unplaced sequences and merge to marker tilings alignments
		hits_1_by_chromosome = {}
		# unplaced_unique_hits[Tid] =
		# 	[ 	... ,
		# 		[ Qid , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ] ,
	 	#		... ]

		# path1_component_alignment[Chr] =
		# 	[ 	... ,
		# 		[ sed_id|orintation , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ] ,
		#		... ]

		for chr_id in unplaced_unique_hits.keys() :
			hits_1_by_chromosome[chr_id] = unplaced_unique_hits[chr_id]
		for chr_id in path1_component_alignment :
			if chr_id not in hits_1_by_chromosome :
				hits_1_by_chromosome[chr_id] = []
			hits_1_by_chromosome[chr_id] += path1_component_alignment[chr_id]


		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Generating tiling paths"
		print >> sys.stderr, "# Generating tiling paths"
		## make 1st tiling path
		if not options.use1 :
			preferred_db = {}
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap1 "
			print >> sys.stderr, "## Hap1 "
			new_best_1_paths = {}
			new_best_1_paths_edges = {}
			for Target_sequence in sorted(hits_1_by_chromosome.keys()):
				hits_1 = hits_1_by_chromosome[Target_sequence]
				hits_1.append(["ChrStart" , 0 , 0 , 0 , 0 , 0 , 0 ])
				hits_1.append(["ChrStop" , reference_len[Target_sequence] , reference_len[Target_sequence] , 0 , 0 , 0 , 0 ])
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Graphing: " + Target_sequence
				print >> sys.stderr, Target_sequence

				distance1 = int(options.distance1)
				unwanted_to_remove_list = blacklist_1[Target_sequence]
				if Target_sequence in unwanted_sequences_alternative_to_forced :
					unwanted_to_remove_list += unwanted_sequences_alternative_to_forced[Target_sequence]
				unwanted_to_remove_list = list(set(unwanted_to_remove_list))
				unwanted_to_reuse_list = []
				best_found = False
				new_forced = new_forced_list_1[Target_sequence][:]
				while not best_found :
					hit_graph = make_forced_graph(hits_1, distance1 , new_forced , unwanted_to_remove_list  , "markers_vs_alignment.hap1.incompatibility.txt" )
					#print >> sys.stderr, hit_graph.edges.data()
					while not nx.has_path(hit_graph , source=0 , target=reference_len[Target_sequence] ) :
						distance1 += 100000
						print >> sys.stderr, "#### Rebuilding increased allowed gap to " + str(distance1)
						hit_graph = make_forced_graph(hits_1, distance1 , new_forced_list_1[Target_sequence] , unwanted_to_remove_list  , "markers_vs_alignment.hap1.incompatibility.txt" )
					# Select first
					best_1_nodes = nx.dag_longest_path(hit_graph, weight='align')
					best_1_edges , best_1_edges_names = get_subgraph_from_path( hit_graph , best_1_nodes )
					unwanted_to_remove = search_unwanted( best_1_edges_names , unwanted_pairs , all_seq_length , new_forced_list_1[Target_sequence] + forced_list_1[Target_sequence])

					if unwanted_to_remove == "" :
						all_preferred_used = True
						to_reintegrate = []
						preferred_to_remove = []
						# Check if the best path contains all preferred sequences >> otherwise remove them and try improving integrating the incompatible sequences

						for seq_id in preferred_db.keys() :
							if not seq_id in best_1_edges_names :
								print >> sys.stderr, "##### Found equal optimal path without " + seq_id + ", trying to improve by reintegrating the incompatible sequences"
								#print >> sys.stderr, "###### To reintegrate: " + str(preferred_db[seq_id])
								all_preferred_used = False
								unwanted_to_reuse_list += three_orientation_list(seq_id , True)
								preferred_to_remove.append(seq_id)
								to_reintegrate += preferred_db[seq_id]

						#print >> sys.stderr, "to_reintegrate: " + str(to_reintegrate)
						if all_preferred_used :
							# Cannot not be improved as best uses all preferred >> best found
							best_found = True
						else :
							# Pathway do not uses all preferred, it is equal to one using preffered.
							# Unused preferred >> unwanted_to_remove_list (will not reduce the performance not using that sequence
							# Incompatible with unused preferred >> reintegrate to try to improve the pathway
							new_forced = best_1_edges_names
							to_reintegrate = list(set(to_reintegrate))
							for reintegrate_id in to_reintegrate :
								if (reintegrate_id in unwanted_to_remove_list) and (not reintegrate_id in unwanted_to_reuse_list) :
									unwanted_to_remove_list.remove(reintegrate_id)
									print >> sys.stderr, "###### Reintegrating: " + str(reintegrate_id)

							for seq_id in preferred_to_remove :
								del preferred_db[seq_id]
								unwanted_to_remove_list += three_orientation_list(seq_id , True)

							unwanted_to_remove_list = list(set(unwanted_to_remove_list))

					else :
						# unwanted_to_remove["keep"] = keep_seq_id|+
						# unwanted_to_remove["remove"] = [ reject_seq_id|+ , reject_seq_id|- , reject_seq_id|. ]
						print >> sys.stderr, "#### Found pair of sequences in the path that should not be placed in the same haplotype: " + str(unwanted_to_remove["keep"]) + " (to keep) | " + str(unwanted_to_remove["remove"])  + " (to remove)"

						# Check if the one we are going to remove was a preferred before, eventually reconsider the removal decision
						to_reintegrate = []
						if unwanted_to_remove["remove"][0] in preferred_db :
							to_reintegrate += preferred_db[unwanted_to_remove["remove"][0]]
							del preferred_db[unwanted_to_remove["remove"][0]]
						if unwanted_to_remove["remove"][1] in preferred_db :
							to_reintegrate += preferred_db[unwanted_to_remove["remove"][1]]
							del preferred_db[unwanted_to_remove["remove"][1]]
						if unwanted_to_remove["remove"][2] in preferred_db :
							to_reintegrate += preferred_db[unwanted_to_remove["remove"][2]]
							del preferred_db[unwanted_to_remove["remove"][2]]
							# The sequence to be removed was once preferred to other(s) sequence(s) in the path.
							# Reintegrate the deleted one(s), if they are compatible

						if to_reintegrate == [] :
							unwanted_to_remove_list += unwanted_to_remove["remove"]
						else :
							unwanted_to_remove_list = list(set(unwanted_to_remove_list))
							# Reintegrate the old sequences by removing them from the unwanted list
							for reintegrate_id in to_reintegrate :
								if (reintegrate_id in unwanted_to_remove_list) and (reintegrate_id not in unwanted_to_reuse_list) :
									unwanted_to_remove_list.remove(reintegrate_id)
									print >> sys.stderr, "###### Reintegrating: " + reintegrate_id
							unwanted_to_remove_list += unwanted_to_remove["remove"]
						unwanted_to_remove_list = list(set(unwanted_to_remove_list))

						if unwanted_to_remove["keep"] not in preferred_db :
							preferred_db[unwanted_to_remove["keep"]] = []
						preferred_db[unwanted_to_remove["keep"]] += unwanted_to_remove["remove"]

				print >> sys.stderr ,"### Best tiling @ 1st iter (Haplotype 1): [" + "] -> [".join( ",".join(str(r) for r in x) for x in best_1_edges) + "]"
				new_best_1_paths[Target_sequence] = best_1_edges

				used , new_best_1_paths_edges[Target_sequence] = make_list_from_path( new_best_1_paths[Target_sequence] )
				map_used += used
				used_by_chr_hap1[Target_sequence] = used

				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Used for " + Target_sequence + ": " + ",".join(new_best_1_paths_edges[Target_sequence])
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Used so far: " + ",".join(map_used)

		# 7: Generate network of second haplotype
		# 	Respect forced_list_1[chr_id] , blacklist_1[chr_id] , forced_list_2[chr_id] , blacklist_2[chr_id]
		else :
			new_best_1_paths_edges = best_1_paths_edges
			new_best_1_paths = best_1_paths

		# Find all sequences apired to the sequences used in the best_1_paths
		paired_to_used = []
		if not options.rearrangements :
			print >> sys.stderr , "## Setting constrains based on sequences used for Hap 1"
			for Target_sequence in new_best_1_paths_edges.keys() :
				for seq_id in new_best_1_paths_edges[Target_sequence] :
					if seq_id in known_groups_by_seqid :
						paired_to_used += known_groups_by_seqid[seq_id]
						#print >> sys.stderr ,  "#### To discard: " + str(known_groups_by_seqid[seq_id])
			paired_to_used = list(set(paired_to_used))

		# Blacklist all sequences alternate to best_1_paths_edges[Target_sequence]
		unwanted_sequences_alternative_to_used = {}
		for Target_sequence in best_1_paths_edges.keys() :
			unwanted_sequences_alternative_to_used[Target_sequence] = []
			for seq_id in best_1_paths_edges[Target_sequence] :
				if seq_id in alternative_sequences :
					alternative_to_seq_id = alternative_sequences[seq_id]
					for Target_sequence_2 in best_1_paths_edges.keys() :
						if not Target_sequence_2 == Target_sequence :
							if Target_sequence_2 not in unwanted_sequences_alternative_to_used :
								unwanted_sequences_alternative_to_used[Target_sequence_2] = []
							unwanted_sequences_alternative_to_used[Target_sequence_2]+=alternative_to_seq_id

		unwanted_sequences_alternative_to_forced = {}
		for Target_sequence in new_forced_list_2.keys() :
			for seq_id in new_forced_list_2[Target_sequence] :
				if seq_id in alternative_sequences :
					alternative_to_seq_id = alternative_sequences[seq_id]
					for Target_sequence_2 in new_forced_list_2.keys() :
						if Target_sequence_2 not in unwanted_sequences_alternative_to_forced :
							unwanted_sequences_alternative_to_forced[Target_sequence_2] = []
						unwanted_sequences_alternative_to_forced[Target_sequence_2]+=alternative_to_seq_id

		# Add path_1 sequences from new_best_1_paths_edges to blacklists, check if any forced sequence was used
		# new_best_1_paths_edges[chr_id] = [ ... , seq_id|orientation , ... ]
		for chr_id in new_best_1_paths_edges :
			for used_seq_id in new_best_1_paths_edges[chr_id]:
				for used_name in three_orientation_list(used_seq_id , True) :
					if used_name in forced_list_2_id_oriented :
						print >> sys.stdout , "[ERROR] Haplotype 1 reconstruction uses a sequence required dor Haplotype 2"
						print >> sys.stderr , "[ERROR] Haplotype 1 reconstruction uses a sequence required dor Haplotype 2 (" + used_name + ")"
						sys.exit(1)
				for chr_id_2 in blacklist_2.keys() :
					blacklist_2[chr_id_2] += three_orientation_list(used_seq_id , True)
		for chr_id in blacklist_2.keys() :
			blacklist_2[chr_id] = list(set(blacklist_2[chr_id]))


		# Split by chromosome,
		#	filter unusable unplaced sequences and sequences already used for path_1, merge to marker tilings alignments
		hits_2_by_chromosome = {}
		# unplaced_unique_hits[Tid] =
		# 	[ 	... ,
		# 		[ Qid , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ] ,
	 	#		... ]
		# path2_component_alignment[Chr] =
		# 	[ 	... ,
		# 		[ sed_id|orientation , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ] ,
		#		... ]

		for chr_id in unplaced_unique_hits.keys() :
			hits_2_by_chromosome[chr_id] = []
			for hit in unplaced_unique_hits[chr_id] :
				seq_id = hit[0]
				if (seq_id in map_used) or (seq_id in blacklist_2[chr_id]) :
					continue
				else :
					hits_2_by_chromosome[chr_id].append(hit)
		for chr_id in path2_component_alignment :
			if chr_id not in hits_2_by_chromosome :
				hits_2_by_chromosome[chr_id] = []
			hits_2_by_chromosome[chr_id] += path2_component_alignment[chr_id]

		## make 2nd tiling path
		if not options.use2 :
			if not options.No2 :
				preferred_db = {}
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap2 "
				print >> sys.stderr, "## Hap2 "
				new_best_2_paths = {}
				new_best_2_paths_edges = {}
				for Target_sequence in sorted(hits_2_by_chromosome.keys()):
					hits_2 = hits_2_by_chromosome[Target_sequence]
					hits_2.append(["ChrStart" , 0 , 0 , 0 , 0 , 0 , 0 ])
					hits_2.append(["ChrStop" , reference_len[Target_sequence] , reference_len[Target_sequence] , 0 , 0 , 0 , 0 ])
					print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Graphing: " + Target_sequence
					print >> sys.stderr, Target_sequence

					distance2 = int(options.distance2)
					unwanted_to_remove_list = blacklist_2[Target_sequence]
					unwanted_to_remove_list += paired_to_used
					unwanted_to_remove_list += unwanted_sequences_alternative_to_used[Target_sequence]
					if Target_sequence in unwanted_sequences_alternative_to_forced :
						unwanted_to_remove_list += unwanted_sequences_alternative_to_forced[Target_sequence]
					# Blacklist the sequences that are known to be in hap1 along to the one sued
					for seq_id in new_best_1_paths_edges[Target_sequence] :
						if seq_id in known_groups_by_seqid :
							unwanted_to_remove_list += known_groups_by_seqid[seq_id]
							print >> sys.stderr ,  "#### To discard because of " + seq_id + " in hap1: " + str(known_groups_by_seqid[seq_id])
					for forced_seq_id in new_forced_list_2[Target_sequence] :
						if forced_seq_id in unwanted_to_remove_list :
							print >> sys.stdout , "[ERROR] Conflicting sequences constrains: Requested sequences should belong to Haplotype 1"
							print >> sys.stderr , "[ERROR] Requested sequence " + forced_seq_id + " is in conflict with Haplotype 1 reconstruction as known to belong to the same haplotype."
							sys.exit(1)

					unwanted_to_remove_list = list(set(unwanted_to_remove_list))
					unwanted_to_reuse_list = []
					new_forced = new_forced_list_2[Target_sequence][:]
					best_found = False
					while not best_found :
						hit_graph_2 = make_forced_graph(hits_2, distance2 , new_forced , unwanted_to_remove_list  , "markers_vs_alignment.hap2.incompatibility.txt" )
						json.dump( sorted(hits_2) , open("last_hits.json",'w') )
						print >> sys.stderr, sorted(hit_graph_2.edges.data())
						json.dump( sorted(hit_graph_2.edges.data()) , open("last_graph.json",'w') )
						while not nx.has_path(hit_graph_2 , source=0 , target=reference_len[Target_sequence] ) :
							distance2 += 100000
							print >> sys.stderr, "#### Rebuilding increased allowed gap to " + str(distance2)
							hit_graph_2 = make_forced_graph(hits_2, distance2 , new_forced_list_2[Target_sequence] , unwanted_to_remove_list  , "markers_vs_alignment.hap2.incompatibility.txt" )
						# Select first
						best_2_nodes = nx.dag_longest_path(hit_graph_2, weight='align')
						best_2_edges , best_2_edges_names = get_subgraph_from_path( hit_graph_2 , best_2_nodes )
						unwanted_to_remove = search_unwanted( best_2_edges_names , unwanted_pairs , all_seq_length , new_forced_list_2[Target_sequence] + forced_list_2[Target_sequence])

						if unwanted_to_remove == "" :
							all_preferred_used = True
							to_reintegrate = []
							preferred_to_remove = []
							# Check if the best path contains all preferred sequences >> otherwise remove them and try improving integrating the incompatible sequences

							for seq_id in preferred_db.keys() :
								if not seq_id in best_2_edges_names :
									print >> sys.stderr, "##### Found equal optimal path without " + seq_id + ", trying to improve by reintegrating the incompatible sequences"
									#print >> sys.stderr, "###### To reintegrate: " + str(preferred_db[seq_id])
									all_preferred_used = False
									unwanted_to_reuse_list += three_orientation_list(seq_id , True)
									preferred_to_remove.append(seq_id)
									to_reintegrate += preferred_db[seq_id]

							#print >> sys.stderr, "to_reintegrate: " + str(to_reintegrate)
							if all_preferred_used :
								# Cannot not be improved as best uses all preferred >> best found
								best_found = True
							else :
								# Pathway do not uses all preferred, it is equal to one using preffered.
								# Unused preferred >> unwanted_to_remove_list (will not reduce the performance not using that sequence
								# Incompatible with unused preferred >> reintegrate to try to improve the pathway
								new_forced = best_2_edges_names
								to_reintegrate = list(set(to_reintegrate))
								for reintegrate_id in to_reintegrate :
									if (reintegrate_id in unwanted_to_remove_list) and (not reintegrate_id in unwanted_to_reuse_list) :
										unwanted_to_remove_list.remove(reintegrate_id)
										print >> sys.stderr, "###### Reintegrating: " + str(reintegrate_id)

								for seq_id in preferred_to_remove :
									del preferred_db[seq_id]
									unwanted_to_remove_list += three_orientation_list(seq_id , True)

								unwanted_to_remove_list = list(set(unwanted_to_remove_list))

						else :
							# unwanted_to_remove["keep"] = keep_seq_id|+
							# unwanted_to_remove["remove"] = [ reject_seq_id|+ , reject_seq_id|- , reject_seq_id|. ]
							print >> sys.stderr, "#### Found pair of sequences in the path that should not be placed in the same haplotype: " + str(unwanted_to_remove["keep"]) + " (to keep) | " + str(unwanted_to_remove["remove"])  + " (to remove)"

							# Check if the one we are going to remove was a preferred before, eventually reconsider the removal decision
							to_reintegrate = []
							if unwanted_to_remove["remove"][0] in preferred_db :
								to_reintegrate += preferred_db[unwanted_to_remove["remove"][0]]
								del preferred_db[unwanted_to_remove["remove"][0]]
							if unwanted_to_remove["remove"][1] in preferred_db :
								to_reintegrate += preferred_db[unwanted_to_remove["remove"][1]]
								del preferred_db[unwanted_to_remove["remove"][1]]
							if unwanted_to_remove["remove"][2] in preferred_db :
								to_reintegrate += preferred_db[unwanted_to_remove["remove"][2]]
								del preferred_db[unwanted_to_remove["remove"][2]]
								# The sequence to be removed was once preferred to other(s) sequence(s) in the path.
								# Reintegrate the deleted one(s), if they are compatible

							if to_reintegrate == [] :
								unwanted_to_remove_list += unwanted_to_remove["remove"]
							else :
								unwanted_to_remove_list = list(set(unwanted_to_remove_list))
								# Reintegrate the old sequences by removing them from the unwanted list
								for reintegrate_id in to_reintegrate :
									if (reintegrate_id in unwanted_to_remove_list) and (reintegrate_id not in unwanted_to_reuse_list) :
										unwanted_to_remove_list.remove(reintegrate_id)
										print >> sys.stderr, "###### Reintegrating: " + reintegrate_id
								unwanted_to_remove_list += unwanted_to_remove["remove"]
							unwanted_to_remove_list = list(set(unwanted_to_remove_list))

							if unwanted_to_remove["keep"] not in preferred_db :
								preferred_db[unwanted_to_remove["keep"]] = []
							preferred_db[unwanted_to_remove["keep"]] += unwanted_to_remove["remove"]

					print >> sys.stderr ,"### Best tiling @ 2nd iter (Haplotype 2): [" + "] -> [".join( ",".join(str(r) for r in x) for x in best_2_edges) + "]"
					new_best_2_paths[Target_sequence] = best_2_edges

					used , new_best_2_paths_edges[Target_sequence] = make_list_from_path( new_best_2_paths[Target_sequence] )
					map_used += used
					used_by_chr_hap1[Target_sequence] = used

					print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Used for " + Target_sequence + ": " + ",".join(new_best_2_paths_edges[Target_sequence])
					print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Used so far: " + ",".join(map_used)
		else :
			new_best_2_paths_edges = best_2_paths_edges
			new_best_2_paths = best_2_paths

		# Report sequences that were placed with markers but removed from path because not mapping
		# Parse best_1_paths_edges|best_1_paths_edges vs map_used -> return the seq_id before and after
		with_makers_unmapped_file = open( options.out + ".with_markers_not_mapping.txt" , 'w')
		print >> with_makers_unmapped_file , "\t".join(["component_id" , "orientation" ,  "chr_id" , "haplotype" , "prev_component_id" , "next_component_id" ])
		for chr_id in sorted(best_1_paths_edges.keys()) :
			best_1_paths_edges_list = best_1_paths_edges[chr_id]
			for i in range(0,len(best_1_paths_edges_list)) :
				prev = i-1
				next = i+1
				component_id = best_1_paths_edges_list[i]
				if component_id not in map_used :
					if prev<0 :
						prev_id = "."
					else :
						prev_id = best_1_paths_edges_list[prev]
					if next >= len(best_1_paths_edges_list) :
						next_id = "."
					else :
						next_id = best_1_paths_edges_list[next]
					print >> with_makers_unmapped_file , "\t".join([component_id[:-2] , component_id[-1] , chr_id , "hap1" , prev_id , next_id])
		for chr_id in sorted(best_2_paths_edges.keys()) :
			best_2_paths_edges_list = best_2_paths_edges[chr_id]
			for i in range(0,len(best_2_paths_edges_list)) :
				prev = i-1
				next = i+1
				component_id = best_2_paths_edges_list[i]
				if component_id not in map_used :
					if prev<0 :
						prev_id = "."
					else :
						prev_id = best_2_paths_edges_list[prev][:-2]
					if next >= len(best_2_paths_edges_list) :
						next_id = "."
					else :
						next_id = best_2_paths_edges_list[next][:-2]
					print >> with_makers_unmapped_file , "\t".join([component_id[:-2] , component_id[-1] , chr_id , "hap2" , prev_id , next_id])
		with_makers_unmapped_file.close()

		# Update all files from intermediate to final for postprocessing
		all_used = map_used
		best_1_paths = new_best_1_paths
		best_1_paths_edges = new_best_1_paths_edges
		best_2_paths = new_best_2_paths
		best_2_paths_edges = new_best_2_paths_edges

	# Tiling path creation completed.
	# data available:
	# all_used -> list of used sequences (both orientations reported) ([ "seq_id|+" , "seq_id|-" , ... ,])
	# best_1_paths
	# best_1_paths_edges  sorted list of used sequences
	# best_2_paths
	# best_2_paths_edges

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Find leftover query sequences"
	#### Extract list of unused sequences

	all_unused = []

	for id in sorted(query_len.keys()) :
		if str(id) + "|+" in all_used :
			# Query sequence, forward or reversed, used for assembling
			continue
		else :
			all_unused.append(id)

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == " + str(len(all_unused)) + " mapping sequences left unused"
	print >> sys.stderr ,"------------------------------"
	print >> sys.stderr ,"## Leftover sequences: "+ " , ".join(all_unused)
	print >> sys.stderr ,"------------------------------"

	#### Print output files

	gap_length = int(options.gap)

	#### Print lists

	list_1_file = open(options.out + ".1" + ".list" , "w")
	list_Un_file = open(options.out + ".Un" + ".list" , "w")
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Printing sorted list of sequences in tiling paths"
	for Target_sequence in sorted(best_1_paths_edges.keys()):
		print >> list_1_file , Target_sequence + "\t" + ",".join(best_1_paths_edges[Target_sequence])
	if not options.No2 :
		list_2_file = open(options.out + ".2" + ".list" , "w")
		for Target_sequence in sorted( best_2_paths_edges.keys() ):
			print >> list_2_file , Target_sequence + "\t" + ",".join(best_2_paths_edges[Target_sequence])
		list_2_file.close()
	print >> list_Un_file , ",".join(all_unused)
	list_1_file.close()
	list_Un_file.close()

	rejected_list_name = options.out + ".unused_sequences.list"
	rejected_list = open(rejected_list_name , 'w')
	for chr_id in sorted(unused_by_chr.keys()) :
		print >> rejected_list , chr_id + "\t" + ",".join(unused_by_chr[chr_id])
	rejected_list.close()


	##### Print path in AGP format
	if options.agp :
		### Query list element format
		# Gap:	[61252	,	(0:61252)		,	gap		,	61252	,	0]
		#		[length	,	(T_start:Tstop)	, 	"gap" 	, 	length 	, 	0]
		# Object:	[b40-14.HS_iter4_seq3733|+	,	(61252:6463804)	,	93612:7595148	,	-6402552			,	4526208]
		#			[ID|strand					,	(T_start:Tstop)	,	Q_start:Q_stop	,	-(alignment length)	,	matches]
		agp_1_file_name = options.out + ".1" + ".agp"
		agp_2_file_name = options.out + ".2" + ".agp"
		agp_Un_file_name = options.out + ".Un" + ".agp"
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Printing the tiling path in AGP file (gap size " + str(gap_length) + "bp)"
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap1"
	else :
		agp_1_file_name = "." + options.out + ".1" + ".agp"
		agp_Un_file_name = "." + options.out + ".Un" + ".agp"
		agp_2_file_name = "." + options.out + ".2" + ".agp"

	agp_1_file = open( agp_1_file_name, "w")
	agp_Un_file = open(agp_Un_file_name , "w")

	for Target_sequence in sorted(best_1_paths_edges.keys()):
		Obj_name = options.prefix + "_Hap1_" + Target_sequence
		hap1_ids.append(Obj_name)
		ref_to_hap1[Target_sequence] = Obj_name
		hap1_to_ref[Obj_name] = Target_sequence
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === " + Obj_name
		make_agp_from_list( best_1_paths_edges[Target_sequence] , query_len , gap_length , Obj_name , agp_1_file)
	agp_1_file.close()

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap2"

	if not options.No2 :
		agp_2_file = open( agp_2_file_name , "w")
		for Target_sequence in sorted(best_2_paths_edges.keys()):
			Obj_name = options.prefix + "_Hap2_" + Target_sequence
			hap2_ids.append(Obj_name)
			ref_to_hap2[Target_sequence] = Obj_name
			hap2_to_ref[Obj_name] = Target_sequence
			hap1_id = ref_to_hap1[Target_sequence]
			hap2_to_hap1[Obj_name] = hap1_id
			hap1_to_hap2[hap1_id] = Obj_name
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === " + Obj_name
			make_agp_from_list( best_2_paths_edges[Target_sequence] , query_len , gap_length , Obj_name , agp_2_file)
		agp_2_file.close()

	if not options.conc == "" :
		conc_gap = int(options.conc)
		Obj_name = options.prefix + "_Un"
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Unplaced sequences will be concatenated (gap size" + options.conc + "bp)"
		make_agp_from_list( all_unused , query_len , conc_gap , Obj_name , agp_Un_file)
	else :
		id = 0
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Unplaced sequences"
		for comp in all_unused :
			#	print >> sys.stderr , comp
			id += 1
			Obj_name = options.prefix + "_Un_" + str(id)
			make_agp_from_list( [comp + "|+"] , query_len , 0 , Obj_name , agp_Un_file)
	agp_Un_file.close()

	agp_db_1 = read_agp(agp_1_file_name)
	agp_db = dict(agp_db_1)
	if not options.No2 :
		agp_db_2 = read_agp(agp_2_file_name)
		agp_db.update(agp_db_2)
	agp_db_un = read_agp(agp_Un_file_name)
	agp_db.update(agp_db_un)

	# Translate marker coordinates
	if options.markers_hits :
		bed_regions = read_bed_sorted_list(options.markers_hits)
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Translating marker coordinates"
		new_bed = translate_bed_sorted_list( bed_regions , agp_db )
		out_bed_file_name = options.out + ".markers.bed"
		out_bed_file = open(out_bed_file_name , 'w')
		for line in sorted(new_bed) :
			print >> out_bed_file, "\t".join([str(x) for x in line])
		out_bed_file.close()

	# Port legacy AGP information
	if options.input_agp :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Translating input sequence structure from AGP file " + options.input_agp
		old_legacy_agp = read_agp(options.input_agp)
		legacy_agp = agp_translate_agp( agp_db , old_legacy_agp )
		legacy_agp_file_name = options.out + ".legacy_structure.agp"
		legacy_agp_file_name = write_agp( legacy_agp , legacy_agp_file_name )
	else :
		old_legacy_agp = ""
		legacy_agp = ""

	##### Generate fasta sequences
	fasta_db_1 = {}
	fasta_1_len = {}
	chr_to_fasta_1 = {}
	fasta_to_chr = {}
	fasta_chr_1 = {}
	fasta_db_2 = {}
	fasta_2_len = {}
	fasta_chr_2 = {}
	chr_to_fasta_2 = {}
	for Target_sequence in sorted(best_1_paths.keys()):
		Obj_name = options.prefix + "_Hap1_" + Target_sequence
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === " + Obj_name
		fasta_db_1[Obj_name] = make_fasta_from_list( best_1_paths_edges[Target_sequence] , query_fasta_db ,  gap_length , Obj_name )
		fasta_1_len[Obj_name] = len(fasta_db_1[Obj_name])
		fasta_chr_1[Obj_name] = Target_sequence
		fasta_to_chr[Obj_name] = Target_sequence
		chr_to_fasta_1[Target_sequence] = Obj_name
	if not options.No2 :
		for Target_sequence in sorted(best_2_paths.keys()):
			Obj_name = options.prefix + "_Hap2_" + Target_sequence
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === " + Obj_name
			fasta_db_2[Obj_name] = make_fasta_from_list( best_2_paths_edges[Target_sequence] , query_fasta_db ,  gap_length , Obj_name)
			fasta_2_len[Obj_name] = len(fasta_db_2[Obj_name])
			fasta_chr_2[Obj_name] = Target_sequence
			fasta_to_chr[Obj_name] = Target_sequence
			chr_to_fasta_2[Target_sequence] = Obj_name
	fasta_db_un = {}
	fasta_un_len = {}
	if not options.conc == "" :
		conc_gap = int(options.conc)
		Obj_name = options.prefix + "_Un"
		all_unused_path = [ [x+"|+" , 0 ,0 ,0 ,0] for x in all_unused ]
		fasta_db_un[Obj_name] = make_fasta_from_list( all_unused , query_fasta_db , conc_gap , Obj_name)
		fasta_un_len[Obj_name] = fasta_db_un[Obj_name]
	else :
		id = 0
		for comp in all_unused :
			id += 1
			Obj_name = options.prefix + "_Un_" + str(id)
			fasta_db_un[Obj_name] = make_fasta_from_list( [ comp+"|+" ] , query_fasta_db , 0 , Obj_name)
			fasta_un_len[Obj_name] = fasta_db_un[Obj_name]



	# Perform relationship check pseudomolecule vs. pseudomolecule and pseudomolecule vs. unplaced
	# Report input sequence groups to pseudomolecules and add them to input information (known_groups >> groups_by_seqid)
	# search for alternative sequences to pseudomolecule (alternative_sequences)
	groups_by_seqid = {}
	for group_id in known_groups.keys() :
		for seq_id in known_groups[group_id] :
			if seq_id not in groups_by_seqid :
				groups_by_seqid[seq_id] = []
			groups_by_seqid[seq_id].append(group_id)

	groups_by_pseudomolecule = {}
	pseudomolecule_by_group = {}
	alternative_by_pseudomolecule = {}
	pseudomolecule_by_alternative = {}
	for Target_sequence in sorted(best_1_paths_edges.keys()):
		new_id = options.prefix + "_Hap1_" + Target_sequence
		groups_by_pseudomolecule[new_id] = []
		for seq_id in best_1_paths_edges[Target_sequence] :
			clean_name = seq_id[:-2]
			if clean_name in groups_by_seqid :
				seq_id_groups = groups_by_seqid[clean_name]
				groups_by_pseudomolecule[new_id] += seq_id_groups
				for group_id in seq_id_groups :
					if group_id not in pseudomolecule_by_group :
						pseudomolecule_by_group[group_id] = []
					pseudomolecule_by_group[group_id].append(new_id)
			if clean_name in alternative_sequences :
				alternative_ids = alternative_sequences[clean_name]
				alternative_by_pseudomolecule[new_id] = alternative_ids
				alternative_by_pseudomolecule[new_id] = list(set(alternative_by_pseudomolecule[new_id]))
				for alt_id in alternative_ids :
					if alt_id not in pseudomolecule_by_alternative :
						pseudomolecule_by_alternative[alt_id] = []
					pseudomolecule_by_alternative[alt_id].append(new_id)
					pseudomolecule_by_alternative[alt_id] = list(set(pseudomolecule_by_alternative[alt_id]))

	for Target_sequence in sorted(best_2_paths_edges.keys()):
		new_id = options.prefix + "_Hap2_" + Target_sequence
		groups_by_pseudomolecule[new_id] = []
		for seq_id in best_2_paths_edges[Target_sequence] :
			clean_name = seq_id[:-2]
			if clean_name in groups_by_seqid :
				seq_id_groups = groups_by_seqid[clean_name]
				groups_by_pseudomolecule[new_id] += seq_id_groups
				for group_id in seq_id_groups :
					if group_id not in pseudomolecule_by_group :
						pseudomolecule_by_group[group_id] = []
					pseudomolecule_by_group[group_id].append(new_id)
			if clean_name in alternative_sequences :
				alternative_ids = alternative_sequences[clean_name]
				alternative_by_pseudomolecule[new_id] = alternative_ids
				alternative_by_pseudomolecule[new_id] = list(set(alternative_by_pseudomolecule[new_id]))
				for alt_id in alternative_ids :
					if alt_id not in pseudomolecule_by_alternative :
						pseudomolecule_by_alternative[alt_id] = []
					pseudomolecule_by_alternative[alt_id].append(new_id)
					pseudomolecule_by_alternative[alt_id] = list(set(pseudomolecule_by_alternative[alt_id]))

	#json.dump( groups_by_pseudomolecule , open("groups_by_pseudomolecule.json" , 'w') )
	#json.dump( pseudomolecule_by_group , open("pseudomolecule_by_group.json" , 'w') )

	# Pseudomolecules vs Pseudomolecules
	conflicting_pseudomolecules_file = open( options.out + ".conflicting_pseudomolecules.txt" , 'w')
	print >> conflicting_pseudomolecules_file , "Chr_id\tConflicting_pseudomolecules"
	for Target_sequence in sorted(groups_by_pseudomolecule.keys()) :
		pseudomolecule_groups = groups_by_pseudomolecule[Target_sequence]
		conflicting_pseudomolecules = []
		for group_id in pseudomolecule_groups :
			conflicting_pseudomolecules += pseudomolecule_by_group[group_id]
		conflicting_pseudomolecules = list(set(conflicting_pseudomolecules))
		if Target_sequence in conflicting_pseudomolecules :
			conflicting_pseudomolecules.remove(Target_sequence)
		print >> conflicting_pseudomolecules_file , Target_sequence + "\t" + ",".join(conflicting_pseudomolecules)
	conflicting_pseudomolecules_file.close()


	# unused_by_chr[chr_id] -> dict on seqids
	unused_by_seq_id = {}
	all_unused_by_seq_id = {}
	for chr_id in unused_by_chr :
		unused_list = unused_by_chr[chr_id]
		for seq_id_oriented in unused_list :
			seq_id , orientation = seq_id_oriented.split("|")
			unused_by_seq_id[seq_id] = [chr_id , orientation]
			all_unused_by_seq_id[seq_id] = [chr_id , orientation]
	
	conflicting_unplaced_file = open( options.out + ".unplaced_to_pseudomolecule.txt" , 'w')
	print >> conflicting_unplaced_file , "Seq_id\tConflicting_marker_Vs_associated\tMultiple_associated_pseudomolecules\tConflicting_marker_Vs_alternative\tMultiple_alternative_pseudomolecules\tExpected_chr\tOrientation\tAssociated_pseudomolecules\tAlternative_pseudomolecules"
	for seq_id in all_unused :
		if seq_id in unused_by_seq_id :
			chr_id , orientation = unused_by_seq_id[seq_id]
		else :
			chr_id , orientation = [ "." , "."]

		associated_pseudomolecules = []
		if seq_id in groups_by_seqid :
			groups = groups_by_seqid[seq_id]
			for group_id in groups :
				if group_id in pseudomolecule_by_group :
					associated_pseudomolecules += pseudomolecule_by_group[group_id]
			associated_pseudomolecules = list(set(associated_pseudomolecules))

		if len(associated_pseudomolecules) == 0 :
			err_chr_assoc = "."
			mult_assoc = "."
			assoc_chromosomes = [ ["."] ]
		elif len(associated_pseudomolecules) > 1 :
			mult_assoc = "TRUE"
			assoc_chromosomes = list(set([ fasta_to_chr[x] for x in associated_pseudomolecules ]))
			if len(assoc_chromosomes) > 1 :
				err_chr_assoc = "TRUE"
			else :
				if assoc_chromosomes[0] == chr_id:
					err_chr_assoc = "FALSE"
				else :
					if chr_id == "." :
						err_chr_assoc = "."
					else :
						err_chr_assoc = "TRUE"
		else :
			#len == 1
			mult_assoc = "FALSE"
			assoc_chromosomes = [ fasta_to_chr[associated_pseudomolecules[0]] ]
			if assoc_chromosomes == [ chr_id ] :
				err_chr_assoc = "FALSE"
			else :
				if chr_id == "." :
					err_chr_assoc = "."
				else :
					err_chr_assoc = "TRUE"

		alternative_pseudomolecules = []
		if seq_id in pseudomolecule_by_alternative :
			alternative_pseudomolecules = list(set(pseudomolecule_by_alternative[seq_id]))

		if len(alternative_pseudomolecules) == 0 :
			err_chr_alt = "."
			mult_alt = "."
			alt_chromosomes = [ ["."] ]
		elif len(alternative_pseudomolecules) > 1 :
			mult_alt = "TRUE"
			alt_chromosomes = list(set([ fasta_to_chr[x] for x in alternative_pseudomolecules ]))
			if len(alt_chromosomes) > 1 :
				err_chr_alt = "TRUE"
			else :
				if alt_chromosomes[0] == chr_id:
					err_chr_alt = "FALSE"
				else :
					if chr_id == "." :
						err_chr_alt = "."
					else :
						err_chr_alt = "TRUE"
		else :
			mult_alt = "FALSE"
			alt_chromosomes = [ fasta_to_chr[alternative_pseudomolecules[0]] ]
			if alt_chromosomes == [ chr_id ] :
				err_chr_alt = "FALSE"
			else :
				if chr_id == "." :
					err_chr_alt = "."
				else :
					err_chr_alt = "TRUE"

		print >> conflicting_unplaced_file , seq_id + "\t" + err_chr_assoc + "\t" + mult_assoc + "\t" + err_chr_alt + "\t" + mult_alt + "\t" + chr_id + "\t" + orientation + "\t" + ",".join(associated_pseudomolecules) + "\t" + ",".join(alternative_pseudomolecules)

		if seq_id not in all_unused_by_seq_id :
			# no chr assigned based on markers
			#print >> conflicting_unplaced_file , "assoc_chromosomes:" + str(assoc_chromosomes)
			#print >> conflicting_unplaced_file , "alt_chromosomes:" + str(alt_chromosomes)

			if ( len(assoc_chromosomes) == 1 ) and ( len(alt_chromosomes) == 1) :
				#print >> conflicting_unplaced_file , "1, only 1 chr"
				# alternative/associated pseudomolecules info lead to one chr only
				if assoc_chromosomes[0] == alt_chromosomes[0] :
					#print >> conflicting_unplaced_file , "2, equal"
					# alternative/associated pseudomolecules info lead to the same chromosome
					if not assoc_chromosomes[0] == ["."] :
						#print >> conflicting_unplaced_file , "3, known"
						# valid chr, add the info of chr association
						#print >> conflicting_unplaced_file , "*"
						all_unused_by_seq_id[seq_id] = [assoc_chromosomes[0] , "+"]
				else :
					#print >> conflicting_unplaced_file , "2, different"
					if (not assoc_chromosomes[0] == ["."] ) and ( alt_chromosomes[0] == ["."] ) :
						#print >> conflicting_unplaced_file , "3, 1 known"
						#print >> conflicting_unplaced_file , "*"
						all_unused_by_seq_id[seq_id] = [assoc_chromosomes[0] , "+"]
					elif (assoc_chromosomes[0] == ["."] ) and (not alt_chromosomes[0] == ".") :
						#print >> conflicting_unplaced_file , "3, 1 known"
						#print >> conflicting_unplaced_file , "*"
						all_unused_by_seq_id[seq_id] = [alt_chromosomes[0] , "+"]
					#else :
					#	print >> conflicting_unplaced_file , "3, different chr"
			#else :
			#	print >> conflicting_unplaced_file , "1, multiple chr"

	conflicting_unplaced_file.close()

	#### Print sequence FASTA
	if options.sequence :
		fasta_1_file = options.out + ".1" + ".fasta"
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Printing sequences of tiling path in fasta files (gap size " + str(gap_length) + "bp)"
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap1"
		write_fasta_from_db( fasta_db_1 , fasta_1_file )

		if not options.No2 :
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap2"
			fasta_2_file = options.out + ".2" + ".fasta"
			write_fasta_from_db( fasta_db_2 , fasta_2_file )

		fasta_Un_file = options.out + ".Un" + ".fasta"
		if not options.conc == "" :
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Unplaced sequences will be concatenated (gap size" + options.conc + "bp)"
		else :
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Unplaced sequences"
		write_fasta_from_db( fasta_db_un , fasta_Un_file )

	#### Convert annotation, if given
	if options.gff3 :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Converting annotation"

		#annotation_dir = options.out + ".dedup_dir"
		#mkdir( "./" + annotation_dir )

		### Read AGP files and conver them to offset values
		translation_db = translate_from_AGP_whole_genome(read_agp(options.out + ".1" + ".agp"))
		if not options.No2 :
			translation_db.update(translate_from_AGP_whole_genome(read_agp(options.out + ".2" + ".agp")))
		translation_db.update(translate_from_AGP_whole_genome(read_agp(options.out + ".Un" + ".agp")))

		### Convert coordinates
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Converting coordinates"
		new_gff3 = translate_gff3(annotation_gff3 , translation_db , options.out + ".broken_genes.txt" )

		### Write GFF3
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Writing update GFF3"
		output_seq_lengths = get_fasta_lengths_from_file(options.out + ".1" + ".fasta")
		if not options.No2 :
			output_seq_lengths.update(get_fasta_lengths_from_file(options.out + ".2" + ".fasta"))
		output_seq_lengths.update(get_fasta_lengths_from_file(options.out + ".Un" + ".fasta"))

		write_gff3( new_gff3 , options.out + ".annotation.gff3" , output_seq_lengths )
		mRNA_to_gene_db = get_gene2mRNA(options.out + ".annotation.gff3")


	if options.haplodup :
	#### Generate HaploDup reports and plots
		coord_tables = {}
		plot_files = {}
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Performing HaploDup"
		haplodup_dir = options.out + ".HaploDup_dir"
		mkdir( "./" + haplodup_dir )

		# Convert agp structure in table
		structure_file = "structure.tsv"
		structure_file_fullpath = haplodup_dir + "/structure.tsv"
		agp_table = []
		for seq_id in agp_db :
			seq_agp = agp_db[seq_id]
			for start in sorted(seq_agp.keys()) :
				Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = seq_agp[start]
				if Compnt_Type == "W" :
					agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation ])
		structure_file_fullpath = write_table(agp_table, structure_file_fullpath)

		markers_db = {}
		# Prepare files from markers and agps
		if options.markers_hits :
			# Copy the file in haplodup_dir
			all_markers_file_name = "all_markers.bed"
			dup_markers_file_name = "duplicated_markers.bed"
			all_markers_file_name_fullpath = haplodup_dir + "/all_markers.bed"
			dup_markers_file_name_fullpath = haplodup_dir + "/duplicated_markers.bed"
			seq_all_markers_file = open(all_markers_file_name_fullpath , 'w')
			seq_duplicated_markers_file = open(dup_markers_file_name_fullpath , 'w')
			## Write markers files
			for line in open(options.out + ".markers.bed") :
				chr_id , start , stop , marker_id = line.rstrip().split("\t")
				if chr_id not in markers_db :
					markers_db[chr_id] = {}
				if marker_id not in markers_db[chr_id] :
					markers_db[chr_id][marker_id] = []
				markers_db[chr_id][marker_id].append([chr_id , start , stop , marker_id])

			# Read translated marker coordinates
			for chr_id in sorted(markers_db.keys()) :
				for marker_id in sorted(markers_db[chr_id].keys()) :
					if len(markers_db[chr_id][marker_id]) == 1 :
						# Unique hit
						print >> seq_all_markers_file, "\t".join(markers_db[chr_id][marker_id][0])
					else :
						# Multiple hits, report all in both files
						for hit in sorted(markers_db[chr_id][marker_id]) :
							print >> seq_all_markers_file, "\t".join(hit)
							print >> seq_duplicated_markers_file, "\t".join(hit)
			seq_all_markers_file.close()
			seq_duplicated_markers_file.close()
		else :
			all_markers_file_name = ""
			all_markers_file_name_fullpath = ""
			dup_markers_file_name = ""
			dup_markers_file_name_fullpath = ""

		if not legacy_agp == "" :
			# convert the db in table and write it down
			legacy_structure_file = "legacy_structure.tsv"
			legacy_structure_file_fullpath = haplodup_dir + "/legacy_structure.tsv"
			legacy_agp_table = []
			for seq_id in legacy_agp :
				seq_agp = legacy_agp[seq_id]
				for start in sorted(seq_agp.keys()) :
					Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = seq_agp[start]
					if Compnt_Type == "W" :
						legacy_agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation ])
			legacy_structure_file_fullpath = write_table(legacy_agp_table, legacy_structure_file_fullpath)
		else :
			legacy_structure_file = ""

		# Write ID correspondence table
		# Tab separated, Columns: Hap1	[Hap2]	Ref

		corr_file_name = "correspondence.tsv"
		corr_file = open( haplodup_dir + "/" + corr_file_name , 'w')
		if not options.No2 :
			# Include Hap2 IDs
			for id in sorted(hap1_to_ref.keys()):
				print >> corr_file , str(id) + "\t" + str(hap1_to_hap2[id]) + "\t" + str(hap1_to_ref[id])
		else :
			for id in sorted(hap1_to_ref.keys()):
				print >> corr_file , str(id) + "\t" + str(hap1_to_ref[id])
		corr_file.close()

		showcoords_path = paths["show-coords"]
		if showcoords_path == "" :
			minimap2_search=subprocess.Popen( "which show-coords" , shell=True, stdout=subprocess.PIPE )
			command_line , error = minimap2_search.communicate()
			command_line = command_line.rstrip()
		else :
			command_line = showcoords_path + "/show-coords"
		if not os.path.exists(command_line) :
			print >> sys.stderr , "[ERROR] wrong or no path to show-coords (Mummer4)"
			sys.exit(1)

		# Perform nucmer alignments
		hap1_ids = ",".join(sorted(fasta_db_1.keys()))
		if options.reference :
			reference_ids = ",".join(sorted(reference.keys()))
		if not options.No2 :
			hap2_ids = ",".join(sorted(fasta_db_2.keys()))

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Mapping coordinates"
		query_1_file = options.out + ".1.fasta"

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Hap1 vs Hap1"
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Mapping"
		if options.debug :
			outfile_prefix = "Hap1.on.Hap1"
		else :
			outfile_prefix = map_nucmer_dotplot("Hap1" , query_1_file , "Hap1" , query_1_file , haplodup_dir , options.cores , paths , False )
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Converting files"
		coord_tables["Hap1_vs_Hap1"] = make_coords_table( outfile_prefix , haplodup_dir , command_line)
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Generating dotplot"
		plot_files["Hap1_vs_Hap1"] = {}
		if options.debug :
			plot_files["Hap1_vs_Hap1"]["Whole"] , plot_files["Hap1_vs_Hap1"]["All_Dotplots"] = whole_genome_dotplot( hap1_ids , hap1_ids , outfile_prefix, haplodup_dir, "Hap1_vs_Hap1", coord_tables["Hap1_vs_Hap1"] , True)
		else :
			plot_files["Hap1_vs_Hap1"]["Whole"] , plot_files["Hap1_vs_Hap1"]["All_Dotplots"] = whole_genome_dotplot( hap1_ids , hap1_ids , outfile_prefix, haplodup_dir, "Hap1_vs_Hap1", coord_tables["Hap1_vs_Hap1"])

		if not options.No2 :
			query_2_file = options.out + ".2.fasta"

			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Hap2 vs Hap2"
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Mapping"
			if options.debug :
				outfile_prefix = "Hap2.on.Hap2"
			else :
				outfile_prefix = map_nucmer_dotplot("Hap2" , query_2_file , "Hap2" , query_2_file , haplodup_dir , options.cores , paths , False  )
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Converting files"
			coord_tables["Hap2_vs_Hap2"] = make_coords_table( outfile_prefix , haplodup_dir , command_line)
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Generating dotplot"
			plot_files["Hap2_vs_Hap2"] = {}
			if options.debug :
				plot_files["Hap2_vs_Hap2"]["Whole"] , plot_files["Hap2_vs_Hap2"]["All_Dotplots"]  = whole_genome_dotplot( hap2_ids, hap2_ids , outfile_prefix, haplodup_dir, "Hap2_vs_Hap2", coord_tables["Hap2_vs_Hap2"] , True)
			else :
				plot_files["Hap2_vs_Hap2"]["Whole"] , plot_files["Hap2_vs_Hap2"]["All_Dotplots"]  = whole_genome_dotplot( hap2_ids, hap2_ids , outfile_prefix, haplodup_dir, "Hap2_vs_Hap2", coord_tables["Hap2_vs_Hap2"])

			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Hap2 vs Hap1"
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Mapping"
			if options.debug :
				outfile_prefix = "Hap2.on.Hap1"
			else :
				outfile_prefix = map_nucmer_dotplot("Hap1" , query_1_file , "Hap2" , query_2_file , haplodup_dir , options.cores , paths , False  )
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Converting files"
			coord_tables["Hap2_vs_Hap1"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , coord_tables["Hap1_vs_Hap1"] ]
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Generating dotplot"
			plot_files["Hap2_vs_Hap1"] = {}
			if options.debug :
				plot_files["Hap2_vs_Hap1"]["Whole"] , plot_files["Hap2_vs_Hap1"]["All_Dotplots"] = whole_genome_dotplot( hap1_ids , hap2_ids , outfile_prefix, haplodup_dir, "Hap2_vs_Hap1", coord_tables["Hap2_vs_Hap1"][0] , True)
			else :
				plot_files["Hap2_vs_Hap1"]["Whole"] , plot_files["Hap2_vs_Hap1"]["All_Dotplots"] = whole_genome_dotplot( hap1_ids , hap2_ids , outfile_prefix, haplodup_dir, "Hap2_vs_Hap1", coord_tables["Hap2_vs_Hap1"][0])

			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Hap1 vs Hap2"
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Mapping"
			if options.debug :
				outfile_prefix = "Hap1.on.Hap2"
			else :
				outfile_prefix = map_nucmer_dotplot("Hap2" , query_2_file , "Hap1" , query_1_file , haplodup_dir , options.cores , paths , False  )
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Converting files"
			coord_tables["Hap1_vs_Hap2"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , coord_tables["Hap2_vs_Hap2"] ]
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Generating dotplot"
			plot_files["Hap1_vs_Hap2"] = {}
			if options.debug :
				plot_files["Hap1_vs_Hap2"]["Whole"], plot_files["Hap1_vs_Hap2"]["All_Dotplots"] = whole_genome_dotplot( hap2_ids , hap1_ids, outfile_prefix, haplodup_dir, "Hap1_vs_Hap2", coord_tables["Hap1_vs_Hap2"][0] , True)
			else :
				plot_files["Hap1_vs_Hap2"]["Whole"], plot_files["Hap1_vs_Hap2"]["All_Dotplots"] = whole_genome_dotplot( hap2_ids , hap1_ids, outfile_prefix, haplodup_dir, "Hap1_vs_Hap2", coord_tables["Hap1_vs_Hap2"][0])


		if options.reference :
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Hap1 vs Reference"
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Mapping"
			if options.debug :
				outfile_prefix = "Hap1.on.Ref"
			else :
				outfile_prefix = map_nucmer_dotplot("Ref" , options.reference , "Hap1" , query_1_file , haplodup_dir , options.cores , paths , False )
				# outfile_prefix.delta and outfile_prefix.coords (show-cords -c) are generated
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Converting files"
			coord_tables["Hap1_vs_Reference"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , "" ]
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Generating dotplot"
			plot_files["Hap1_vs_Reference"] = {}
			if options.debug :
				plot_files["Hap1_vs_Reference"]["Whole"] , plot_files["Hap1_vs_Reference"]["All_Dotplots"] = whole_genome_dotplot( reference_ids , hap1_ids , outfile_prefix, haplodup_dir, "Hap1_vs_Reference", coord_tables["Hap1_vs_Reference"][0] , True)
			else :
				plot_files["Hap1_vs_Reference"]["Whole"] , plot_files["Hap1_vs_Reference"]["All_Dotplots"] = whole_genome_dotplot( reference_ids , hap1_ids , outfile_prefix, haplodup_dir, "Hap1_vs_Reference", coord_tables["Hap1_vs_Reference"][0])

			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Reference vs Hap1 "
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Mapping"
			if options.debug :
				outfile_prefix = "Ref.on.Hap1"
			if options.debug :
				outfile_prefix = map_nucmer_dotplot( "Hap1" , query_1_file , "Ref" , options.reference , haplodup_dir , options.cores , paths , False )
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Converting files"
			coord_tables["Reference_vs_Hap1"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , coord_tables["Hap1_vs_Hap1"] ]
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Generating dotplot"
			plot_files["Reference_vs_Hap1"] = {}
			if options.debug :
				plot_files["Reference_vs_Hap1"]["Whole"] , plot_files["Reference_vs_Hap1"]["All_Dotplots"] = whole_genome_dotplot( hap1_ids , reference_ids , outfile_prefix, haplodup_dir, "Reference_vs_Hap1", coord_tables["Reference_vs_Hap1"][0] , True)
			else :
				plot_files["Reference_vs_Hap1"]["Whole"] , plot_files["Reference_vs_Hap1"]["All_Dotplots"] = whole_genome_dotplot( hap1_ids , reference_ids , outfile_prefix, haplodup_dir, "Reference_vs_Hap1", coord_tables["Reference_vs_Hap1"][0])

		if not options.No2 :
			if options.reference :
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Hap2 vs Reference"
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Mapping"
				if options.debug :
					outfile_prefix = "Hap2.on.Ref"
				else:
					outfile_prefix = map_nucmer_dotplot("Ref" , options.reference , "Hap2" , query_2_file , haplodup_dir , options.cores , paths , False )
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Converting files"
				coord_tables["Hap2_vs_Reference"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , "" ]
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ===== Generating dotplot"
				plot_files["Hap2_vs_Reference"] = {}
				if options.debug :
					plot_files["Hap2_vs_Reference"]["Whole"] , plot_files["Hap2_vs_Reference"]["All_Dotplots"] = whole_genome_dotplot( reference_ids , hap2_ids , outfile_prefix, haplodup_dir, "Hap2_vs_Reference", coord_tables["Hap2_vs_Reference"][0]  , True)
				else :
					plot_files["Hap2_vs_Reference"]["Whole"] , plot_files["Hap2_vs_Reference"]["All_Dotplots"] = whole_genome_dotplot( reference_ids , hap2_ids , outfile_prefix, haplodup_dir, "Hap2_vs_Reference", coord_tables["Hap2_vs_Reference"][0])

				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Reference vs Hap2"
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Mapping"
				if options.debug :
					outfile_prefix = "Ref.on.Hap2"
				else:
					outfile_prefix = map_nucmer_dotplot("Hap2" , query_2_file , "Ref" , options.reference , haplodup_dir , options.cores , paths , False )
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Converting files"
				coord_tables["Reference_vs_Hap2"] = [make_coords_table( outfile_prefix , haplodup_dir , command_line) , coord_tables["Hap2_vs_Hap2"] ]
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ===== Generating dotplot"
				plot_files["Reference_vs_Hap2"] = {}
				if options.debug :
					plot_files["Reference_vs_Hap2"]["Whole"] , plot_files["Reference_vs_Hap2"]["All_Dotplots"] = whole_genome_dotplot( hap2_ids , reference_ids , outfile_prefix, haplodup_dir, "Reference_vs_Hap2", coord_tables["Reference_vs_Hap2"][0] , True)
				else :
					plot_files["Reference_vs_Hap2"]["Whole"] , plot_files["Reference_vs_Hap2"]["All_Dotplots"] = whole_genome_dotplot( hap2_ids , reference_ids , outfile_prefix, haplodup_dir, "Reference_vs_Hap2", coord_tables["Reference_vs_Hap2"][0])

		#### QC unused sequences for each chromosome build
		if not options.No2 :
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Structure comparison analysis"
			print >> sys.stderr , "# Structure comparison analysis"
			structure_comparison_dir = options.out + ".structure_comparison"
			mkdir(structure_comparison_dir)
			associated_input_seqid = []
			associated_input_seqid_file = ""
			associated_legacy_seqid = []
			associated_legacy_seqid_file = ""
			# Note: Coordinates are on pseudomolecules for hap1 and hap2 regions, on HS for "Query_legacy" == begin-to-end
			# Note: Where group_id is taken from the grouping file and defines the univocally the cluster >> color coding
			if options.legacy_groups :
				associated_legacy_seqid_file = structure_comparison_dir + "/legacy_components.association.tsv"
				# all_agp_db = legacy_agp + old_legacy_agp >> all sequences related to the components >> grouping on components
				all_agp_db = dict(legacy_agp)
				old_legacy_agp = read_agp(options.input_agp)
				all_agp_db.update(old_legacy_agp)
				associated_input_seqid_file , seq_group_db = make_seq_pair_from_groups( associated_legacy_seqid_file , options.legacy_groups , all_agp_db , hap1_ids.split(",") , hap2_ids.split(",") , query_fasta_db.keys() , {} , "legacy" )
				# associated_input_seqid_file >> [ ... , [ "hap1_legacy" , Tid , Tstart , Tstop , "Query_legacy" , Qid , Qstart , Qstop , group_id] , ... ]
			else :
				# No legacy info available -> use what is known the input sequences
				associated_input_seqid_file = structure_comparison_dir + "/input.association.tsv"
				if options.input_groups :
					# override with user option options.input_groups
					if options.input_groups == 0 :
						# force mapping,
						associated_input_seqid_file = ""
					else :
						# file name is given -> read as [ id , group ] table
						# agp_db + feed all input sequences length >> whole sequence relationship for unplaced >> grouping on input sequences
						associated_input_seqid_file , seq_group_db = make_seq_pair_from_groups( associated_input_seqid_file ,options.input_groups , agp_db , hap1_ids.split(",") , hap2_ids.split(",") , query_fasta_db.keys() , query_len , "input" )
						# associated_input_seqid_file >> [ ... , [ "hap1_HS" , Tid , Tstart , Tstop , "Query_HS" , Qid , Qstart , Qstop , group_id] , ... ]
				#else :
				#	# parse known_groups and unwanted_pairs
				#	associated_input_seqid_file = make_seq_pair_from_constrains(associated_input_seqid_file, known_groups, unwanted_pairs, alternative_sequences , agp_db, hap1_ids.split(","), hap2_ids.split(","), query_fasta_db.keys(), query_len, "input")
				#	# associated_input_seqid_file >> [ ... , [ "hap1_HS" , Tid , Tstart , Tstop , "Query_HS" , Qid , Qstart , Qstop , group_id] , ... ]

			# Merge all positions on input sequences and output sequences
			# Expected dict format: markers_db[seq_id][marker_id] = [ ... , [seq_id , start , stop , marker_id] , ... ]
			all_markers_db = dict(markers_db)
			# add marker_hits_by_seq info
			# marker_hits_by_seq[seq_id] = [ ... , [ int(start) , int(stop) , marker_id , marker_chr , int(marker_pos) ] , ... ]
			for seq_id in marker_hits_by_seq.keys() :
				marker_hits= marker_hits_by_seq[seq_id]
				if seq_id not in all_markers_db:
					all_markers_db[seq_id] = {}
				for element in marker_hits :
					start , stop , marker_id , marker_chr , marker_pos = element
					if marker_id not in all_markers_db[seq_id] :
						all_markers_db[seq_id][marker_id] = []
					all_markers_db[seq_id][marker_id].append([ seq_id , start , stop , marker_id ])


			# TODO: Perform Hap1 vs Hap2 pairwise comparison of known features (without unplaced)
			# Activate creation of associations between sequences using constrains
			# if not associated_input_seqid_file == "" :

			if (not options.avoidrejectedqc) and (not associated_input_seqid_file == "" ) :
				if (not options.only_markers) or (not options.markers_hits) :
					unused_to_check = dict(all_unused_by_seq_id)
					print >> sys.stderr , "## Unplaced sequences to check: " + str(len(unused_to_check.keys()))
				else :
					unused_to_check = {}
					for seq_id in all_unused_by_seq_id.keys() :
						if seq_id in marker_hits_by_seq :
							unused_to_check[seq_id] = all_unused_by_seq_id[seq_id]

					print >> sys.stderr , "## Unplaced sequences with markers to check: " + str(len(unused_to_check.keys())) + " out of " + str(len(all_unused_by_seq_id.keys())) + " unplaced in total"

				all_unused_by_chr = {}
				for seq_id in unused_to_check.keys() :
					chr_id , orientation = unused_to_check[seq_id]
					if chr_id not in all_unused_by_chr :
						all_unused_by_chr[chr_id] = []
					all_unused_by_chr[chr_id].append(seq_id+"|"+orientation)

				if not options.markers_hits :
					clean_marker_set_by_seq = {}
				structure_plot_db = {"Rejected" : {}}
				for hap1_id in sorted(fasta_db_1.keys()) :
					chr_id = fasta_chr_1[hap1_id]
					print >> sys.stdout , '[' + str(datetime.datetime.now()) + "] == Chr: " + chr_id
					print >> sys.stderr , "## Chr: " + chr_id
					if not chr_id in all_unused_by_chr :
						continue
					else :
						# DO QC
						qc_out_dir = structure_comparison_dir + "/" + chr_id
						mkdir(qc_out_dir)
						structure_plot_db["Rejected"][chr_id] = {}
						for seq_id in all_unused_by_chr[chr_id] :
							print >> sys.stderr , "### seq_id: " + seq_id
							structure_plot_db["Rejected"][chr_id][seq_id] = rejected_QC( structure_comparison_dir , seq_id , query_fasta_db , chr_id , fasta_db_1 , fasta_db_2 , chr_to_fasta_1 , chr_to_fasta_2 , coord_tables["Hap2_vs_Hap1"][0] , associated_input_seqid_file , associated_legacy_seqid_file , agp_db , legacy_agp , old_legacy_agp , seq_group_db , all_markers_db , clean_marker_set_by_seq, marker_map_by_seq , options.cores , paths)

				rejected_index_file_name = "index.rejected_sequences.html"
				rejected_index_file_full_path = structure_comparison_dir + "/" + rejected_index_file_name
				rejected_index_file_full_path = make_index_from_report_db(rejected_index_file_name , "." , structure_comparison_dir ,  structure_plot_db  )


		if options.gff3 and not options.No2:
			# Map genes, generate copy number counts, render by chromosome reports
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Generating files about gene map count for fusion dedup"
			print >> sys.stderr , "# Generating files about gene map count for fusion dedup"
			haplodup_dir = options.out + ".HaploDup_dir"

			## Generate CDS sequences from first haplotype
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Generating CDS sequences"
			print >> sys.stderr , "## Generating CDS sequences"
			new_Seq_db = read_fasta(options.out + ".1" + ".fasta")
			if not options.No2 :
				new_Seq_db.update(read_fasta(options.out + ".2" + ".fasta"))
			new_Seq_db.update(read_fasta(options.out + ".Un" + ".fasta"))

			CDS_file = get_sequence( new_gff3 , new_Seq_db , haplodup_dir + "/new" , "CDS")

			index_dir = haplodup_dir + "/gmap_index"
			mkdir(index_dir)
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Mapping on Hap1"
			print >> sys.stderr , "## Mapping on Hap1"
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Indexing Hap1"
			print >> sys.stderr , "### Indexing Hap1"
			# Hap1
			## index
			indexing_out_file = open( haplodup_dir + "/gmap_index.1.log" ,"w" )
			indexing_err_file = open( haplodup_dir + "/gmap_index.1.err" ,"w" )
			hap1_name = options.out + ".1" + ".fasta"
			indexing_command = "gmap_build -D " + index_dir + " -d hap1.fasta " + hap1_name
			indexProcess = subprocess.Popen(indexing_command, shell=True, stdout=indexing_out_file , stderr=indexing_err_file)
			output, error = indexProcess.communicate()
			indexing_out_file.close()
			indexing_err_file.close()

			# Gmap CDS on results
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Mapping with Gmap "
			gmap_results_1 = haplodup_dir + "/CDS.on.hap1.gmap.gff3"
			gmap_1_gff3 = open( gmap_results_1 , "w" )
			gmap_err = open( gmap_results_1 + ".err" , "w" )
			gmapCommand = "gmap -D " + index_dir + " -d hap1.fasta -f 2 -n 500 -t " + str(options.cores) + " " + CDS_file

			gmapProcess = subprocess.Popen(gmapCommand, shell=True, stdout=gmap_1_gff3 , stderr=gmap_err)
			output, error = gmapProcess.communicate()
			gmap_1_gff3.close()
			gmap_err.close()

			# Hap2
			## Index
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Mapping on Hap2"
			print >> sys.stderr , "## Mapping on Hap2"
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Indexing Hap2"
			print >> sys.stderr , "### Indexing Hap2"
			indexing_out_file = open( haplodup_dir + "/gmap_index.2.log" ,"w" )
			indexing_err_file = open( haplodup_dir + "/gmap_index.2.err" ,"w" )
			hap2_name = options.out + ".2" + ".fasta"
			indexing_command = "gmap_build -D " + index_dir + " -d hap2.fasta " + hap2_name
			indexProcess = subprocess.Popen(indexing_command, shell=True, stdout=indexing_out_file , stderr=indexing_err_file)
			output, error = indexProcess.communicate()
			indexing_out_file.close()
			indexing_err_file.close()
			# Gmap CDS on results
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Mapping with Gmap "
			gmap_results_2 = haplodup_dir + "/CDS.on.hap2.gmap.gff3"
			gmap_2_gff3 = open( gmap_results_2 , "w" )
			gmap_err = open( gmap_results_2 + ".err" , "w" )
			gmapCommand = "gmap -D " + index_dir + " -d hap2.fasta -f 2 -n 500 -t " + str(options.cores) + " " + CDS_file

			gmapProcess = subprocess.Popen(gmapCommand, shell=True, stdout=gmap_2_gff3 , stderr=gmap_err)
			output, error = gmapProcess.communicate()
			gmap_2_gff3.close()
			gmap_err.close()



			# Concatenate Hap1 and Hap2 results
			gmap_results = haplodup_dir + "/CDS.on.genome.gmap.gff3"
			gmap_gff3 = open( gmap_results , "w" )
			for line in open( gmap_results_1 ) :
				print >> gmap_gff3 , line.rstrip()
			for line in open( gmap_results_2 ) :
				print >> gmap_gff3 , line.rstrip()
			gmap_gff3.close()






			# Extract valid alignments per locus
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Extracting valid alignments (identity > 95% , coverage> 95%)"
			gmap_hits_hap1 , gmap_hits_hap2 = read_gmap_results_Hap(gmap_results, 95, 95 , "by_locus", mRNA_to_gene_db)
			test_1 = open( haplodup_dir + "/gmap_hits_hap1.txt" ,'w')
			for chr in sorted(gmap_hits_hap1.keys()) :
				for locus in sorted(gmap_hits_hap1[chr].keys()) :
					for hit in sorted(gmap_hits_hap1[chr][locus]) :
						print >> test_1 , chr + "\t" + "\t".join([ str(x) for x in hit ]) + "\t" + locus
			test_1.close()
			## Make Hap1 gene table
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Extracting loci positions on Hap1"
			hap1_genes = gff3_filter2table_Hap(new_gff3, "gene", "_Hap1_")
			#print >> sys.stderr, "### hap1 gene"
			#print >> sys.stderr, hap1_genes

			test_2 = open( haplodup_dir + "/gmap_hits_hap2.txt" ,'w')
			for chr in sorted(gmap_hits_hap2.keys()) :
				for locus in sorted(gmap_hits_hap2[chr].keys()) :
					for hit in sorted(gmap_hits_hap2[chr][locus]) :
						print >> test_2 , chr + "\t" + "\t".join([ str(x) for x in hit ]) + "\t" + locus
			test_2.close()
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Extracting loci positions on Hap2"
			hap2_genes = gff3_filter2table_Hap(new_gff3, "gene", "_Hap2_")
			#print >> sys.stderr, "### hap2 gene"
			#print >> sys.stderr, hap2_genes

			# Join results
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Counting intra-chromosome hits"
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Hap1"
			hit_counts_1 = do_count_hits_Hap(hap1_genes, gmap_hits_hap1, gmap_hits_hap2, "_Hap1_" , {} )
			hit_file_1 = print_hit_counts(hit_counts_1, haplodup_dir + "/diploid_gene_count_trace.hap1.txt")

			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Hap2"
			hit_counts_2 = do_count_hits_Hap(hap2_genes, gmap_hits_hap1, gmap_hits_hap2, "_Hap2_" , {} )
			hit_file_2 = print_hit_counts(hit_counts_2, haplodup_dir + "/diploid_gene_count_trace.hap2.txt")

			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Generating reports"

			for comparison in coord_tables.keys() :
				plot_files[comparison]["Reports"] = {}
				if isinstance(coord_tables[comparison], list):
					if len(coord_tables[comparison]) == 2 :
						coords_file , coords_file_self = coord_tables[comparison]
					else :
						print >> sys.stdout , "[ERROR] QC comparison with unexpected data content"
						print >> sys.stderr , "[ERROR] QC comparison with unexpected data content: " + comparison + " >>> " + coord_tables[comparison]
						sys.exit(1)
				else :
					if not coord_tables[comparison] == ""  :
						coords_file = coord_tables[comparison]
						coords_file_self = ""
					else :
						print >> sys.stdout , "[ERROR] QC comparison with unexpected data content"
						print >> sys.stderr , "[ERROR] QC comparison " + comparison + " has no data content associated"
						sys.exit(1)
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == " + comparison
				outdir_name = haplodup_dir + "/" + comparison
				# make_pair_html_report usage: make_pair_html_report( coords_file , workdir , output_dir , queryID  , refID  , hap1ID  , hap2ID , "diploid_gene_count_trace.hap1.txt" , "diploid_gene_count_trace.hap2.txt" , min_align = "3000" , similarity = "90" , ratio="0.33")
				if comparison == "Hap1_vs_Reference" :
					for queryID in sorted(hap1_ids.split(",")) :
						refID = hap1_to_ref[queryID]
						hap1ID = queryID
						hap2ID = hap1_to_hap2[hap1ID]
						hap1Len = fasta_1_len[hap1ID]
						hap2Len = fasta_2_len[hap2ID]
						plot_files[comparison]["Reports"][queryID] = {}
						plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, "" , "" , "" ,  "" , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", "3000", "90", "0.33")
						#																			coords, 		coords_self, 			workdir, 				output_dir, 				                                     queryID, refID, structure , legacy , markers, dup_markers , hap1ID , hap2ID , hap1Len , hap2Len , counts_hap1 ="diploid_gene_count_trace.hap1.txt", counts_hap2 ="diploid_gene_count_trace.hap2.txt", min_align ="3000", similarity ="90", ratio="0.33") :
						plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, "3000", "90")
				elif comparison == "Reference_vs_Hap1" :
					for queryID in sorted(ref_ids) :
						if queryID in ref_to_hap1 :
							refID = ref_to_hap1[queryID]
							hap1ID = refID
							hap2ID = hap1_to_hap2[hap1ID]
							hap1Len = fasta_1_len[hap1ID]
							hap2Len = fasta_2_len[hap2ID]
							plot_files[comparison]["Reports"][queryID] = {}
							plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", "3000", "90", "0.33")
							plot_files[comparison]["Reports"][queryID]["pdf"] = make_pair_pdf_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file   , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , "diploid_gene_count_trace.hap1.txt", "3000", "90", "0.33")
							#													make_pair_pdf_report(coords						  , coords_self                       , workdir						  , output_dir				     , queryID, refID, structure        , legacy                , markers               , dup_markers           , hap1ID , hap2ID , counts_hap1                        , min_align , similarity , ratio ) :																																												                                    refID, structure = "" , legacy = "" , markers = "" , dup_markers = "" , hap1ID
						else :
							print >> sys.stderr , "[WARNING] Reference sequence " + queryID + " has no related sequence in Hap1 "
				elif comparison == "Hap1_vs_Hap1" :
					for queryID in sorted(hap1_ids.split(",")) :
						refID = queryID
						hap1ID = queryID
						hap2ID = hap1_to_hap2[hap1ID]
						hap1Len = fasta_1_len[hap1ID]
						hap2Len = fasta_2_len[hap2ID]
						plot_files[comparison]["Reports"][queryID] = {}
						plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name ,  dup_markers_file_name , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", "3000", "90", "0.33")
						#																			coords, 					coords_self, 						workdir, 							output_dir, 				queryID, refID, structure 	   , legacy                , markers,                dup_markers           ,  hap1ID , hap2ID , hap1Len , hap2Len  , counts_hap1                         , counts_hap2                      , min_align , similarity , ratio) :
						plot_files[comparison]["Reports"][queryID]["pdf"] = make_pair_pdf_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , "diploid_gene_count_trace.hap1.txt", "3000", "90", "0.33")
						#													make_pair_pdf_report(coords						  , coords_self                       , workdir						  , output_dir				     , queryID, refID, structure      , legacy                , markers               , dup_markers           , hap1ID , hap2ID , counts_hap1                        , min_align , similarity , ratio ) :
				elif comparison == "Hap2_vs_Reference" :
					for queryID in sorted(hap2_ids.split(",")) :
						refID = hap2_to_ref[queryID]
						hap1ID = hap2_to_hap1[queryID]
						hap2ID = queryID
						hap1Len = fasta_1_len[hap1ID]
						hap2Len = fasta_2_len[hap2ID]
						plot_files[comparison]["Reports"][queryID] = {}
						plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, "" , "" , "" ,  "" , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", "3000", "90", "0.33")
						plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, "3000", "90")
				elif comparison == "Reference_vs_Hap2" :
					for queryID in sorted(ref_ids) :
						if queryID in ref_to_hap2 :
							refID = ref_to_hap2[queryID]
							hap1ID = hap2_to_hap1[refID]
							hap2ID = refID
							hap1Len = fasta_1_len[hap1ID]
							hap2Len = fasta_2_len[hap2ID]
							plot_files[comparison]["Reports"][queryID] = {}
							plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name ,  dup_markers_file_name , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", "3000", "90", "0.33")
							plot_files[comparison]["Reports"][queryID]["pdf"] = make_pair_pdf_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , "diploid_gene_count_trace.hap2.txt", "3000", "90", "0.33")
							#													make_pair_pdf_report(coords						  , coords_self                       , workdir						  , output_dir				     , queryID, refID, structure      , legacy                , markers               , dup_markers           , hap1ID , hap2ID , counts_hap1                        , min_align , similarity , ratio ) :
						else :
							print >> sys.stderr , "[WARNING] Reference sequence " + queryID + " has no related sequence in Hap2 "
				elif comparison == "Hap2_vs_Hap1" :
					for hap2ID in sorted(hap2_ids.split(",")) :
						hap1ID = hap2_to_hap1[hap2ID]
						refID = hap1ID
						queryID = hap2ID
						hap1Len = fasta_1_len[hap1ID]
						hap2Len = fasta_2_len[hap2ID]
						plot_files[comparison]["Reports"][queryID] = {}
						plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name ,  dup_markers_file_name , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", "3000", "90", "0.33")
						plot_files[comparison]["Reports"][queryID]["pdf"] = make_pair_pdf_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , "diploid_gene_count_trace.hap1.txt", "3000", "90", "0.33")
						#													make_pair_pdf_report(coords						  , coords_self                       , workdir						  , output_dir				     , queryID, refID, structure      , legacy                , markers               , dup_markers           , hap1ID , hap2ID , counts_hap1                        , min_align , similarity , ratio ) :
				elif comparison == "Hap1_vs_Hap2" :
					for hap1ID in sorted(hap1_ids.split(",")) :
						if hap1ID in hap1_to_hap2:
							hap2ID = hap1_to_hap2[hap1ID]
							queryID = hap1ID
							refID = hap2ID
							hap1Len = fasta_1_len[hap1ID]
							hap2Len = fasta_2_len[hap2ID]
							plot_files[comparison]["Reports"][queryID] = {}
							plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name ,  dup_markers_file_name ,  hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", "3000", "90", "0.33")
							plot_files[comparison]["Reports"][queryID]["pdf"] = make_pair_pdf_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , "diploid_gene_count_trace.hap2.txt", "3000", "90", "0.33")
							#													make_pair_pdf_report(coords						  , coords_self                       , workdir						  , output_dir				     , queryID, refID, structure      , legacy                , markers               , dup_markers           , hap1ID , hap2ID , counts_hap1                        , min_align , similarity , ratio ) :
						else :
							print >> sys.stderr , "[WARNING] Hap1 sequence " + hap1ID + " has no related sequence in Hap2 "
				elif comparison == "Hap2_vs_Hap2" :
					for hap2ID in sorted(hap2_ids.split(",")) :
						hap1ID = hap2_to_hap1[hap2ID]
						refID = hap2ID
						queryID = hap2ID
						hap1Len = fasta_1_len[hap1ID]
						hap2Len = fasta_2_len[hap2ID]
						plot_files[comparison]["Reports"][queryID] = {}
						plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name ,  dup_markers_file_name , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", "3000", "90", "0.33")
						plot_files[comparison]["Reports"][queryID]["pdf"] = make_pair_pdf_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , "diploid_gene_count_trace.hap2.txt", "3000", "90", "0.33")
						#													make_pair_pdf_report(coords						  , coords_self                       , workdir						  , output_dir				     , queryID, refID, structure      , legacy                , markers               , dup_markers           , hap1ID , hap2ID , counts_hap1                        , min_align , similarity , ratio ) :
				else :
					print >> sys.stderr, "[ERROR] Report required for unknown comparison: " + comparison
					sys.exit(1)
		else:
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Generating reports"

			for comparison in coord_tables.keys() :
				coords_file , coords_file_self = coord_tables[comparison]
				plot_files[comparison]["Reports"] = {}
				print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == " + comparison
				outdir_name = haplodup_dir + "/" + comparison
				mkdir(outdir_name)
				# make_no_genes_html_report( coords_file , haplodup_dir , outdir_name , queryID  , refID , "3000" , "90" )

				if comparison == "Hap1_vs_Reference" :
					for queryID in sorted(hap1_ids) :
						refID = hap1_to_ref[queryID]
						plot_files[comparison]["Reports"][queryID] = {}
						plot_files[comparison]["Reports"][queryID]["html"] = make_no_genes_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), haplodup_dir, outdir_name, queryID, refID, "3000", "90")
						plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, "3000", "90")
				elif comparison == "Reference_vs_Hap1" :
					for queryID in sorted(ref_ids) :
						if queryID in ref_to_hap1 :
							refID = ref_to_hap1[queryID]
							plot_files[comparison]["Reports"][queryID] = {}
							plot_files[comparison]["Reports"][queryID]["html"] = make_no_genes_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), haplodup_dir, outdir_name, queryID, refID, "3000", "90")
							plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, "3000", "90")
						else :
							print >> sys.stderr , "[WARNING] Reference sequence " + queryID + " has no related sequence in Hap1"
				elif comparison == "Hap1_vs_Hap1" :
					continue
				elif comparison == "Hap2_vs_Reference" :
					for queryID in sorted(hap2_ids) :
						refID = hap2_to_ref[queryID]
						plot_files[comparison]["Reports"][queryID] = {}
						plot_files[comparison]["Reports"][queryID]["html"] = make_no_genes_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), haplodup_dir, outdir_name, queryID, refID, "3000", "90")
						plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, "3000", "90")
				elif comparison == "Reference_vs_Hap2" :
					for queryID in sorted(ref_ids) :
						if queryID in ref_to_hap2 :
							refID = ref_to_hap2[queryID]
							plot_files[comparison]["Reports"][queryID] = {}
							plot_files[comparison]["Reports"][queryID]["html"] = make_no_genes_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), haplodup_dir, outdir_name, queryID, refID, "3000", "90")
							plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, "3000", "90")
						else :
							print >> sys.stderr , "[WARNING] Reference sequence " + queryID + " has no related sequence in Hap2"
				elif comparison == "Hap2_vs_Hap1" :
					for queryID in sorted(hap2_ids) :
						refID = hap2_to_hap1[queryID]
						plot_files[comparison]["Reports"][queryID] = {}
						plot_files[comparison]["Reports"][queryID]["html"] = make_no_genes_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), haplodup_dir, outdir_name, queryID, refID, "3000", "90")
						plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, "3000", "90")
				elif comparison == "Hap1_vs_Hap2" :
					for queryID in sorted(hap1_ids) :
						if queryID in hap1_to_hap2 :
							refID = hap1_to_hap2[queryID]
							plot_files[comparison]["Reports"][queryID] = {}
							plot_files[comparison]["Reports"][queryID]["html"] = make_no_genes_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), haplodup_dir, outdir_name, queryID, refID, "3000", "90")
							plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, "3000", "90")
						else :
							print >> sys.stderr , "[WARNING] Hap1 sequence " + hap1ID + " has no related sequence in Hap2 "
				elif comparison == "Hap2_vs_Hap2" :
					continue
				else :
					print >> sys.stderr, "[ERROR] Report required for unknown comparison: " + comparison
					sys.exit(1)

		# Make Index
		html_index = make_index_from_report_db("index.html" , "." , haplodup_dir ,  plot_files  )

	# Clean up
	if not options.agp :
		# Remove hidden intermediate files
		os.remove("." + options.out + ".1" + ".agp")
		os.remove("." + options.out + ".Un" + ".agp")
		if not options.No2 :
			os.remove("." + options.out + ".2" + ".agp")

	##### Finished

	print >> sys.stdout , "------------------------------"
	print >> sys.stdout,  "- Done"
	print >> sys.stdout , "------------------------------"
	print >> sys.stderr , "##############################"
	print >> sys.stderr , "# Done"
	print >> sys.stderr , "##############################"

if __name__ == '__main__':
	main()
