#!/usr/bin/env python

import argparse
import gzip
import sys
import gc
import os
import datetime
import subprocess
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq

gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)



def make_map_graph( hit_list , max_distance ) :

	hit_graph = nx.DiGraph()
	stop_nodes = []
	start_nodes = []
	#Hit format [ id , Tstart , Tstop , Qstart , Qstop , matches , hitLen ]

	# Make nodes and add mapping edges
	for hit in hit_list :
		id , Tstart , Tstop , Qstart , Qstop , matches , hitLen = hit
		align = int(Tstop)-int(Tstart)

		# make nodes
		hit_graph.add_node( (Tstart , Qstart) )
		hit_graph.add_node( (Tstop , Qstop) )

		stop_nodes.append( (Tstop , Qstop) )
		start_nodes.append( (Tstart , Qstart) )

		# add mapping edge
		hit_graph.add_edge( (Tstart,Qstart) , (Tstop,Qstop) , id=id , align=int(align) , match=int(matches))

	#print >> sys.stderr, hit_graph.edges()

	### Add allowed gap edges between nodes

	stop_nodes_len = len(stop_nodes)
	start_nodes_len = len(start_nodes)

	#print >> sys.stderr, start_nodes
	#print >> sys.stderr, stop_nodes

	for i in range(stop_nodes_len) :
		for j in range(start_nodes_len):
			Added="False"
			f_Tstop , f_Qstop = stop_nodes[i]
			t_Tstart , t_Qstart = start_nodes[j]
			T_gap = int(t_Tstart) - int(f_Tstop)
			Q_gap = int(t_Qstart) - int(f_Qstop)
			if T_gap >= 0 and T_gap <= max_distance and Q_gap >= 0 and Q_gap <= max_distance :
				if not hit_graph.has_edge( (f_Tstop , f_Qstop) , (t_Tstart , t_Qstart) ) :
					if (f_Tstop , f_Qstop) != (t_Tstart , t_Qstart) :
						hit_graph.add_edge( (f_Tstop , f_Qstop) , (t_Tstart , t_Qstart) , id="gap", align=0 , match=0)
						#Added="True"
						#print >> sys.stderr , "(" + str(f_Tstop) + ":" + str(f_Qstop) + ") -> (" + str(t_Tstart) + ":" + str(t_Qstart) + ") | Addedd " + Added + " | Gaps: "  + str(T_gap) + ":" + str(Q_gap)

	return hit_graph


def get_subgraph_from_path_tuples( original_graph , selected_path , Tid , Qid , Tlen , Qlen):
	absQid = Qid[:-2]
	orietantion = Qid[-1]
	matching_regions = []
	subgraph = nx.DiGraph()

	for i in range( 0 , len(selected_path) -1 ):
		from_node = selected_path[i]
		#print >> sys.stderr, "from_node"
		#print >> sys.stderr, from_node

		to_node = selected_path[i+1]
		#print >> sys.stderr, "to_node"
		#print >> sys.stderr, to_node

		match_info = original_graph[from_node][to_node]

		Tstart , Qstart = from_node
		Tstop , Qstop = to_node

		align = match_info["align"]
		match = match_info["match"]
		match_id = match_info["id"]
		if not match_id == "gap" :
			matching_regions.append( [ Tid , Tlen , Tstart , Tstop , absQid , Qlen , Qstart , Qstop , orietantion , align , match ]  )

		subgraph.add_node(from_node)
		subgraph.add_node(to_node)

		subgraph.add_edge(from_node,to_node)
		subgraph[from_node][to_node].update(match_info)

	return subgraph , matching_regions


def read_fasta( filename ) :
	sequences_db = {}
	seq = ""
	try:
		for line in gzip.open( filename ):
			if line[0] == ">" :
				if not seq == "" :
					# close open sequences
					sequences_db[seq_id] = seq
				# open new one
				seq_id = line.rstrip()[1:].split()[0]
				seq = ""
			else :
				# concatenate line to seq
				seq += line.rstrip().upper()
		if not seq == "" :
			# close open sequences
			sequences_db[seq_id] = seq

	except IOError:
		for line in open( filename ) :
			if line[0] == ">" :
				if not seq == "" :
					# close open sequences
					sequences_db[seq_id] = seq
				# open new one
				seq_id = line.rstrip()[1:].split()[0]
				seq = ""
			else :
				# concatenate line to seq
				seq += line.rstrip().upper()
		if not seq == "" :
			# close open sequences
			sequences_db[seq_id] = seq

	return sequences_db


def write_fasta_from_db( fasta_db , out_name , compress = False ) :

	if compress :
		out_fasta_file = gzip.open(out_name, "w")
	else :
		out_fasta_file = open(out_name, "w")

	for seq_id in sorted(fasta_db.keys()) :
		#print >> sys.stderr , "Printing " + seq_id
		print >> out_fasta_file , ">" + str(seq_id)
		seq_len = len(fasta_db[seq_id])
		for i in range( 0 , seq_len , 60 ):
			print >> out_fasta_file , fasta_db[seq_id][i : min( i + 60 , seq_len) ]

	out_fasta_file.close()

	return out_name


def read_nucmer_coords( coords_file ) :
	# coords file as produced with show-coords + " -l -r -T -H "
	original_hits = {}
	id=0
	for line in open(coords_file) :
		id+=1
		if line[0] == "#" or line.rstrip() == "": continue ;
		if line[0] == "" or line.rstrip() == "": continue ;
		Tstart , Tstop , Qstart , Qstop , hitLen , QhitLen , iden , Tlen , Qlen , Tcov , Qcov , Tid , Qid = line.rstrip().split("\t")[0:13]
		matches = int( float(hitLen) * float(iden)/100 )

		if ( Tid , Qid ) not in original_hits :
			# create a new entry
			original_hits[(Tid,Qid)] = []

		# Add a new range to the db
		original_hits[(Tid,Qid)].append([ id , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ])
	return original_hits


def map_nucmer( ref_file , query_file ,  cores ,  out_file_name , nucmer_path , showcoords_path , parameters = "" , filter = "" ) :
	out_file_name_prefix = out_file_name.split(".coords")[0]

	if nucmer_path == "" :
		nucmer_search=subprocess.Popen( "which nucmer" , shell=True, stdout=subprocess.PIPE )
		nucmer_command_line , error = nucmer_search.communicate()
		nucmer_command_line = nucmer_command_line.rstrip()
		if nucmer_command_line == "" :
			print >> sys.stderr , '[ERROR] Nucmer expected to be in $PATH, not found'
			exit(1)
	else :
		nucmer_command_line = nucmer_path + "/nucmer"
	nucmer_command_line += " -p " + out_file_name_prefix + " -t " + str(cores) + " " + parameters + " "
	map_file_err = open( out_file_name_prefix + ".err" , "w")
	print >> sys.stderr, "### Running command line: " + nucmer_command_line + ref_file + " " + query_file + " "
	mapProcess = subprocess.Popen(nucmer_command_line + ref_file + " " + query_file + " ", shell=True, stderr=map_file_err)
	output, error = mapProcess.communicate()
	map_file_err.close()

	if showcoords_path == "" :
		showcoords_search=subprocess.Popen( "which show-coords" , shell=True, stdout=subprocess.PIPE )
		nucmer_command_line , error = showcoords_search.communicate()
		extract_coords_process = nucmer_command_line.rstrip()
		if extract_coords_process == "" :
			print >> sys.stderr , '[ERROR] show-coords expected to be in $PATH, not found'
			exit(1)
	else :
		extract_coords_process = showcoords_path + "/show-coords"

	extract_coords_process += " -c " + filter + " "
	coords_file = open( out_file_name , 'w' )
	input_delta = out_file_name_prefix + ".delta"
	print >> sys.stderr, "### Running command line: " + extract_coords_process + input_delta
	coordsProcess = subprocess.Popen(extract_coords_process + input_delta, shell=True, stdout=coords_file)
	output, error = coordsProcess.communicate()
	coords_file.close()

	return out_file_name


def hit_global( hits_file , max_gap_size , seq_len_db , intermediate_file , ordering_feature ) :
	print >> sys.stderr, "## Merge hits in blocks and select best alignments"
	merged_hits = {}
	unique_hits = {}
	if ordering_feature == "coverage" :
		weighting_feature = 'align'
		sorting_cell = 10
	elif ordering_feature == "identity" :
		weighting_feature = 'match'
		sorting_cell = 9
	else :
		print >> sys.stderr, "[ERROR] Unknown results sorting feature: " + ordering_feature
		exit(4)

	# Write uniquified alignments ian paf-like format
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
	#  10 - Alignment block length
	#  11 - Mapping quality (0-255; 255 for missing)
	#
	########################

	print >> sys.stderr, "### Read mapping hits"
	original_hits = read_nucmer_coords( hits_file )

	# original_hits[(Tid,Qid)] =
	# 	[
	# 		[ Qid , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ] , ...
	#	]
	print >> sys.stderr, "#### " + str(len(original_hits.keys())) + " hits read"
	print >> sys.stderr, "### Merging hits"

	for entry in sorted(original_hits.keys()) :
		new_merged_hit = 13*["-"]
		Tid , Qid = entry

		#print >> sys.stderr, "#### Target: " + Tid + " - Query: " + Qid
		#print >> sys.stderr, "##### Target length: " + str(seq_len_db[Tid])
		#print >> sys.stderr, "##### Query length: " + str(seq_len_db[Qid])

		map_graph = make_map_graph( original_hits[entry], max_gap_size )
		try :
			longest_merged_path = nx.dag_longest_path(map_graph, weight=weighting_feature)
			longest_graph , matching_regions = get_subgraph_from_path_tuples( map_graph , longest_merged_path , Tid , Qid , seq_len_db[Tid] , seq_len_db[Qid] )
			longest_graph_matches = longest_graph.size(weight='match')

		except :
			print >> sys.stderr, "[ERROR] Unable to merge hits"
			print >> sys.stderr, "[ERROR] Unable to merge hits of " + Qid + " on " + Tid + " (" + str(len(original_hits[entry])) + " hits)"
			print >> sys.stderr, sorted(map_graph.edges.data())
			print >> sys.stderr, nx.find_cycle(map_graph)
			exit(5)

		#print >> sys.stderr, "##### longest_merged_path"
		#print >> sys.stderr, str(longest_merged_path)
		#print >> sys.stderr, "##### longest_graph_matches"
		#print >> sys.stderr, longest_graph_matches
		#print >> sys.stderr, "##### matching_regions"
		#print >> sys.stderr, matching_regions

		absQid = Qid[:-2]
		orientation = Qid[-1]

		new_merged_hit[0] = absQid
		new_merged_hit[1] = seq_len_db[Qid]
		new_merged_hit[4] = orientation
		new_merged_hit[5] = Tid
		new_merged_hit[6] = seq_len_db[Tid]

		Tstart , Qstart = longest_merged_path[0]
		Tstop , Qstop = longest_merged_path[-1]

		new_merged_hit[2] = int(Qstart)
		new_merged_hit[3] = int(Qstop)
		new_merged_hit[7] = int(Tstart)
		new_merged_hit[8] = int(Tstop)
		new_merged_hit[9] = int(longest_graph_matches)
		new_merged_hit[10] = int(new_merged_hit[8]) - int(new_merged_hit[7]) + 1
		new_merged_hit[11] = int( float(255 * (float(1) - (float(new_merged_hit[9]) / float(new_merged_hit[1]) ) ) ) )

		#print >> merged_temp_file , "\t".join(str(x) for x in new_merged_hit)

		# Add the new merged hit to the database
		absQid = Qid[:-2] # absolute id of the sequence, without |+ or |- at the end
		if absQid not in merged_hits :
			merged_hits[absQid] = []
		new_merged_hit[12] = matching_regions
		merged_hits[absQid].append(new_merged_hit)
		#print >> sys.stderr, "#### merged_hits[" + absQid + "] (number of hits: " + str(len(merged_hits[absQid])) + ")"
		#print >> sys.stderr, merged_hits[absQid]

	for Qid in sorted(merged_hits.keys()) : # For each absolute id find the best alignment
		merged_list = sorted(merged_hits[Qid] , key=lambda item: int(item[sorting_cell]) , reverse = True )
		for hit in merged_list :
			hit_line = "\t".join([ str(x) for x in hit ] )
			print >> intermediate_file, hit_line

		best_hit = merged_list[0]
		#print >> sys.stderr, "merged_list:"
		#print >> sys.stderr, merged_list
		unique_hits[Qid] = {}
		unique_hits[Qid]["region"] = best_hit[0:12]
		unique_hits[Qid]["hits"] = best_hit[12]

	return unique_hits



def main() :

	###### Options and help ######
	parser = argparse.ArgumentParser()

	parser.add_argument("-t", "--target", dest="target",
					help="FASTA file of the target sequences. [REQUIRED]", metavar="target.fasta")
	parser.add_argument("-q", "--query", dest="query",
					help="FASTA file of the query sequences. [REQUIRED]", metavar="query.fasta")
	parser.add_argument("-o", "--output", dest="out", default="out",
					help="Output files prefix [Default: out]", metavar="out_prefix")
	parser.add_argument("-c", "--cores", dest="cores", default=8,
					help="Cores used in mapping process [default: 8]", metavar="N")
	parser.add_argument("--hitgap" , dest="hitgap", default="500000",
					help="Maximum gap size between hits allowed in path extension [default: 500000]", metavar="N")
	parser.add_argument("--tool_path" , dest="tool_path" , default="" ,
					help="Path to nucmer installation folder [default: system installation]", metavar="/path/to/mummer_folder")
	parser.add_argument("--maximize" , dest="order", default="coverage" ,
						help="Select the best alignment by maximizing the coverage or the identity. Accepted values 'coverage' or 'identity'. [Default: coverage]")
	parser.add_argument("--reuse_mappings" , dest="reuse_mappings", default=False, action="store_true",
						help="If set, alignments present in the output folder are reused and not overwritten by performing alignments again [Default: overwrite]")

	# Sanity Check

	if len(sys.argv) < 3:
		parser.print_help()
		sys.exit(1)

	options = parser.parse_args()

	ordering = str(options.order).lower()

	if not ( options.target and options.query ) :
		print >> sys.stderr , "[ERROR] Missing FASTA sequences files(s)"
		parser.print_help()
		sys.exit(1)

	tool_path = options.tool_path
	nucmer_path = tool_path
	if tool_path == "" :
		nucmer_search=subprocess.Popen( "which nucmer" , shell=True, stdout=subprocess.PIPE )
		nucmer_command_line , error = nucmer_search.communicate()
		nucmer_command_line = nucmer_command_line.rstrip()
		if nucmer_command_line == "" :
			print >> sys.stderr , '[ERROR] Nucmer expected to be in $PATH, not found'
			exit(1)
	else :
		nucmer_search=subprocess.Popen( "which " + nucmer_path + "/nucmer" , shell=True, stdout=subprocess.PIPE )
		nucmer_command_line , error = nucmer_search.communicate()
		nucmer_command_line = nucmer_command_line.rstrip()
		if nucmer_command_line == "" :
			print >> sys.stderr , '[ERROR] Nucmer expected to be in '+ nucmer_path + ', not found'
			exit(1)

	showcoords_path = tool_path

	if showcoords_path == "" :
		showcoords_search=subprocess.Popen( "which show-coords" , shell=True, stdout=subprocess.PIPE )
		nucmer_command_line , error = showcoords_search.communicate()
		extract_coords_process = nucmer_command_line.rstrip()
		if extract_coords_process == "" :
			print >> sys.stderr , '[ERROR] show-coords expected to be in $PATH, not found'
			exit(1)
	else :
		showcoords_search=subprocess.Popen( "which " + showcoords_path + "/show-coords" , shell=True, stdout=subprocess.PIPE )
		nucmer_command_line , error = showcoords_search.communicate()
		extract_coords_process = nucmer_command_line.rstrip()
		if extract_coords_process == "" :
			print >> sys.stderr , '[ERROR] show-coords expected to be in ' + showcoords_path + ', not found'
			exit(1)


	# Test presence of nucmer

	tmp_dir = options.out + ".temp_dir"
	if not os.path.isdir( tmp_dir ) :
		os.mkdir( tmp_dir )

	fasta_len_dict = {}

	# Make sequences
	## Target
	# options.target
	target_fasta = tmp_dir + "/target.fasta"
	target_temp_fasta_db = read_fasta(options.target)
	target_fasta_db = {}
	for seq_id in sorted(target_temp_fasta_db.keys()) :
		target_fasta_db[seq_id] = target_temp_fasta_db[seq_id].upper()
		fasta_len_dict[seq_id] = len(target_temp_fasta_db[seq_id])
	write_fasta_from_db( target_fasta_db , target_fasta )

	## Read query -> + and - sequences
	query_fasta = tmp_dir + "/query.fasta"
	query_fasta_db = {}
	query_temp_fasta_db =read_fasta(options.query)

	for seq_id in sorted(query_temp_fasta_db.keys()) :
		query_fasta_db[seq_id + "|+"] = query_temp_fasta_db[seq_id].upper()
		query_fasta_db[seq_id + "|-"] = str(Seq(query_temp_fasta_db[seq_id]).reverse_complement()).upper()
		fasta_len_dict[seq_id + "|+"] = len(query_temp_fasta_db[seq_id])
		fasta_len_dict[seq_id + "|-"] = len(query_temp_fasta_db[seq_id])

	write_fasta_from_db( query_fasta_db , query_fasta )

	# Map
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Mapping query sequences"
	print >> sys.stderr, "# Mapping query sequences"
	coords_file = tmp_dir + "/query.on.target.coords"
	if options.reuse_mappings :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Reading pre-calculated results"
		print >> sys.stderr, "## Reading pre-calculated results"
	else :
		coords_file = map_nucmer( target_fasta , query_fasta ,  int(options.cores) ,  coords_file , nucmer_path , showcoords_path , " --forward " , " -l -r -T -H ")

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Sorting hits on " + ordering
	print >> sys.stderr, "## Sorting hits on " + ordering
	# Extract best alignment region
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Extractign best global alignments" 
	print >> sys.stderr, "# Extractign best global alignments"
	temp_results_file = open(options.out + ".all_tiles.tsv" , 'w')
	global_results = hit_global(coords_file , int(options.hitgap) , fasta_len_dict , temp_results_file , ordering )
	temp_results_file.close()
	regions_file = open(options.out + ".regions.tsv" , 'w')
	hits_file = open(options.out + ".hits.tsv" , 'w')
	for id in sorted(global_results.keys()) :
		region_paf = global_results[id]["region"][:]
		global_hits = global_results[id]["hits"][:]
		# [ ... , [ Tid , Tlen , Tstart , Tstop , absQid , Qlen , Qstart , Qstop , orietantion , align , match ] , ... ]
		print >> regions_file , "\t".join([str(x) for x in region_paf])
		for single_hit in global_hits :
			print >> hits_file , "\t".join([str(x) for x in single_hit])

	regions_file.close()
	hits_file.close()

if __name__ == '__main__':
	main()