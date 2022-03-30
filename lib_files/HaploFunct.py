#!/usr/bin/env python

import subprocess
import sys
import operator
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
import gc
import os
import pipes
import re
import numpy as np
import scipy
from scipy.signal import butter, filtfilt
import gzip
import json
import shutil
import toml
import yaml
import pickle
import datetime
import itertools
from GFF_lib import *
from AGP_lib import *
from FASTA_lib import *
from map_lib import *
from collections import Counter

gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
scriptDirectory = os.path.dirname(os.path.realpath(__file__)) + "/../support_scripts"


###### Functions ######

def get_version() :
	version_file = os.path.dirname(os.path.realpath(__file__))  + "/version.txt"
	for line in open( version_file ) :
		version = line.rstrip().split("\t")[1]
		break
	return version


def isDigit(x):
	try:
		float(x)
		return True
	except ValueError:
		return False


def touch(filename) :
	file=open(filename, "w")
	file.close()
	return filename


def JSON_key_to_int( dict_from_json ):
	correctedDict = {}

	for dict_key, item_value in dict_from_json.items():
		if isinstance(item_value, list) :
			item_value = [JSON_key_to_int(item) if isinstance(item, dict) else item for item in item_value]
		elif isinstance( item_value, dict ) :
			item_value = JSON_key_to_int(item_value)

		try:
			dict_key = int(dict_key)
		except :
			dict_key = dict_key
		correctedDict[dict_key] = item_value

	return correctedDict


def get_filename( x ):
	base=os.path.basename('/root/dir/sub/file.ext')
	return str(os.path.splitext(base)[0])


def remove_extension( file , ext , keep=True):
	if keep :
		return str(re.sub(str(ext) + "$", "", file))
	else :
		base=os.path.basename('/root/dir/sub/file.ext')
		return str(re.sub(str(ext) + "$", "", base))


def set_paths( conf_file ) :
	path_db=toml.load(conf_file)["paths"]
	return path_db


def get_subgraph_from_path( original_graph , selected_path ):

	edges = []
	edges_names=[]

	for i in range( 0 , len(selected_path) -1 ):
		from_node = selected_path[i]
		to_node = selected_path[i+1]
		element = original_graph[from_node][to_node]
		if "name" in element :
			edges_names.append(element["name"])
			edges.append([ element["name"] , "(" + str(from_node) + ":" + str(to_node) + ")" ,  element["region"] , element["align"] , element["match"] ])
		else :
			try :
				edges.append([ element["length"] , "(" + str(from_node) + ":" + str(to_node) + ")" ,  element["region"] , element["align"] , element["match"] ])
			except :
				print >> sys.stderr , element
				sys.exit(1)

	return edges , edges_names


def get_subgraph_from_path_tuples( original_graph , selected_path ):

	subgraph = nx.DiGraph()

	for i in range( 0 , len(selected_path) -1 ):
		from_node = selected_path[i]
		to_node = selected_path[i+1]

		subgraph.add_node(from_node)
		subgraph.add_node(to_node)

		subgraph.add_edge(from_node,to_node)
		subgraph[from_node][to_node].update(original_graph[from_node][to_node])

	return subgraph


def get_node_matchLength( original_graph, selected_path) :
	align_length_list = []
	for i in range( 0 , len(selected_path) ):
		#print >> sys.stderr, original_graph.node[selected_path[i]]['align_length']
		align_length_list.append(int(original_graph.node[selected_path[i]]['align_length']))
	return( align_length_list )


def make_graph( hit_list , max_distance , blacklist) :

	hit_graph = nx.DiGraph()
	stop_nodes = []
	start_nodes = []
	#Hit format [ Qid , Tstart , Tstop , Qstart , Qstop , matches , hitLen ]

	if blacklist != [] :
		print >> sys.stderr , "### Blacklist: " + " ".join(blacklist)

	# Make nodes and add mapping edges
	for hit in hit_list :
		Qid , Tstart , Tstop , Qstart , Qstop , matches , hitLen = hit
		align = int(Tstop)-int(Tstart)
		Q_region=str(Qstart)+":"+str(Qstop)

		if Qid not in blacklist :
			# make nodes
			hit_graph.add_node(int(Tstart))
			hit_graph.add_node(int(Tstop))

			# add mapping edge
			if not ( Tstart == Tstop ) :
				# Exclude ChrStart an ChrStop self
				hit_graph.add_edge( int(Tstart) , int(Tstop) , align=int(align) , match=int(matches), name=Qid, region=Q_region)
				stop_nodes.append(int(Tstop))
				start_nodes.append(int(Tstart))
			else :
				if Tstart == 0 :
					stop_nodes.append(int(Tstop))
				else :
					start_nodes.append(int(Tstart))

	### Add allowed gap edges between nodes

	stop_nodes_len = len(stop_nodes)
	start_nodes_len = len(start_nodes)
	stop_nodes.sort()
	start_nodes.sort()

	print >> sys.stderr, sorted(hit_graph.edges())
	print >> sys.stderr, sorted(start_nodes)
	print >> sys.stderr, sorted(stop_nodes)

	for i in range(stop_nodes_len) :
		for j in range(start_nodes_len):
			T_from = stop_nodes[i]
			T_to = start_nodes[j]
			T_gap = int(T_to) - int(T_from)

			#if (i == 0) or (j == start_nodes_len-1) :
			#	if not hit_graph.has_edge( T_from , T_to ) and not ( (i == 0) and (j == start_nodes_len-1) ) :
			#		hit_graph.add_edge( T_from , T_to , align=0 , match=0, name=T_gap, region="gap")
			#		print >> sys.stderr , "(" + str(T_from) + ") -> (" + str(T_to) + ") | Gap: "  + str(T_gap)

			if T_gap > 0  and  T_gap <= max_distance and not hit_graph.has_edge( T_from , T_to ):
				hit_graph.add_edge( T_from , T_to , align=0 , match=0, length=T_gap, region="gap")
				print >> sys.stderr , "(" + str(T_from) + ") -> (" + str(T_to) + ") | Gap: "  + str(T_gap)

	return( hit_graph )


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


def make_forced_graph( hit_list , max_distance , forced_sorted_list, blacklist, err_filename) :

	hit_graph = nx.DiGraph()
	stop_nodes = []
	start_nodes = []
	forced_nodes = {}
	forced_nodes_ranges = []

	if blacklist != [] :
		print >> sys.stderr , "### Blacklist: " + " ".join(blacklist)

	#Hit format [ Qid , Tstart , Tstop , Qstart , Qstop , matches , hitLen ]

	# Make nodes and add mapping edges from forced list
	for hit in hit_list :
		Qid , Tstart , Tstop , Qstart , Qstop , matches , hitLen = hit

		if Qid in forced_sorted_list :
			if Qid not in blacklist :
				print >> sys.stderr , "###### Adding " + Qid + " to graph : Forced"
				align = int(Tstop)-int(Tstart)
				Q_region=str(Qstart)+":"+str(Qstop)

				forced_nodes_ranges.append([int(Tstart),int(Tstop)])
				forced_nodes[Qid] = [int(Tstart),int(Tstop)]

				# make nodes
				hit_graph.add_node(int(Tstart))
				hit_graph.add_node(int(Tstop))

				# add mapping edge
				if not ( Tstart == Tstop ) :
					# Exclude ChrStart an ChrStop self
					hit_graph.add_edge( int(Tstart) , int(Tstop) , align=int(align) , match=int(matches), name=Qid, region=Q_region)
					stop_nodes.append(int(Tstop))
					start_nodes.append(int(Tstart))
				else :
					if Tstart == 0 :
						stop_nodes.append(int(Tstop))
					else :
						start_nodes.append(int(Tstart))

					#edge = hit_graph[int(Tstart)][int(Tstop)]
					#print >> sys.stderr , "####### edge: " + str(Tstart) + " " + str(Tstop) + " - " + ", ".join( str(key) + ": " + str(edge[key]) for key in sorted(edge.keys()) )

		#	else:
		#		print >> sys.stdout , "----- Error: " + Qid + " present both as required and blacklisted."
		#		print >> sys.stderr , "##### Error: " + Qid + " present both as required and blacklisted."
		#		sys.exit(3)

	# Forced Qids mappings sanity check
	print >> sys.stdout , "---- Testing required query sequences mapping positions order"
	prev_id , prev_start , prev_stop = [ "" , "" ,"" ]

	for id in forced_sorted_list :

		if id not in forced_nodes :
			print >> sys.stdout , "----- " + str(id) + " not mapped, removing from required"
			print >> sys.stderr , "##### " + str(id) + " not mapped, removing from required"
			#print >> sys.stderr , "##### forced_sorted_list: " + str(forced_sorted_list)
			continue

		start = int(forced_nodes[id][0])
		stop = int(forced_nodes[id][1])

		if not prev_id == "" :
			# Not the first element, test regions
			if (stop <= prev_stop) or (start <= prev_start) :
				# Mapping order not compatible
				print >> sys.stdout , "[ERROR] Incompatibility issue: (" + prev_id + " -> " + id + ") order not coherent with mapping results."
				print >> sys.stderr , "[ERROR] Incompatibility issue: (" + prev_id + " -> " + id + ") order not coherent with mapping results. Check " +  err_filename + " for further information"
				# Generate error tracking file
				errorfile = open(err_filename , "w")
				print >> errorfile , "ID\tAlignment_Start\tAlignment_Stop"
				for id in forced_sorted_list :
					if id in forced_nodes :
						print >> errorfile , id + "\t" + str(forced_nodes[id][0]) + "\t" + str(forced_nodes[id][1])
				# Exit with error
				sys.exit(2)
			else :
				if ( start < prev_stop ) :
					# mapping regions do overlap. Update node and edge in the network to allow the path, keep alignment length and matches unmodified
					#print >> sys.stdout , "----- (" + prev_id + " -> " + id + ") overlap in mapping results."
					new_start = prev_stop
					hit_graph.add_edge( new_start , stop)
					hit_graph[new_start][stop].update(hit_graph[start][stop])
					hit_graph.remove_node(start)

		# Update
		prev_id = id
		prev_start = start
		prev_stop = stop

	#print >> sys.stderr , "##### Forced query sequences mappings"
	#print >> sys.stderr , "##### ID\tAlignment_Start\tAlignment_Stop"
	#for id in forced_sorted_list :
	#	if id in forced_nodes :
	#		print >> sys.stderr ,"##### " + id + "\t" + str(forced_nodes[id][0]) + "\t" + str(forced_nodes[id][1])

	# Make nodes and add mapping edges of Qid compatible with forced list
	for hit in hit_list :
		Qid , Tstart , Tstop , Qstart , Qstop , matches , hitLen = hit
		#print >> sys.stderr , "#### " + Qid + " Required: " +  str(Qid in forced_sorted_list) + " ; blacklisted: " + str(Qid in blacklist)
		start_nodes.sort()
		stop_nodes.sort()

		if Qid not in forced_sorted_list:
			if Qid not in blacklist :
				# Compatibility test
				add = True
				for query_range in forced_nodes_ranges :
					if ( query_range[0] <= Tstart <= query_range[1] ) or (query_range[0] <= Tstop <= query_range[1]) :
						add = False
				if add :
					print >> sys.stderr , "###### Adding " + Qid + " to graph: Usable"
					align = int(Tstop)-int(Tstart)
					Q_region=str(Qstart)+":"+str(Qstop)

					# make nodes
					hit_graph.add_node(Tstart)
					hit_graph.add_node(Tstop)

					# add mapping edge
					if not ( Tstart == Tstop ) :
						# Exclude ChrStart an ChrStop self
						hit_graph.add_edge( int(Tstart) , int(Tstop) , align=int(align) , match=int(matches), name=Qid, region=Q_region)
						stop_nodes.append(int(Tstop))
						start_nodes.append(int(Tstart))
						#print >> sys.stderr , "####### edge: " + str(Tstart) + " " + str(Tstop) + " - " + ", ".join( str(key) + ": " + str(edge[key]) for key in sorted(edge.keys()) )

					else :
						#print >> sys.stderr , "###### Chromosome extremity"
						if int(Tstart) == 0 :
							#print >> sys.stderr , "####### Chromosome start"
							#print >> sys.stderr, sorted(stop_nodes)
							stop_nodes.append(int(Tstop))
							stop_nodes.sort()
							#print >> sys.stderr, sorted(stop_nodes)
						else :
							#print >> sys.stderr , "####### Chromosome end"
							#print >> sys.stderr, sorted(start_nodes, reverse=True)
							start_nodes.append(int(Tstart))
							start_nodes.sort()
							#print >> sys.stderr, sorted(start_nodes, reverse=True)
				#else:
				#	print >> sys.stderr , "###### Refused " + Qid + " to graph, incompatible"

	### Add allowed gap edges between nodes

	stop_nodes_len = int(len(stop_nodes))
	start_nodes_len = int(len(start_nodes))

	#print >> sys.stderr, "##### Graph mapping contigs edges"
	#print >> sys.stderr, sorted(hit_graph.edges())
	#print >> sys.stderr, sorted(start_nodes)
	#print >> sys.stderr, start_nodes_len
	#print >> sys.stderr, sorted(stop_nodes)
	#print >> sys.stderr, stop_nodes_len

	print >> sys.stderr, "##### Add allowed gaps"
	for i in range(stop_nodes_len) :
		for j in range(start_nodes_len):
			T_from = stop_nodes[i]
			T_to = start_nodes[j]
			T_gap = int(T_to) - int(T_from)

			#if (i == 0) or (j == start_nodes_len-1) :
			#	if not hit_graph.has_edge( T_from , T_to ) and not ( (i == 0) and (j == start_nodes_len-1) ) :
			#		hit_graph.add_edge( T_from , T_to , align=0 , match=0, name=T_gap, region="gap")
			#		print >> sys.stderr , "(" + str(T_from) + ") -> (" + str(T_to) + ") | Gap: "  + str(T_gap)

			if T_gap > 0  and  T_gap <= max_distance :
				if not hit_graph.has_edge( T_from , T_to ):
					hit_graph.add_edge( T_from , T_to , align=0 , match=0, length=T_gap, region="gap")
					#print >> sys.stderr , "(" + str(T_from) + ") -> (" + str(T_to) + ") | Gap: "  + str(T_gap)

	return( hit_graph )


def make_tiling_paths( new_graph ) :

	print >> sys.stdout , "- Generating all tiling paths"
	graph_order = new_graph.order()
	print >> sys.stdout , "-- Graph order: " + str(graph_order)
	print >> sys.stdout , "-- Graph size: " + str(new_graph.size())
	print >> sys.stdout , "-- Graph traversal and tiling path identification"
	tilings = nx.all_simple_paths(new_graph, source="Chr_start" , target="Chr_end" , cutoff=round(graph_order/2) )
	paths = []

	print >> sys.stdout , "-- Calculating length, identity and coverage for the tiling paths"
	counter = 0
	for t_path in tilings :
		#print >> sys.stdout, counter
		#print >> sys.stderr, Tpath
		t_path_graph = get_subgraph_from_path(new_graph, t_path)
		t_path_matches = get_node_matchLength(new_graph, t_path)
		#print >> sys.stderr, t_path_graph.edges(data=True)
		#print >> sys.stderr, t_path_matches
		totalweight = t_path_graph.size(weight='weight')
		totalmatches = sum(t_path_matches)
		matches2gapratio = int(100*float(totalmatches)/float(totalweight))
		averagealignlength = float(totalmatches)/float(len(t_path)-2)
		#print >> sys.stderr, "[" + " -> ".join(t_path) + "] | total_weight = " + str(totalweight) + " | total_matches = " + str(totalmatches) + " | matches_to_gap_ratio = " + str(matches2gapratio) + " | average_align_length = " + str(averagealignlength)
		#print >> sys.stderr, t_path_graph.edges(data=True)
		#print >> sys.stderr, ""
		paths.append( [ t_path , totalweight , totalmatches , matches2gapratio , averagealignlength ] )
		counter += 1

	print >> sys.stdout , "-- Found " + str(counter) + " tiling paths"
	print >> sys.stdout , "-- Sorting tiling paths by length"
	sorted_paths = [ paths.sort(key=lambda x: float(x[1])) ]

	return( sorted_paths )


def make_list_from_path( querylist ) :

	### Query list element format
	# Gap:	[61252	,	(0:61252)		,	gap		,	61252	,	0]
	#		[length	,	(T_start:Tstop)	, 	"gap" 	, 	length 	, 	0]
	# Object:	[b40-14.HS_iter4_seq3733|+	,	(61252:6463804)	,	93612:7595148	,	-6402552			,	4526208]
	#			[ID|strand					,	(T_start:Tstop)	,	Q_start:Q_stop	,	-(alignment length)	,	matches]

	id_list = []
	seq_list = []

	for CompntId in querylist :
		Id, T_range , Q_range , alignment , matches = CompntId

		if not Q_range == "gap" :
			# Skip alignment gaps in the list
			Id_clean = Id[:-1]
			seq_list.append(Id)
			id_list.append(Id_clean+"+")
			id_list.append(Id_clean+"-")

	return id_list , seq_list


def hit_mu( hits_file , input_format , max_gap_size , rlen , qlen) :
	print >> sys.stderr, "## Merge hits in blocks and select best alignments"
	merged_hits = {}
	unique_hits = {}

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
	if input_format == "paf" :
		original_hits = read_paf( hits_file )
	elif input_format == "coords" :
		original_hits = read_nucmer_coords( hits_file )
	else :
		print >> sys.stdout , "[ERROR] Unknown alignment format"
		print >> sys.stderr , "[ERROR] Unknown alignment format: " + input_format
		sys.exit(1)
	# original_hits[(Tid,Qid)] =
	# 	[
	# 		[ Qid , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ] , ...
	#	]
	print >> sys.stderr, "#### " + str(len(original_hits.keys())) + " hits read"
	print >> sys.stderr, "### Merging hits"

	for entry in sorted(original_hits.keys()) :
		new_merged_hit = 12*["-"]
		Tid , Qid = entry
		#print >> sys.stderr, entry

		map_graph = make_map_graph( original_hits[entry], max_gap_size )
		try :
			longest_merged_path = nx.dag_longest_path(map_graph, weight='match')
			longest_graph = get_subgraph_from_path_tuples( map_graph , longest_merged_path )
			longest_graph_matches = longest_graph.size(weight='match')

		except :
			print >> sys.stderr, "[ERROR] Unable to merge hits"
			print >> sys.stderr, "[ERROR] Unable to merge hits of " + Qid + " on " + Tid + " (" + str(len(original_hits[entry])) + " hits)"
			print >> sys.stderr, sorted(map_graph.edges.data())
			print >> sys.stderr, nx.find_cycle(map_graph)
			exit(5)

		#print >> sys.stdout, "---- Longest Path : " + "\t".join(str(x) for x in longest_merged_path)

		new_merged_hit[0] = Qid
		new_merged_hit[1] = qlen[Qid]
		new_merged_hit[4] = "+"
		new_merged_hit[5] = Tid
		new_merged_hit[6] = rlen[Tid]

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
		if absQid not in merged_hits : merged_hits[absQid]=[]
		merged_hits[absQid].append(new_merged_hit)

	unique_temp_file = open( hits_file + ".uniq" , "w")

	print >> sys.stderr, "### Selecting best alignment"

	for Qid in sorted(merged_hits.keys()) : # For each absolute id find the best alignment
		merged_list = sorted(merged_hits[Qid] , key=lambda item: int(item[9]) , reverse = True )
		#print >> sys.stderr, merged_list
		unique_hits[Qid] = merged_list[0]
		print >> unique_temp_file , "\t".join(str(x) for x in merged_list[0])

	return unique_hits


def get_gap_bed( fasta_file ) :
	#print >> sys.stderr, fasta_file
	#print >> sys.stderr, os.path.exists(fasta_file)
	out_file_name = remove_extension(fasta_file , ".fasta") + ".gap.bed"
	out_file = open(out_file_name , 'w')
	with open(fasta_file) as handle:
		for record in SeqIO.parse(handle, "fasta"):
			#print >> sys.stderr, record
			#print >> sys.stderr, str(record.seq)
			for match in re.finditer('N+', str(record.seq)):
				start = str(match.start())
				stop = str(match.end())
				if int(stop) - int(start) >= 20 :
					print >> out_file , "\t".join([ record.id, start , stop ])
	out_file.close()
	return out_file_name


def invert_regions( original_dict , seq_length_db ) :
	inverted_dict = {}
	for seq_id in original_dict :
		inverted_dict[seq_id] = []
		start = 0
		for chunk in sorted(original_dict[seq_id]) :
			if not int(chunk[1]) == 0 :
				inverted_dict[seq_id].append([ seq_id , start , int(chunk[1]) -1 ])
				start = int(chunk[2]) + 1
		if not start == seq_length_db[seq_id] :
			inverted_dict[seq_id].append([ seq_id , start , int(seq_length_db[seq_id]) ])

	return inverted_dict


def nearest_to_gap( start , stop , gap ) :
	type, gap_start , gap_stop = gap
	start_quad_distance = min( ( (int(gap_start) - int(start)) ** 2 ) , ( (int(gap_stop) - int(start)) ** 2 ) )
	stop_quad_distance = min( ( (int(gap_start) - int(stop)) ** 2 ) , ( (int(gap_stop) - int(stop)) ** 2 ) )
	if start_quad_distance <= stop_quad_distance :
		return start
	else:
		return stop


def find_nearest_gap( pos , gap_list , seq_len , direction , max_distance ) :
	# Search if falling in an ungapped region or in the gap between 2 of them
	#	-> if corresponding to start or end of a sequence -> report extremity
	#	-> if falling in gap -> report gap
	#	-> if falling sequence
	#		-> shorter distance than max_distance -> report gap
	#		-> longer distance than max_distance -> report sequence

	print >> sys.stderr , "# Searching a gap for: " + str( pos )

	if direction == "upstream" :
		sort_direction = True
		target = ""
	elif direction == "downstream" :
		sort_direction = False
		target = ""
	else :
		prev_sq_distance = int(pos) ** 2
		target = ""
		sort_direction = False

	available_gaps = []
	available_gaps.append(["extremity" , 0 , 0 ])
	available_gaps.append(["extremity" , int(seq_len) + 1 , int(seq_len) + 1 ])
	for gap in gap_list :
		seq_id , gap_start , gap_stop = gap
		available_gaps.append([ "gap" , int(gap_start) , int(gap_stop) ])
	#print >> sys.stderr, "## Available gaps"
	#print >> sys.stderr, available_gaps

	for gap in sorted( available_gaps , key = lambda x: x[1] , reverse=sort_direction ) :
		#search the nearest gap in the list, if not found in the list use extremity
		#print >> sys.stderr, gap
		type , gap_start , gap_stop = gap
		if direction == "upstream" :
			# Gap sorted from the end of the sequence to the beginning
			#print >> sys.stderr, "## Search Upstream"
			if int(gap_start) <= int(pos)  :
				print >> sys.stderr , "## Nearest gap found: [ " + str(gap_start) + " , " + str(gap_stop) + "]"
				# the first gap that do start before pos is the nearest
				if int(gap_start) <= int(pos) <= int(gap_stop) :
					# pos falls within a gap
					# save the gap as target and break the loop
					target = [ "gap" , int(gap_start) , int(gap_stop) ]
				else :
					# Test distance
					if ( (int(gap_start) - int(pos)) ** 2 >= max_distance ** 2 ) and ( (int(gap_stop) - int(pos)) ** 2 >= max_distance ** 2 ) :
						# nearest gap is too distant
						# Report sequence and break
						target = [ "sequence" , int(gap_stop) +1 , int(pos) ]
					else :
						# save the gap as target and break the loop
						target = [ "gap" , int(gap_start) , int(gap_stop) ]
				break
			else :
				continue

		elif direction == "downstream" :
			# Gap sorted from the start of the sequence to the end
			#print >> sys.stderr, "## Search Downstream"
			try :
				int(gap_stop) >= int(pos)
			except :
				print >> sys.stderr, gap_stop
				print >> sys.stderr, pos
			if int(gap_stop) >= int(pos)  :
				print >> sys.stderr , "## Nearest gap found: [ " + str(gap_start) + " , " + str(gap_stop) + "]"
				# the first gap that do ends after pos is the nearest
				if int(gap_start) <= int(pos) <= int(gap_stop) :
					# pos falls within a gap
					# save the gap as target and break the loop
					target = [ "gap" , int(gap_start) , int(gap_stop) ]
				else :
					# Test distance
					if ( (int(gap_start) - int(pos)) ** 2 >= max_distance ** 2 ) and ( (int(gap_stop) - int(pos)) ** 2 >= max_distance ** 2 ) :
						# nearest gap is too distant
						# Report sequence and break
						target = [ "sequence" , int(gap_stop) +1 , int(pos) ]
					else :
						# save the gap as target and break the loop
						target = [ "gap" , int(gap_start) , int(gap_stop) ]
				break
			else :
				continue
		else :
			# direction == "both"
			#print >> sys.stderr, "## Search both directions"
			# Evaluate the distance from each gap, starting from the beginning of the sequence, stop when found the nearest
			sq_distance = min( (int(gap_stop) - int(pos))**2 , (int(gap_start) - int(pos))**2 )
			#print >> sys.stderr, "## Previous distance: " + str(prev_sq_distance) + " - Actual gap distance: " + str(sq_distance)

			if prev_sq_distance >= sq_distance :
				target = [ "gap" , int(gap_start) , int(gap_stop) ]
				prev_sq_distance = sq_distance
			else :
				# found the nearest gap
				# test the distance
				#print >> sys.stderr, "Max square distance: " + str(max_distance ** 2)
				#print >> sys.stderr, "Min gap start square distance: " + str((target[1] - int(pos)) ** 2)
				#print >> sys.stderr, ( (target[1] - int(pos)) ** 2 >= max_distance ** 2 )
				#print >> sys.stderr, "Min gap stop square distance: " + str((target[2] - int(pos)) ** 2)
				#print >> sys.stderr, ( (target[2] - int(pos)) ** 2 >= max_distance ** 2 )
				if ( (target[1] - int(pos)) ** 2 >= max_distance ** 2 ) and ( (target[2] - int(pos)) ** 2 >= max_distance ** 2 ) :
					# nearest gap is too distant
					# Report sequence instead
					target = [ "sequence" , int(pos) , int(gap_start)]
				break

	print >> sys.stderr , "# Scan of gaps completed"
	if ( target == "" ) :
		print >> sys.stderr , "## No good region found, test extremity"
		# Search for nearest gap got to the end of the sequence without finding a proper region gap -> test if extremity is good
		if direction == "upstream" :
			gap_start , gap_stop = 0 , 0
			# test distance
			if ( (int(gap_start) - int(pos)) ** 2 >= max_distance ** 2 ) and ( (int(gap_stop) - int(pos)) ** 2 >= max_distance ** 2 ) :
				target = [ "sequence" , int(gap_stop) + 1 , int(pos) ]
			else :
				target = [ "extremity" , int(gap_start) , int(gap_stop) ]
		else :
			gap_start , gap_stop = int(seq_len) + 1 , int(seq_len) + 1
			# test distance
			if ( (int(gap_start) - int(pos)) ** 2 >= max_distance ** 2 ) and ( (int(gap_stop) - int(pos)) ** 2 >= max_distance ** 2 ) :
				target = [ "sequence" , int(pos) , int(gap_stop) ]
			else :
				target = [ "extremity" , int(gap_start) , int(gap_stop) ]

	else :
		if ( target[1] <= 1 ) or ( target[1] >= int(seq_len) ) :
			if ( (target[1] - int(pos)) ** 2 >= max_distance ** 2 ) and ( (target[2] - int(pos)) ** 2 >= max_distance ** 2 ) :
				target = [ "sequence" , min(target[1] , int(pos)) , max( target[2] ,int(pos) ) ]
			else :
				target[0] = "extremity"

	print >> sys.stderr , "# Final target"
	print >> sys.stderr , target
	print >> sys.stderr , "#####"
	return target


def translate_gap( target , translation_dict , old_seq_len_db ) :
		print >> sys.stderr, "Translate gap coordinates"
		# convert to old coordinates and return both
		target_type, target_start , target_stop = target
		old_target_start_seq_id , old_target_start = translate_coords( int(target_start) , translation_dict )
		old_target_stop_seq_id , old_target_stop = translate_coords( int(target_stop) , translation_dict )
		old_seq = "None"
		target_old = "None"
		# Check if both hit inside the same sequence
		if old_target_start_seq_id == old_target_stop_seq_id :
			if old_target_start_seq_id == "NONE" :
				print >> sys.stderr, "## No corrensponding gap"
			else :
				old_seq = old_target_start_seq_id
				target_old = ["gap" , min( int(old_target_start) , int(old_target_stop)) , max( int(old_target_start) , int(old_target_stop)) ]
				print >> sys.stderr, "## corrensponding gap:"
				print >> sys.stderr, [old_seq , target_old]
		else :
			print >> sys.stderr, "[WARNING] The gap selected seems to span two legacy sequences. Check your inputs"
			print >> sys.stderr, [old_target_start_seq_id , old_target_start , old_target_stop_seq_id , old_target_stop]

		return old_seq , target_old


def translate_coords( pos , offset_dict ) :
	new_chr = ""
	offset = ""
	direction = ""

	# Find correct interval
	for element in sorted(offset_dict.keys()) :
		if int(element[0]) <= int(pos) <= int(element[1]) :
			new_chr , offset , direction = offset_dict[element]

	if direction == "+" :
		new_pos = int(pos) + offset
	elif direction == "-" :
		new_pos = -int(pos) + offset
	else :
		new_chr = "NONE"
		new_pos = "NONE"

	return new_chr , new_pos


def get_adjacent_gap( actual_gap , gap_list , direction ) :
	new_gap = ""
	get_next = False
	if direction == "" :
		for element in sorted(gap_list) :
			if actual_gap == element :
				get_next = True
			else :
				if get_next :
					new_gap = element
					break
	else :  # direction == "upstream"
		for element in sorted(gap_list , reverse=True) :
			if actual_gap == element :
				get_next = True
			else :
				if get_next :
					new_gap = element
					break
	return new_gap


def longest_hit_path( paf_file , max_gap_size) :
	print >> sys.stdout, "- Merge hits in blocks and select best alignments"
	original_hits = {}

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

	id=0

	print >> sys.stdout, "-- Read mapping hits"
	for line in open(paf_file) :
		id+=1
		if line[0] == "#" or line.rstrip() == "": continue ;
		if line[0] == "" or line.rstrip() == "": continue ;
		# subset only the necessary columns 1 -> 11
		Qid , Qlen , Qstart , Qstop , strand , Tid , Tlen , Tstart , Tstop , matches , hitLen = line.rstrip().split("\t")[0:11]

		if ( Tid , Qid ) not in original_hits :
			# create a new entry
			original_hits[(Tid,Qid)] = []

		# Add a new range to the db
		original_hits[(Tid,Qid)].append([ id , Tstart , Tstop , Qstart , Qstop , matches , hitLen ])
	print >> sys.stdout, "--- " + str(id) + " hits read, " + str(len(original_hits.keys())) + " Query to reference couples"

	merged_temp_file = open( paf_file + ".merged" , "w")
	print >> sys.stdout, "-- Merging hits"

	for entry in sorted(original_hits.keys()) :
		new_merged_hit = 12*["-"]
		Tid , Qid = entry
		#print >> sys.stderr, entry

		map_graph = make_map_graph( original_hits[entry], max_gap_size )
		try :
			longest_merged_path = nx.dag_longest_path(map_graph, weight='align')
			longest_graph = get_subgraph_from_path_tuples( map_graph , longest_merged_path )
			#longest_graph_matches = longest_graph.size(weight='match')
			return longest_graph
		except :
			print >> sys.stderr, "#### Error merging hits of " + Qid + " on " + Tid + " (" + str(len(original_hits[entry])) + " hits)"
			print >> sys.stderr, sorted(map_graph.edges.data())
			print >> sys.stderr, nx.find_cycle(map_graph)
			exit(5)


def gfx2bed_print( gfx_file , outfile_name ) :
	outfile=open(outfile_name , 'w')
	for line in open( gfx_file ) :
		if line[0] == "#" or line.rstrip()=="" :
			continue
		el = line.rstrip().split("\t")
		print >> outfile , "\t".join([ str(el[0]) , str( int(el[3]) - 1 ) , str(el[4]) ] )
	outfile.close()
	return(outfile_name)


def write_bed(bed_db, outfile_name , compressed = True) :
	if compressed :
		out_file = gzip.open(outfile_name, 'wb')
	else :
		out_file = open(outfile_name, 'w')

	for chr in sorted(bed_db.keys()):
		for start in sorted(bed_db[chr].keys()):
			for stop in sorted(bed_db[chr][start].keys()):
				print >> out_file , "\t".join(str(x) for x in bed_db[chr][start][stop])
	out_file.close()


def read_bed( bed_file ) :
	bed_db = {}

	try:
		for line in gzip.open( bed_file ) :
			el = line.rstrip().split("\t")
			if el[0] not in bed_db :
				bed_db[el[0]] = {}
			if int(el[1]) not in bed_db[el[0]] :
				bed_db[el[0]][int(el[1])] = {}
			bed_db[el[0]][int(el[1])][int(el[2])] = [ str(el[0]) , int(el[1]) , int(el[2]) ]
			# bed_db[chr][start][stop] = [ chr , start , stop ]

	except IOError:
		for line in open( bed_file ) :
			el = line.rstrip().split("\t")
			if el[0] not in bed_db :
				bed_db[el[0]] = {}
			if int(el[1]) not in bed_db[el[0]] :
				bed_db[el[0]][int(el[1])] = {}
			bed_db[el[0]][int(el[1])][int(el[2])] = [ str(el[0]) , int(el[1]) , int(el[2]) ]

	return bed_db


def read_bed_sorted_list( bed_file ) :
	bed_db = {}
	component = 0
	try:
		for line in gzip.open( bed_file ) :
			component += 1
			el = line.rstrip().split("\t")
			bed_db[component] = el

	except IOError:
		for line in open( bed_file ) :
			component += 1
			el = line.rstrip().split("\t")
			bed_db[component] = el
	# Output format:
	# bed_db[num][ chrom , chromStart(0->) , chromEnd(1->) , name , score , strand , ... ]
	return bed_db


def read_bed_by_sequence(bed_file) :
	bed_db = {}

	try:
		for line in gzip.open( bed_file ) :
			el = line.rstrip().split("\t")
			if el[0] not in bed_db :
				bed_db[el[0]] = []
			bed_db[el[0]] = el

	except IOError:
		for line in open( bed_file ) :
			el = line.rstrip().split("\t")
			if el[0] not in bed_db :
				bed_db[el[0]] = {}
			bed_db[el[0]] = el

	return bed_db


def bed_db_to_range_list( bed_db , seq_id ) :
	range_list = []
	if not bed_db == {} :
		for start in sorted(bed_db[seq_id].keys()) :
			for stop in sorted(bed_db[seq_id][start].keys()) :
				range_list.append( [ start , stop ] )

	return sorted(range_list)


def write_coverage_bed( bam_file , chunk_db , seq_list , bedtools_path , samtools_path ) :
	new_chunk_db = dict(chunk_db)

	if bedtools_path == "" :
		bedtools_search=subprocess.Popen( "which bedtools" , shell=True, stdout=subprocess.PIPE )
		bedtools_command , error = bedtools_search.communicate()
		bedtools_command = bedtools_command.rstrip()
	else :
		bedtools_command = bedtools_path + "/bedtools"

	if samtools_path == "" :
		samtools_search=subprocess.Popen( "which samtools" , shell=True, stdout=subprocess.PIPE )
		samtools_command , error = samtools_search.communicate()
		samtools_command = samtools_command.rstrip()
	else :
		samtools_command = samtools_path + "/samtools"

	if not os.path.exists(bedtools_command) :
		print >> sys.stderr , "[ERROR] wrong or no path to bedtools"
		sys.exit(1)
	if not os.path.exists(samtools_command) :
		print >> sys.stderr , "[ERROR] wrong or no path to samtools"
		sys.exit(1)

	for chr in sorted( seq_list ) :
		coverage_bam_file = new_chunk_db["sequences"][chr]["folder"] + "/" + chr + ".cov.bam"
		samtools_view_err_file_1 = new_chunk_db["sequences"][chr]["folder"] + "/." + chr + ".samtools_view.err1"
		samtools_view_err_file_2 = new_chunk_db["sequences"][chr]["folder"] + "/." + chr + ".samtools_view.err2"
		samtools_view_out_file_2 = new_chunk_db["sequences"][chr]["folder"] + "/." + chr + ".samtools_view.out2"
		coverage_bed_file = new_chunk_db["sequences"][chr]["folder"] + "/" + chr + ".cov.single_base.bed.gz"
		coverage_err_file = new_chunk_db["sequences"][chr]["folder"] + "/." + chr + ".cov.single_base.err"
		coverage_signal_file = new_chunk_db["sequences"][chr]["folder"] + "/" + chr + ".cov.txt.gz"
		coverage_ranges_bed_file = new_chunk_db["sequences"][chr]["folder"] + "/" + chr + ".cov.bed.gz"
		new_chunk_db["sequences"][chr]["coverage_bam"] = coverage_bam_file
		new_chunk_db["sequences"][chr]["coverage_bed.single_base"] = coverage_bed_file
		new_chunk_db["sequences"][chr]["coverage_file"] = coverage_signal_file
		new_chunk_db["sequences"][chr]["coverage_bed"] = coverage_ranges_bed_file
		print >> sys.stderr , "#### " + chr + " -> " + new_chunk_db["sequences"][chr]["coverage_file"]
		# Sequence length
		genome_length_file_name = new_chunk_db["sequences"][chr]["folder"] + "/" + chr + ".len"
		new_chunk_db["sequences"][chr]["length_file"] = genome_length_file_name
		genome_length_file = open(genome_length_file_name , 'w')
		print >> genome_length_file , chr + "\t" + str(chunk_db["sequences"][chr]["length"])
		genome_length_file.close()
		# Subset bam
		command_line = samtools_command + " view " + bam_file + " " + chr + " 2> " + samtools_view_err_file_1 + " | " + samtools_command + " view -b -o " + coverage_bam_file + " -T " + chunk_db["sequences"][chr]["fasta_file"] + " 2> " + samtools_view_err_file_2 + " > " + samtools_view_out_file_2
		print >> sys.stderr, "##### Running command line: " + command_line
		subset_process = subprocess.Popen( command_line , shell=True, stdout=subprocess.PIPE)
		output, error = subset_process.communicate()
		# Indexing
		index_bam( coverage_bam_file , samtools_path)
		# Run coverage
		command_line = bedtools_command + " genomecov -d -ibam " + coverage_bam_file + " -g " + genome_length_file_name + " 2> " + coverage_err_file + " | gzip -9 "
		print >> sys.stderr, "##### Running command line: " + command_line
		coverage_bed = open(coverage_bed_file , 'w')
		coverage_process = subprocess.Popen( command_line , shell=True, stdout=coverage_bed)
		output, error = coverage_process.communicate()
		coverage_bed.close()
		# load coverage
		print >> sys.stderr, "##### Converting coverage BED file to signal"
		coverage_db = read_coverage_bed( coverage_bed_file , chunk_db )
		# Convert
		write_signal_file( coverage_db[chr] , new_chunk_db["sequences"][chr]["coverage_file"] )
		# Write ranged bed
		signal2category_range_bed( coverage_db[chr] , chr, new_chunk_db["sequences"][chr]["coverage_bed"] )

	return new_chunk_db


def split_coverage_bed( coverage_bed_file , chunk_db ) :
	new_chunk_db = dict(chunk_db)

	print >> sys.stderr , '### Loading coverage file, this may take a while'
	coverage_db = read_coverage_bed( coverage_bed_file , new_chunk_db )
	print >> sys.stderr , '### Splitting coverage by sequence'
	for chr in coverage_db.keys() :
		new_chunk_db["sequences"][chr]["coverage_file"] = new_chunk_db["sequences"][chr]["folder"] + "/" + chr + ".cov.txt.gz"
		print >> sys.stderr , chr + " -> " + new_chunk_db["sequences"][chr]["coverage_file"]
		write_signal_file( coverage_db[chr] , new_chunk_db["sequences"][chr]["coverage_file"] )

	return new_chunk_db


def read_coverage_bed( coverage_bed , chunk_db) :
	coverage_signal = {}

	for line in gzip.open( coverage_bed ) :
		el = line.rstrip().split("\t")
		if el[0] not in coverage_signal :
			coverage_signal[el[0]] = [0]*int(chunk_db["sequences"][el[0]]["length"])
		coverage_signal[el[0]][int(el[1])-1] = el[2]

	return coverage_signal


def bed2signal( bed_db_info , sequence_length ):
	signal = {}

	for chr in sorted(bed_db_info.keys()):
		if chr not in signal :
			chr_len = sequence_length[chr]
			signal[chr] = [0]*chr_len
		for start in bed_db_info[chr].keys():
			stop = start+1
			signal[chr][1][start]=int(bed_db_info[chr][start][stop][3])

	return signal


def signal2single_base_bed( signal_db , bed_file ) :
	with gzip.open(bed_file , 'wb') as out_file :
		for chr in sorted(signal_db.keys()):
			start = -1
			for value in sorted(signal_db[chr]):
				start += 1
				stop = start+1
				out_file.write("\t".join(str(x) for x in [ chr, start, stop, value ]))


def range_bed_file2signal( bed_file_gz ) :
	signal = {}

	for line in gzip.open(bed_file_gz , 'rb') :
		chr, start, stop, value = line.rstrip().split("\t")
		new_seq = value * ( int(stop) - int(start) )
		if chr not in signal:
			signal[chr] = new_seq
		else :
			signal[chr] += new_seq

	return signal


def signal2category_range_bed(signal, chr, bed_file) :
	out_file = gzip.open(bed_file , 'wb')
	value = ""
	for pos in range(len(signal)) :
		call = signal[pos]
		if value == "" :
			# Start
			value = call
			start = 0
			stop = 1
		elif not value == call :
			# New value
			# Clone previous range, print it and open a new one
			print >> out_file, "\t".join( str(x) for x in [ chr, start, stop, value ] )
			value = call
			start = pos
			stop = pos + 1
		else :
			stop = pos + 1
	print >> out_file, "\t".join( str(x) for x in [ chr, start, stop, value ] )
	out_file.close()

	return bed_file


def signal2numeric_range_bed(signal, chr, bed_file, tolerance = 0) :
	out_file = gzip.open(bed_file , 'wb')
	value = ""
	for pos in range(len(signal)) :
		call = signal[pos]
		if value == "" :
			# Start
			value = call
			start = 0
			stop = 1
		elif not ( float(value) - tolerance ) <= float(call) <= ( float(value) + tolerance ) :
			# New value
			# Clone previous range, print it and open a new one
			print >> out_file, "\t".join( str(x) for x in [ chr, start, stop, value ] )
			value = call
			start = pos
			stop = pos + 1
		else :
			stop = pos + 1
	out_file.close()

	return bed_file


def moving_average_filter( value_list , window_size ) :
	avg_value_list = [0]*len(value_list)

	for pos in range(len(value_list)) :
		start = max( 0 ,pos-(window_size/2))
		stop = min( len(value_list) , pos + (window_size/2) + 1 )
		values=value_list[start:stop]
		avg_value_list[pos] = int(round(float(sum(values))/float(len(values)),0))

	return avg_value_list


def savitzky_golay_filter( value_list , window_size , order ) :
	new_value = []
	for el in np.round( scipy.signal.savgol_filter( value_list , window_size , order ) , decimals=0 ) :
		new_value.append(int(el))
		#if int(el) > 0 :
		#	new_value.append(int(el))
		#else :
		#	new_value.append(0)
	return new_value


def smooth_coverage( coverage_signal , chr , method , chunk_db ) :
	new_chunk_db = dict(chunk_db)
	smoothed_signal = {}
	# Smooth coverage signal
	if method == "moving_average" :
		smoothed_signal = moving_average_filter(coverage_signal , 5001 )
	elif method == "savitzky_golay" :
		smoothed_signal = savitzky_golay_filter(coverage_signal , 25001 , 3 )
	else :
		print >> sys.stdout , "[ERROR] Smoothing function unknown"

	# Convert to single position bed and write file
	#smoothed_bed_file = new_chunk_db["sequences"][chr]["folder"] + "/" + chr + ".smooth_cov.single_base.bed.gz"
	#smoothed_bed = gzip.open(smoothed_bed_file , 'wb')
	#for pos in range(len(smoothed_signal)) :
	#	print >> smoothed_bed , "\t".join( [ str(chr) , str(pos) , str(pos+1) , str(smoothed_signal[pos]) ] )
	#smoothed_bed.close()
	#new_chunk_db["sequences"][chr]["smooth_cov.single_base"] = smoothed_bed_file
	# Convert to range bed
	range_bed_file_name = new_chunk_db["sequences"][chr]["folder"] + "/" + chr + ".smooth_cov.bed.gz"
	new_chunk_db["sequences"][chr]["smooth_coverage"] = signal2category_range_bed(smoothed_signal, chr, range_bed_file_name)

	return smoothed_signal , new_chunk_db


def mkdir( path ):
	if not os.path.isdir( path ) :
		os.mkdir( path )


def organize_fasta( fasta_file , sequence_db , dir_path , hap ) :

	new_sequence_db = dict(sequence_db)

	for record in SeqIO.parse(open(fasta_file), "fasta"):
		sequence_id = str(record.id)
		sequence_string = remove_islets_fasta(str(record.seq).upper())
		if not hap == "U" :
			### make folder
			folder_path = dir_path + "/" + str(sequence_id)
			mkdir(folder_path)
			### write new file
			file_path = folder_path + "/" + str(sequence_id) + ".fasta"
			f = open(file_path , 'w')
			print >> f , ">" + sequence_id
			print >> f , sequence_string
			f.close()
		else :
			folder_path = dir_path + "/unplaced"
			mkdir(folder_path)
			file_path = folder_path + "/" + str(sequence_id) + ".fasta"
			f = open(file_path , 'w')
			print >> f , ">" + sequence_id
			print >> f , sequence_string
			f.close()
		### add info to sequence_db
		new_sequence_db["sequences"][sequence_id] = {}
		new_sequence_db["sequences"][sequence_id]["folder"] = folder_path
		new_sequence_db["sequences"][sequence_id]["fasta_file"] = file_path
		new_sequence_db["sequences"][sequence_id]["length"] = len(record)
		new_sequence_db["sequences"][sequence_id]["hap"] = hap

	return new_sequence_db


def organize_bed( bed_file , chunk_db , dbtype , clean=True ):
	new_chunk_db = dict(chunk_db)

	db_name = dbtype + "_file"

	bed_db = {}
	#### read bed file
	for line in open(bed_file) :
		record = line.rstrip().split("\t")
		chr = record[0]
		start = int(record[1])
		stop = int(record[2])

		if chr not in bed_db : bed_db[chr] = {}
		if start not in bed_db[chr] : bed_db[chr][start] = {}

		bed_db[chr][start][stop]=record

	#### write separate sorted files
	for fasta_id in sorted(bed_db.keys()) :
		### Get path
		fasta_path = new_chunk_db["sequences"][fasta_id]["fasta_file"]
		### write bed file
		file_path = str(fasta_path.rstrip(".fasta")) + "." + str(dbtype) + ".bed.gz"
		bed_db_to_print = {fasta_id: bed_db[fasta_id]}
		write_bed(bed_db_to_print, file_path)

		### update chunk db with files path
		new_chunk_db["sequences"][fasta_id][db_name] = file_path

	#if clean :
	#	os.remove(bed_file)

	return new_chunk_db


def save_status(chunk_db, chunk_db_file, pairing, pairing_file, status_db, status_file):
	json.dump( chunk_db, open( chunk_db_file , 'w' ) , indent=4)
	json.dump( pairing, open( pairing_file , 'w' ) , indent=4 )
	json.dump( status_db, open( status_file , 'w' ) , indent=4 , sort_keys=True)


def compress_file( in_file , out_file="" ) :
	if out_file == "" : out_file = in_file + ".gz"
	with open( in_file , 'rb') as f_in, gzip.open( out_file , 'wb') as f_out:
		shutil.copyfileobj(f_in, f_out)


def decompress_file( in_file , out_file="" ) :
	if out_file == "" : out_file = remove_extension(in_file , ".gz")
	with gzip.open( in_file , 'rb') as f_in, open( out_file , 'wb') as f_out:
		shutil.copyfileobj(f_in, f_out)


def write_signal_file( signal , signal_file) :
	f = gzip.open( signal_file , 'wb' )
	print >> f , "\t".join(str(x) for x in signal )
	f.close()
	return signal_file


def read_signal_file( signal_file , type = "string" ):
	with gzip.open( signal_file , 'rb' ) as f :
		if type == "string" :
			signal = [ str(x) for x in f.read().split("\t") ]
		elif type == "float" :
			signal = [ float(x) for x in f.read().split("\t") ]
		elif type == "int" :
			signal = [ int(round(x)) for x in f.read().split("\t") ]
		elif type == "mixed_int" :
			# Number as int
			# other as string
			signal = []
			for x in f.read().split("\t") :
				if x.isdigit() :
					signal.append(int(x))
				else :
					signal.append(x)
		elif type == "mixed_float" :
			# Number as int
			# other as string
			signal = []
			for x in f.read().split("\t") :
				if x.isdigit() :
					signal.append(float(x))
				else :
					signal.append(x)
		else :
			print >> sys.stderr , "[ERROR] Unknown type of signal [" + type + "]. Recognized types: [string|int|float|mixed_int|mixed_float]"
			sys.exit(1)

	return signal


def read_coverage_file( coverage_file ) :
	with gzip.open( coverage_file , 'rb') as f :
		coverage_signal = f.read().rstrip().split("\t")

	return coverage_signal


def copy_file( in_file , out_file="" ) :
	if out_file == "" : out_file = os.getcwd() + "/" + os.path.basename(in_file)

	file_copy_command=subprocess.Popen( "cp " + in_file + " " + out_file , shell=True, stdout=subprocess.PIPE )
	out , error = file_copy_command.communicate()


def calculate_clean_median( chunk_db ) :
	rep_and_gap_masked_signal = {}
	rep_and_gap_masked_smoothed = {}
	todo_list = sorted(chunk_db["inputs"]["1_list"] + chunk_db["inputs"]["2_list"])
	for chr in todo_list :
		cov_signal = read_signal_file(chunk_db["sequences"][chr]["coverage_file"], "float")
		smoothed_coverage , chunk_db = smooth_coverage( cov_signal , chr , "savitzky_golay" , chunk_db )

		gap_db = read_bed( chunk_db["sequences"][chr]["gap_file"] )
		cov_signal = mask_regions( cov_signal , gap_db , "G" )
		smoothed_coverage = mask_regions( smoothed_coverage , gap_db , "G" )

		repeats_db = read_bed( chunk_db["sequences"][chr]["repeat_file"] )
		cov_signal = mask_regions( cov_signal , repeats_db , "R" )
		smoothed_coverage = mask_regions( smoothed_coverage , repeats_db , "R" )

		rep_and_gap_masked_signal[chr] = cov_signal
		rep_and_gap_masked_smoothed[chr] = smoothed_coverage

	values = []

	for chr in rep_and_gap_masked_signal.keys() :
		values += [x for x in rep_and_gap_masked_signal[chr] if not (x == "G" or x == "R" ) ]

	median_value = np.median(values)

	return [ median_value , rep_and_gap_masked_signal , rep_and_gap_masked_smoothed ]


def mask_regions( signal_list , bed_db , code ):
	edited_signal = signal_list[:]
	#print >> sys.stderr , signal_list
	for chr in sorted(bed_db.keys()) :
		for start in sorted(bed_db[chr].keys()) :
			for stop in sorted(bed_db[chr][start].keys()) :
				for pos in range(start , stop , 1) :
					edited_signal[pos] = code

	return edited_signal


def write_masked_signal(chunk_db, masked_raw_signal_db, masked_smoothed_signal_db) :
	new_chunk_db = dict(chunk_db)

	for chr in sorted(masked_smoothed_signal_db.keys()) :
		masked_raw_signal_file = chunk_db["sequences"][chr]["folder"] + "/" + chr + ".cov.masked.txt.gz"
		masked_smoothed_signal_file = chunk_db["sequences"][chr]["folder"] + "/" + chr + ".smooth_cov.masked.txt.gz"
		write_signal_file(masked_raw_signal_db[chr], masked_raw_signal_file)
		write_signal_file(masked_smoothed_signal_db[chr], masked_smoothed_signal_file)
		new_chunk_db["sequences"][chr]["masked_raw_signal_file"] = masked_raw_signal_file
		new_chunk_db["sequences"][chr]["masked_signal_file"] = masked_smoothed_signal_file

	return new_chunk_db


def load_smoothed_masked( chunk_db ) :
	masked_signal_db = {}
	smoothed_masked_signal_db = {}
	todo_list = sorted(chunk_db["inputs"]["1_list"] + chunk_db["inputs"]["2_list"])
	for chr in todo_list :
		masked_signal_db[chr] = read_signal_file(chunk_db["sequences"][chr]["masked_raw_signal_file"],  "mixed_int")
		smoothed_masked_signal_db[chr] = read_signal_file(chunk_db["sequences"][chr]["masked_signal_file"], "mixed_int")
		print >> sys.stderr , '#### ' + chr + " loaded"
	return masked_signal_db , smoothed_masked_signal_db


def get_category( chunk_db , masked_signal_db , med_cov ):
	new_chunk_db = dict(chunk_db)
	medianCov = str(med_cov)
	if len(medianCov.split(",")) == 1 :
		print >> sys.stderr , "#### Thresholds (based on an expected coverage of :" + str(medianCov) + "):"
		cov = float(medianCov)
		l_threshold = str( cov * 0.1 )
		m_threshold = str( cov * 0.6 )
		h_threshold = str( cov * 2.4 )
	elif len(medianCov.split(",")) == 3 :
		l_threshold , m_threshold , h_threshold = medianCov.split(",")
	else :
		print >> sys.stdout , "[ERROR] Coverage parameter must be either the expected coverage (one value) or a coma separated list of the 3 levels (min,half,high) "
		print >> sys.stderr , "[ERROR] Coverage parameter must be either the expected coverage (one value) or a coma separated list of the 3 levels (min,half,high) "
		exit(321)

	print >> sys.stderr , "##### 0 : value <= " + l_threshold
	print >> sys.stderr , "##### 1 : " + l_threshold + " < value <= " + m_threshold
	print >> sys.stderr , "##### 2 : " + m_threshold + " < value < " + h_threshold
	print >> sys.stderr , "##### H : value >= " + h_threshold

	for chr in sorted(masked_signal_db.keys()) :
		## Categorize
		new_signal = ["0"]*len(masked_signal_db[chr])
		for pos in range(len(new_signal)) :
			position_value = str(masked_signal_db[chr][pos])
			try :
				numeric_count = float(position_value)
			except ValueError :
				#print >> sys.stderr , "literal"
				new_signal[pos] = position_value
			else :
				if numeric_count <= float(l_threshold) :
					new_signal[pos] = "0"
				elif float(l_threshold) < numeric_count <= float(m_threshold) :
					new_signal[pos] = "1"
				elif float(m_threshold) < numeric_count < float(h_threshold) :
					new_signal[pos] = "2"
				elif numeric_count >= float(h_threshold) :
					new_signal[pos] = "H"
				else :
					new_signal[pos] = "N"
			#finally :
			#	print >> sys.stderr , position_value + " -> " + new_signal[pos]

		## Save coverage levels in BED file
		category_signal = chunk_db["sequences"][chr]["folder"] + "/" + chr + ".cat.txt.gz"
		write_signal_file( new_signal , category_signal)
		category_file = chunk_db["sequences"][chr]["folder"] + "/" + chr + ".cat.bed.gz"
		signal2category_range_bed(new_signal, chr, category_file)

		## Update database
		new_chunk_db["sequences"][chr]["category_signal_file"] = category_signal
		new_chunk_db["sequences"][chr]["category_file"] = category_file
		print >> sys.stderr , '#### ' + chr + " coverage categorized"
	return new_chunk_db


def read_categories( category_file ) :
	categories_db = {}
	for line in gzip.open(category_file) :
		chr, start, stop, value = line.rstrip().split("\t")
		categories_db[start] = [ [chr, start, stop] , value ]

	return categories_db


def plot_histogram( list_db , title , x_label , file_name , medianCov = "" , expected_cov = "" ) :
	data = []
	for key in list_db.keys() :
		data += [ float(x) for x in list_db[key] if isDigit(x) ]

	with PdfPages( file_name ) as pdf:

		plt.style.use("ggplot")

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( data , 50, edgecolor='#E6E6E6', color='#444444')
		if medianCov != "" :
			plt.axvline(x=float(medianCov) , label= "Median value" , color='r' , alpha=.7 )
		if expected_cov != "" :
			plt.axvline(x=float(expected_cov) , label= "Expected value" , color='kb' , alpha=.7 )
		plt.title( title , fontsize=8)
		plt.xlabel( x_label , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()


def write_masked_coverage_bed(raw_masked_db, smoothed_masked_db, chunk_db) :
	new_chunk_db = dict(chunk_db)
	to_check = chunk_db["inputs"]["1_list"] + chunk_db["inputs"]["2_list"]

	for chr in to_check :
		raw_file_name = new_chunk_db["sequences"][chr]["folder"] + "/" + chr + ".cov.masked.bed.gz"
		smoothed_file_name = new_chunk_db["sequences"][chr]["folder"] + "/" + chr + ".smooth_cov.masked.bed.gz"
		y1 = raw_masked_db[chr]
		y2 = smoothed_masked_db[chr]
		y1_clean = [0]*len(y1)
		y2_clean = [0]*len(y2)

		for pos in range(len(y1)) :
			if isDigit(y1[pos]) :
				y1_clean[pos] = float(y1[pos])
				y2_clean[pos] = float(y2[pos])

		new_chunk_db["sequences"][chr]["coverage_clean_signal"] = signal2category_range_bed(y1_clean, chr, raw_file_name)
		new_chunk_db["sequences"][chr]["smoothed_coverage_clean_signal"] = signal2category_range_bed(y2_clean, chr, smoothed_file_name)

		# TODO
		#  run R script in folder to make plot
		#plotCommand = "Rscript " + scriptDirectory + "/clean_coverage.R"
		#plot_log = open(new_chunk_db["sequences"][chr]["folder"] + "/" + chr + ".coverage.clean_plot.log" , "w")
		#plot_err = open(new_chunk_db["sequences"][chr]["folder"] + "/" + chr + ".coverage.clean_plot.err" , "w")
		#plotProcess = subprocess.Popen(plotCommand, shell=True, stdout=plot_log , stderr=plot_err, cwd=new_chunk_db["sequences"][chr]["folder"])
		#output, error = plotProcess.communicate()
		#plot_log.close()
		#plot_err.close()

	return new_chunk_db


def check_all_masked( chunk_db ) :
	if "calculated" not in chunk_db["coverage"] :
		return False
	else :
		all_done = True
		to_check = chunk_db["inputs"]["1_list"] + chunk_db["inputs"]["2_list"]
		for chr in to_check :
			if ( not os.path.exists(chunk_db["sequences"][chr]["folder"] + "/" + chr + ".cov.masked.txt.gz") ) or ( not "masked_signal_file" in chunk_db["sequences"][chr] ) :
				all_done = False
		return all_done


def check_all_categories( chunk_db ) :
	all_done = True
	to_check = chunk_db["inputs"]["1_list"] + chunk_db["inputs"]["2_list"]
	for chr in to_check :
		if ( not os.path.exists(chunk_db["sequences"][chr]["folder"] + "/" + chr + ".smooth_cov.masked.bed.gz") ) or (not "smoothed_coverage_clean_signal" in chunk_db["sequences"][chr]) :
			all_done = False
			print >> sys.stderr , "[ERROR] File for " + chr + " is missing"
	return all_done


def uniquify_paf(paf_file, chunk_db) :
	unique_paf_file = remove_extension(paf_file, ".paf") + ".unique.paf"
	original_hits = []

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

	id=0

	print >> sys.stdout, "-- Read mapping hits"
	for line in open(paf_file, 'r') :
		id+=1
		if line[0] == "#" or line.rstrip() == "": continue ;
		if line[0] == "" or line.rstrip() == "": continue ;
		# subset only the necessary columns 1 -> 11
		Qid , Qlen , Qstart , Qstop , strand , Tid , Tlen , Tstart , Tstop , matches , hitLen = line.rstrip().split("\t")[0:11]

		# Add a new range to the db
		original_hits.append([ str(id) , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ])

	if "sequences" in chunk_db :
		Tid_length = chunk_db["sequences"][Tid]["length"]
		Qid_length = chunk_db["sequences"][Qid]["length"]
	else :
		Tid_length , Qid_length = chunk_db

	original_hits.append(["ChrStart" , -1 , 0 , -1 , 0 , 0 , 0 ])
	original_hits.append(["ChrStop" , Tid_length , Tid_length + 1, Qid_length , Qid_length + 1 , 0 , 0 ])

	map_graph = make_map_graph( original_hits , 100000000 )

	try :
		longest_path_in_map = nx.dag_longest_path(map_graph, weight='align')
		longest_graph = get_subgraph_from_path_tuples( map_graph , longest_path_in_map )

	except :
		print >> sys.stderr, "#### Error uniquifying hits of " + paf_file + ", closed loop found"
		print >> sys.stderr, "#### Map graph edges: "
		print >> sys.stderr, sorted(map_graph.edges.data())
		print >> sys.stderr, "#### Loop: "
		print >> sys.stderr, nx.find_cycle(map_graph)
		exit(5)

	matches = {}

	for edge in longest_graph.edges(data=True) :
		matches[int(edge[0][0])] = [edge[0][0] , edge[1][0] , edge[0][1] , edge[1][1] , edge[2]["id"] , edge[2]["align"] , edge[2]["match"] ]

	f = open( unique_paf_file , 'w' )
	for key in sorted(matches.keys()) :
		if not ( matches[key][4] == "gap" or matches[key][4] == "ChrStart" or matches[key][4] == "ChrStop" ) :
			Tstart , Tstop , Qstart , Qstop , id , align ,  match  = matches[key]
			paf_line = [Qid , Qid_length , Qstart , Qstop , "+" , Tid , Tid_length , Tstart , Tstop , match , align , 255 ]
			print >> f, "\t".join(str(x) for x in paf_line )
	f.close()
	return unique_paf_file


def hits_best_tiling_path(hits_db, seq_length_db) :
	# hits_db[(Tid,Qid)]= [
	#	...
	# 	[ id , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ]
	#	]
	matches = {}
	hit_keys = sorted(hits_db.keys())
	for pair in hit_keys :
		Tid , Qid = pair
		matches[pair] = []
		original_hits = hits_db[pair]
		if Tid in seq_length_db :
			Tid_length = seq_length_db[Tid]
		else :
			print >> sys.stderr , "[ERROR] " + Tid + " sequence is missing in FASTA input file"
			exit(1)
		if Qid in seq_length_db :
			Qid_length = seq_length_db[Qid]
		else :
			print >> sys.stderr , "[ERROR] " + Qid + " sequence is missing in FASTA input file"
			exit(1)

		original_hits.append(["ChrStart" , -1 , 0 , -1 , 0 , 0 , 0 ])
		original_hits.append(["ChrStop" , Tid_length , Tid_length + 1, Qid_length , Qid_length + 1 , 0 , 0 ])

		map_graph = make_map_graph( original_hits , 100000000 )

		try :
			longest_path_in_map = nx.dag_longest_path(map_graph, weight='align')
			longest_graph = get_subgraph_from_path_tuples( map_graph , longest_path_in_map )

		except :
			print >> sys.stderr, "#### Error uniquifying hits, closed loop found"
			print >> sys.stderr, "#### Map graph edges: "
			print >> sys.stderr, sorted(map_graph.edges.data())
			print >> sys.stderr, "#### Loop: "
			print >> sys.stderr, nx.find_cycle(map_graph)
			exit(5)
		edges = []
		for edge in longest_graph.edges(data=True) :
			edges.append([ Tid, int(edge[0][0]) , int(edge[1][0]) , Qid, int(edge[0][1]) , int(edge[1][1]) , int(edge[2]["align"]) , int(edge[2]["match"]) ] )
			# [ Tid , int(Tstart) , int(Tstop) , Qid, int(Qstart) , int(Qstop) , int(align_length) ,  int(match_length) ]
		matches[pair] = sorted(edges)
	return matches


def used_hits_best_tiling_path( hits_file , seq_length_db) :
	print >> sys.stderr, "### Reading Hits"
	hits_db = read_nucmer_coords(hits_file)
	# hits_db[(Tid,Qid)]= [
	#	...
	# 	[ id , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ]
	#	]
	# Like hits_best_tiling_path but filter to keep only the used hits, not the junctions between them
	match_db = {}
	for pair in sorted(hits_db.keys()) :
		Tid , Qid = pair
		print >> sys.stderr, "#### Parsing results for " + Qid + " on " + Tid
		match_db[pair] = []
		original_hits = hits_db[pair]
		hit_regions = []
		for hit in hits_db[pair] :
			id , Tstart , Tstop , Qstart , Qstop , matches , hitLen = hit
			hit_regions.append([int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop)])
		if Tid in seq_length_db :
			Tid_length = seq_length_db[Tid]
		else :
			print >> sys.stderr , "[ERROR] " + Tid + " sequence is missing in FASTA input file"
			exit(1)
		if Qid in seq_length_db :
			Qid_length = seq_length_db[Qid]
		else :
			print >> sys.stderr , "[ERROR] " + Qid + " sequence is missing in FASTA input file"
			exit(1)

		original_hits.append(["ChrStart" , -1 , 0 , -1 , 0 , 0 , 0 ])
		original_hits.append(["ChrStop" , Tid_length , Tid_length + 1, Qid_length , Qid_length + 1 , 0 , 0 ])
		print >> sys.stderr, "##### Found " + str(len(hit_regions)) + " hits"
		print >> sys.stderr, "#### Making the graph"
		map_graph = make_map_graph( original_hits , 100000000 )

		try :
			print >> sys.stderr, "#### Extracting best tiling path"
			longest_path_in_map = nx.dag_longest_path(map_graph, weight='align')
			longest_graph = get_subgraph_from_path_tuples( map_graph , longest_path_in_map )

		except :
			print >> sys.stderr, "#### Error uniquifying hits, closed loop found"
			print >> sys.stderr, "#### Map graph edges: "
			print >> sys.stderr, sorted(map_graph.edges.data())
			print >> sys.stderr, "#### Loop: "
			print >> sys.stderr, nx.find_cycle(map_graph)
			exit(5)
		print >> sys.stderr, "#### Exporting best tiling path hits"
		edges = []
		for edge in longest_graph.edges(data=True) :
			# hit_regions = [ ... , [Tstart , Tstop , Qstart , Qstop] , ... ]
			Tstart = int(edge[0][0])
			Tstop = int(edge[1][0])
			Qstart = int(edge[0][1])
			Qstop = int(edge[1][1])
			align_length = int(edge[2]["align"])
			match_length = int(edge[2]["match"])
			if [Tstart , Tstop , Qstart , Qstop] in hit_regions :
				edges.append([ Tid, int(Tstart) , int(Tstop) , Qid, int(Qstart) , int(Qstop) , int(align_length) , int(match_length) ] )
		match_db[pair] = sorted(edges)
	return match_db


def write_pairs(seq_id, chunk_db, pairing_file) :
	pairs_start_db = {}
	pairs_stop_db = {}
	ranges = []

	pairs_start_db_file = chunk_db["sequences"][seq_id]["folder"] + "/" + seq_id + ".pairing.starts.json.gz"
	pairs_stop_db_file = chunk_db["sequences"][seq_id]["folder"] + "/" + seq_id + ".pairing.stops.json.gz"
	ranges_file = chunk_db["sequences"][seq_id]["folder"] + "/" + seq_id + ".pairing.ranges.json.gz"

	for line in open(pairing_file) :
		Qid, Qstart , Qstop , Tid , Tstart , Tstop , matches , hitLen = line.rstrip().split("\t")
		pairs_start_db[int(Qstart)] = [ [ Qid, int(Qstart) , int(Qstop) ] , [Tid , int(Tstart) , int(Tstop) , int(matches) , int(hitLen) ] ]
		pairs_stop_db[int(Qstop)] = [ [ Qid, int(Qstart) , int(Qstop) ] , [Tid , int(Tstart) , int(Tstop) , int(matches) , int(hitLen) ] ]
		ranges.append([int(Qstart) , int(Qstop)])

	json.dump( pairs_start_db , gzip.open( pairs_start_db_file , 'w' ) , indent=4 )
	json.dump( pairs_stop_db , gzip.open( pairs_stop_db_file , 'w' ) , indent=4 )
	json.dump( ranges , gzip.open( ranges_file , 'w' ) , indent=4 )

	return pairs_start_db_file , pairs_stop_db_file , ranges_file


def read_pairs( chunk_db ) :
	pairs_starts = {}
	pairs_stops = {}
	ranges = {}

	for chr in sorted(chunk_db["sequences"].keys()) :
		if not chunk_db["sequences"][chr]["hap"] == "U" :
			pairs_starts_1 = JSON_key_to_int( json.load( gzip.open( chunk_db["sequences"][chr]["mapping_pairs_starts"] ) ) )
			pairs_stops_1 = JSON_key_to_int( json.load( gzip.open( chunk_db["sequences"][chr]["mapping_pairs_stops"] ) ) )
			ranges_1 = json.load( gzip.open( chunk_db["sequences"][chr]["mapping_ranges"] ) )
			pairs_starts[chr] = pairs_starts_1
			pairs_stops[chr] = pairs_stops_1
			ranges[chr] = ranges_1

	return pairs_starts , pairs_stops , ranges


def check_all_paired( chunk_db ) :
	all_paired = True
	to_check = chunk_db["inputs"]["1_list"] + chunk_db["inputs"]["2_list"]
	for chr in to_check :
		if not os.path.exists( chunk_db["sequences"][chr]["mapping_pairs_starts"] ) : all_paired = False
		if not os.path.exists( chunk_db["sequences"][chr]["mapping_pairs_stops"] )  : all_paired = False
		if not os.path.exists( chunk_db["sequences"][chr]["mapping_ranges"] )       : all_paired = False
	return all_paired


def complement_regions(seq_id , chunk_db , starts_db , ranges_db) :
	complementary_region_db = {}
	mate_edges_to_fill = {}
	mate_id = chunk_db["sequences"][seq_id]["mate_id"]
	seq_len = int(chunk_db["sequences"][seq_id]["length"])
	corr_seq_len = int(chunk_db["sequences"][mate_id]["length"])
	#try :
	#	print >> sys.stderr, ranges_db[seq_id]
	#except :
	#	print >> sys.stderr, ranges_db
	#	exit(51)
	next_start = 0
	start_info = [ [seq_id , 0 , 0 ] , [ mate_id , 0 , 0 , 0 , 0] ]
	for map_range in sorted(ranges_db[seq_id]) :
		#print >> sys.stderr, map_range
		#print >> sys.stderr, starts_db[seq_id][int(map_range[0])]
		#print >> sys.stderr, stops_db[seq_id]
		map_info = starts_db[seq_id][int(map_range[0])]

		if next_start == 0 :
			# first iteration, check if the alignment from 0 on both
			if not ( int(map_range[0]) == 0 and int(map_info[1][1]) == 0 ) :
				if mate_id not in mate_edges_to_fill :
					mate_edges_to_fill[mate_id] = {}
				mate_edges_to_fill[mate_id][0] = 0

		if not int(map_range[0]) == next_start:
			# neither
			#  - the region starts from the beginnig (0)
			#  - the region starts just after the previous one
			# ----> no complementary region of sort exists for those cases
			stop_info = map_info
			#print >> sys.stderr, next_start
			#print >> sys.stderr, start_info
			complementary_region_db[ ( next_start , int(map_range[0]) ) ] = {"start_info":start_info , "stop_info": stop_info}

			if ( int(map_info[1][1]) == 0 ):
				# fist alignment start from 0 on mate
				# add the edge to open database
				if mate_id not in mate_edges_to_fill :
					mate_edges_to_fill[mate_id] = {}
				mate_edges_to_fill[mate_id][0] = 0

		# Update start of next complementary
		next_start = int(map_range[1])
		start_info = map_info


	stop = seq_len
	stop_info = [ [seq_id , seq_len , seq_len ] , [ mate_id , corr_seq_len , corr_seq_len , 0 , 0 ] ]

	if not next_start == stop :
		complementary_region_db[ ( next_start , stop ) ] = {"start_info" : start_info , "stop_info" : stop_info}
		# start_info = region upstream the gap: [ [seq_id , start_upstream , stop_upstream == gap_start] [corr_id , corr_start_upstream , corr_stop_upstream == corr_gap_start , matches , hitLen ] ]
		# stop_info = region downstream the gap: [ [seq_id , start_downstream == gap_stop , stop_downstream] [corr_id , corr_start_downstream == corr_gap_stop , corr_stop_downstream , matches , hitLen ] ]

		if mate_id not in mate_edges_to_fill :
			mate_edges_to_fill[mate_id] = {}
		mate_edges_to_fill[mate_id][corr_seq_len] = corr_seq_len

	return complementary_region_db , mate_edges_to_fill


def translate_bed_sorted_list(bed_sorted_db , agp_db) :
	# bed_sorted_db[id] = [chrom , chromStart(0->) , chromEnd(1->) , name , score , strand , ... ]
	# agp_db[OutSeqID_1]
	# 	#	#	agp_db[OutSeqID_1][int(start_1)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type [=="W"] , CompntId  , CompntStart , CompntEnd , Orientation     ]
	# 	# 	#	agp_db[OutSeqID_1][int(start_2)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type [!="W"] , GapLength , GapType     , Linkage   , LinkageEvidence ]
	# 	#	# start is 1 based!
	translation_db = translate_from_AGP_whole_genome(agp_db)
	# translation_db[old_seq_id][(old_seq_start , old_seq_stop)] = [new_seq_id , offset , direction ]

	translated_bed = []

	for id in bed_sorted_db.keys() :
		feature = bed_sorted_db[id]
		#print >> sys.stderr, "### feature: " + str( feature )
		# feature = [chrom , chromStart(0->) , chromEnd(1->) , name , score , strand , ... ]
		chrom = feature[0]
		chromStart = feature[1]
		chromEnd = feature[2]
		# Check if has more than 3 columns
		if len(feature) > 3:
			# Check if it is a BED6 with strand column
			if (len(feature) > 5 ) and (feature[5] == "+" or feature[5] == "+") :
				strand = feature[5]
				col3_4 = [ feature[3] , feature[4] ]
				if len(feature) > 5 :
					other_cols = feature[6:]
			else:
				strand = ""
				other_cols = feature[3:]
		else :
			strand = ""
			other_cols = []

		# Find the outgoing region
		offset_dict = translation_db[chrom]
		for element in sorted(offset_dict.keys()) :
			if int(chromStart) >= int(element[0]) and int(chromEnd) <= int(element[1]) :
				#print >> sys.stderr, element
				new_chr , offset , direction = offset_dict[element]
				if direction == "+" :
					new_start = min( ( int(chromStart) + offset ) , ( int(chromEnd) + offset ) )
					new_end = max( ( int(chromStart) + offset ) , ( int(chromEnd) + offset ) )
					if not strand == "" :
						new_strand = strand
						translated_bed.append( [ new_chr , new_start , new_end ] + col3_4 + [ strand ] + other_cols )
					else :
						translated_bed.append( [ new_chr , new_start , new_end ] + other_cols )
				elif direction == "-" :
					new_start = min( (-int(chromStart) + offset) , (-int(chromEnd) + offset ) )
					new_end = max( (-int(chromStart) + offset) , (-int(chromEnd) + offset ) )
					if strand == "-" :
						new_strand = "+"
						translated_bed.append( [ new_chr , new_start , new_end ] + col3_4 + [ strand ] + other_cols )
					if strand == "+"  :
						new_strand = "-"
						translated_bed.append( [ new_chr , new_start , new_end ] + col3_4 + [ strand ] + other_cols )
					else :
						# no strand
						#print >> sys.stderr,  [ new_chr , new_start , new_end ] + other_cols
						translated_bed.append( [ new_chr , new_start , new_end ] + other_cols )

	return translated_bed


def translate_region( coords , from_region , to_region , to_seq_len , strand = "+"):
	# by construction ,coords must be inside or match from_region

	start , stop = coords
	from_start , from_stop = from_region
	to_start , to_stop = to_region

	start_delta = start - from_start
	if start_delta < 0 :
		#start_delta = 0
		# Debug error check
		sys.exit(55)
	stop_delta = from_stop - stop
	if stop_delta < 0 :
		#stop_delta = 0
		# Debug error check
		sys.exit(56)

	print >> sys.stderr , "##### Coordinates gap: " + str(start) + ":" + str(stop)
	print >> sys.stderr , "##### Coordinates surrounding region:" + str(from_start) + ":" + str(from_stop)
	print >> sys.stderr , "##### Coordinates destination region:" + str(to_start) + ":" + str(to_stop)
	print >> sys.stderr , "##### Distance from surrounding region start: " + str(start_delta)
	print >> sys.stderr , "##### Distance from surrounding region stop: " + str(stop_delta)

	if to_start == to_stop :
		new_start = to_start
		new_stop = to_stop
		translated_start = to_start
		translated_stop = to_stop

	elif strand == "+" :
		new_start = to_start + start_delta
		try:
			new_stop = to_stop - stop_delta
		except :
			print >> sys.stderr , [ to_stop , stop_delta ]
			exit(521)
		translated_start = min( new_start , to_seq_len )
		translated_stop = max( 0, new_stop)
		print >> sys.stderr , "##### Expected coordinates within destination region: " + str(translated_start) + ":" + str(translated_stop)

	elif strand == "-":
		new_start = to_start + stop_delta
		new_stop = to_stop - start_delta
		translated_start = min( new_start , to_seq_len )
		translated_stop = max( 0, new_stop)
		print >> sys.stderr , "##### Expected coordinates within destination region: " + str(translated_start) + ":" + str(translated_stop)

	else:
		print >> sys.stderr , "### ERROR, unknown strand direction " + strand
		exit(522)



	if translated_start > translated_stop :
		if translated_start > to_stop or translated_stop < to_start :
			print >> sys.stderr , "##### Expected coordinates pointing to an invalid region, discarding information"
			translated_start = ""
			translated_stop = ""
		else :
			half_way = int( (to_stop + to_start) / 2 )
			print >> sys.stderr , "##### Expected coordinates pointing to a negative region. Correcting to " + str(half_way)
			translated_start = half_way
			translated_stop = half_way
	else :
		if (translated_start == translated_stop) and (translated_start == 0 ) and (start == stop) and (start == 0) :
			print >> sys.stderr , "##### Haplotypes beginning regions matching"

	return [ translated_start , translated_stop ]


def check_unmatched_and_open_edges( chunk_db ) :
	passing = True
	for chr in sorted(chunk_db["sequences"].keys()) :
		if not chunk_db["sequences"][chr]["hap"] == "U" :
			if not "open_edges" in chunk_db["sequences"][chr]  :
				passing = False
			else :
				if not os.path.exists(chunk_db["sequences"][chr]["open_edges"]) :
					passing = False
	return passing


def preprocess_gap_list( chr , old_gap_list , open_edges , chr_len , file_out_name ) :
	new_gap_list = []
	chr_gap_info = {}
	# new_gap_list[ int(start) ] = [ chr , start , stop , left_padding , right_padding ]
	for element in old_gap_list :
		new_gap_list.append([ int(element[0]) , int(element[1]) ])
	if not open_edges == {} :
		for edge in open_edges[chr] :
			new_gap_list.append( [ int(edge) , int(edge) ] )
	fileout = open( file_out_name , 'w')
	start = 0
	part_num = -1
	for interval in sorted(new_gap_list) :
		part_num += 1
		if not interval[0] == 0 :
			chr_gap_info[str(part_num)] = [int(start) , int(interval[0]) , "seq"]
			part_num += 1
		chr_gap_info[str(part_num)] = [int(interval[0]) , int(interval[1]) , "gap"]
		print >> fileout , str(chr) + "\t" + str(interval[0]) + "\t" + str(interval[1])
		start = int(interval[1])
	if not int(start) == int(chr_len) :
		part_num += 1
		chr_gap_info[str(part_num)] = [int(start) , int(chr_len), "seq" ]

	fileout.close()
	return new_gap_list , file_out_name , chr_gap_info


def gap_mate_position( seq_id , gap_list , ranges_db , pairs_starts_db , unmatched_regions_db , out_file_name , chunk_db ) :
	gap_db = {}

	# For each gap search for other haplotype matching region and select patching strategy
	for gap in sorted(gap_list) :
		gap_start , gap_stop = gap

		# Search the mapping regions to find if it is mapped or not
		print >> sys.stderr , "### Gap - " + seq_id + ":" + str(gap_start) + "-" + str(gap_stop)
		matching_range = []
		for x in sorted(ranges_db[seq_id]) :
			#print >> sys.stderr , x
			if ( x[0] <= gap_start ) and ( gap_stop <= x[1] ) :
				matching_range.append( [ x[0] , x[1] ]  )

		if matching_range == [] :
			# gap is in an unmatched region
			# retrieve the region and the information of the surrounding matching regions
			print >> sys.stderr , "#### No alignment range is matching"

			unmatched_range = []
			for x in sorted(unmatched_regions_db.keys()) :
				#print >> sys.stderr , x
				if x[0] <= gap_start and x[1] >= gap_stop :
					unmatched_range.append( [ x[0] , x[1] ] )
					break
			#print >> sys.stderr , "#### Unaligned range matching: "
			#print >> sys.stderr , unmatched_range
			#print >> sys.stderr , unmatched_regions_db[x]

			# There should be (one!) unmatched that cover the gap.
			# unmatched_range should be list with just one tuple: unmatching_range==[(start,stop)].
			# TEST uniqueness
			if len(unmatched_range) > 1 :
				print >> sys.stderr , "##### Error: gap in " + seq_id + " extends through multiple unmatching regions. Check pairwise alignment for errors"
				sys.exit(52)
			elif len(unmatched_range) == 1 :
				unmatched_start , unmatched_stop = unmatched_range[0]
			else :
				print >> sys.stderr, "error with unmatched_range"
				exit(521)

			# Extract it and check if patchable
			unmatched_region = unmatched_regions_db[(unmatched_start , unmatched_stop)]
			# Format:
			# unmatched_region == {"start_info" : start_info , "stop_info" : stop_info}
			# [ [Qid, Qstart , Qstop ] , [Tid , Tstart , Tstop , matches , hitLen] ]
			# [ [Qid, Qstart , Qstop ] , [Tid , Tstart , Tstop , matches , hitLen] ]
			corr_unmatched_left_id , corr_unmatched_left_start , corr_unmatched_left_stop , left_match , left_hitLen = unmatched_region["start_info"][1]
			corr_unmatched_right_id , corr_unmatched_right_start , corr_unmatched_right_stop , right_match , right_hitLen = unmatched_region["stop_info"][1]
			corr_unmatched_length = chunk_db["sequences"][corr_unmatched_right_id]["length"]
			# generate paired region
			gap_corr_start , gap_corr_stop = translate_region( (gap_start , gap_stop) , (unmatched_start , unmatched_stop) , (corr_unmatched_left_stop , corr_unmatched_right_start) , corr_unmatched_length )
			#print >> sys.stderr , "#### Unaligned corresponding region: "
			#print >> sys.stderr , [ mate_id , gap_corr_start , gap_corr_stop ]

		else :
			# There should be (one!) match that covers the gap
			# matching_range should be a list with one 2 element list inside: matching_range==[[start,stop]]. TEST
			if len(matching_range) > 1 :
				print >> sys.stderr , "##### Error: gap in " + seq_id + " extends through multiple matching regions. Check pairwise alignment for errors"
				sys.exit(8)
			else :
				match_start , match_stop = matching_range[0]
			# Extract it and check if patchable

			matching_region = pairs_starts_db[seq_id][match_start]
			# Format:
			# matching_region == [ [Qid, Qstart , Qstop ] , [Tid , Tstart , Tstop , matches , hitLen] ]
			match_id , match_start , match_stop = matching_region[0]
			match_id_len = int(chunk_db["sequences"][match_id]["length"])
			corr_matching_id , corr_matching_start , corr_matching_stop , matches , hitLength = matching_region[1]
			gap_corr_start , gap_corr_stop = translate_region( (gap_start , gap_stop) , (match_start , match_stop) , (corr_matching_start , corr_matching_stop) , match_id_len )

		gap_db[int(gap_start)] = [ gap_start , gap_stop, gap_corr_start , gap_corr_stop ]

	# Clean up gap_db from overlapping or crossthreading gaps
	# 	if ( two gaps point to the same region ) OR ( left gap point rightwards that right gap) -> obliterate gap -> gap_corr_start = "" gap_corr_stop = ""
	prev_gap_start = ""
	prev_gap_corr_stop = ""
	for gap_start in sorted(gap_db.keys()) :
		a , b , c , d = gap_db[int(gap_start)]
		gap_start = int(a)
		if c == "" :
			continue
		else :
			gap_corr_start = int(c)
			gap_corr_stop = int(d)

		if prev_gap_start == "" :
			prev_gap_start = gap_start
			prev_gap_corr_stop = gap_corr_stop

		else :
			if gap_corr_start <= prev_gap_corr_stop :
				# ( two gaps point to the same region ) OR ( right gap point to un upstream region than left gap does )
				# Obliterate both gap
				gap_db[gap_start][2] = ""
				gap_db[gap_start][3] = ""
			else :
				prev_gap_start = gap_start
				prev_gap_corr_stop = gap_corr_stop

	json.dump( gap_db , open(out_file_name, 'wb') , indent=4 )
	return gap_db , out_file_name


def get_flanking_region( signal , pos , range , all_gap_db , chr_len , direction = "upstream" ) :
	if direction == "upstream" :
		range_max = min( int(pos) , chr_len )
		range_min = max( 0, (int(pos) - int(range)) )
		for start in sorted(all_gap_db.keys()) :
			gap_start , gap_stop, gap_corr_start , gap_corr_stop = all_gap_db[start]
			if int(gap_stop) > pos :
				break
			else :
				if range_min < int(gap_stop) :
					range_min = int(gap_stop)
	elif direction == "downstream" :
		range_min = max( 0 , int(pos) )
		range_max = min( (int(pos) + int(range)) , chr_len )
		for start in sorted(all_gap_db.keys() , reverse=True) :
			gap_start , gap_stop, gap_corr_start , gap_corr_stop = all_gap_db[start]
			if int(gap_start) < pos :
				break
			else :
				if range_max > int(gap_start) :
					range_max = int(gap_start)
	else :
		print >> sys.stderr, "[ERROR] Unknown direction"
		print >> sys.stdout, "[ERROR] Unknown direction"
		exit(531)
	range_signal = signal[ range_min : range_max ]
	return range_signal , [ range_min , range_max ]


def check_paring_gaps( chunk_db ) :
	complete = True
	for chr in sorted(chunk_db["sequences"].keys()) :
		if not chunk_db["sequences"][chr]["hap"] == "U" :
			if not "processed_gap_db" in chunk_db["sequences"][chr] :
				complete = False
			else :
				if not os.path.exists(chunk_db["sequences"][chr]["processed_gap_db"]) :
					complete = False
	return complete


def extract_sequence_and_signals( seq_id , mate_id , chunk_db , gap_db , mate_gap_db , region_file , seq_fasta , mate_fasta , overhang = 50000) :
	seq_len = int(chunk_db["sequences"][seq_id]["length"])
	mate_seq_len = int(chunk_db["sequences"][mate_id]["length"])

	region_list = []
	print >> sys.stderr , '### Importing category signals'
	category_signal = range_bed_file2signal(chunk_db["sequences"][seq_id]["category_file"])
	mate_category_signal = range_bed_file2signal(chunk_db["sequences"][mate_id]["category_file"])
	print >> sys.stderr , '### Convert coordinates and define surrounding regions'
	for indx in sorted( [ int(x) for x in gap_db.keys() ] ) :
		gap_start = str(indx)
		a , b , c , d = gap_db[gap_start]
		# Postprocess corresponding region and select patching strategy
		gap_start = int(a)
		gap_stop = int(b)
		print >> sys.stderr, "#### Gap: " + str(seq_id) + ":" + str(gap_start) + "-" + str(gap_stop)
		flanking_upstream_content , flanking_upstream_region = get_flanking_region(category_signal[seq_id] , gap_start , overhang , gap_db , seq_len , "upstream")
		flanking_downstream_content , flanking_downstream_region = get_flanking_region(category_signal[seq_id] , gap_stop , overhang , gap_db , seq_len , "downstream")
		if not c == "" :
			gap_corr_start = int(c)
			gap_corr_stop = int(d)
			print >> sys.stderr, "##### Corresponding to >>>> " + str(mate_id) + ":" + str(gap_corr_start) + "-" + str(gap_corr_stop)
			corr_region_content = str(mate_category_signal[mate_id])[int(gap_corr_start):int(gap_corr_stop)]
			upstream_corr_region_content , upstream_corr_region_region = get_flanking_region( mate_category_signal[mate_id] , gap_corr_start , overhang , mate_gap_db , mate_seq_len , "upstream")
			downstream_corr_region_content , downstream_corr_region_region = get_flanking_region( mate_category_signal[mate_id] , gap_corr_stop , overhang , mate_gap_db , mate_seq_len , "downstream")
		else :
			print >> sys.stderr, "##### Corresponding to no region on the alternative haplotype"
			gap_corr_start = ""
			gap_corr_stop = ""
			corr_region_content = ""
			upstream_corr_region_content = ""
			upstream_corr_region_region = [ "" , "" ]
			downstream_corr_region_content = ""
			downstream_corr_region_region = [ "" , "" ]

		#try :
		#	flanking_upstream_content , flanking_upstream_region = get_flanking_region(category_signal[seq_id] , gap_start , overhang , gap_db , seq_len , "upstream")
		#except :
		#	print >> sys.stderr, "Error extracting upstream content of gap"
		#	print >> sys.stderr, gap_start
		#	print >> sys.stderr, max( 0, (gap_start - overhang) )
		#try :
		#	flanking_downstream_content , flanking_downstream_region = get_flanking_region(category_signal[seq_id] , gap_stop , overhang , gap_db , seq_len , "downstream") #category_signal[seq_id][ gap_stop : min( (gap_stop + overhang) , seq_len ) ]
		#except :
		#	print >> sys.stderr, "Error extracting downstream content of gap"
		#	print >> sys.stderr, gap_stop
		#	print >> sys.stderr, min( (gap_stop + overhang) , seq_len )
		#
		#try :
		#	corr_region_content = mate_category_signal[mate_id][gap_corr_start : gap_corr_start]
		#except :
		#	print >> sys.stderr, "Error extracting mate content in place of gap"
		#	print >> sys.stderr, gap_corr_start
		#	print >> sys.stderr, gap_corr_start

		#try :
		#	upstream_corr_region_content , upstream_corr_region_region = get_flanking_region( mate_category_signal[mate_id] , gap_corr_start , overhang , mate_gap_db , mate_seq_len , "upstream") # mate_category_signal[mate_id][max( 0, (gap_corr_start - overhang) ) : gap_corr_start ]
		#except :
		#	print >> sys.stderr, "Error extracting upstream corresponding region"
		#	print >> sys.stderr, gap_corr_start
		#	print >> sys.stderr, max( 0, (gap_corr_start - overhang) )

		#try :
		#	downstream_corr_region_content , downstream_corr_region_region = get_flanking_region( mate_category_signal[mate_id] , gap_corr_stop , overhang , mate_gap_db , mate_seq_len , "downstream") # mate_category_signal[mate_id][ gap_corr_stop : min( (gap_corr_stop + 20000) , corr_seq_len ) ]
		#except :
		#	print >> sys.stderr, "Error extracting downstream corresponding region"
		#	print >> sys.stderr, gap_corr_stop
		#	print >> sys.stderr, min( (gap_corr_stop + 20000) , mate_seq_len )

		print >> sys.stderr, "##### Regions: "
		print >> sys.stderr, "###### Gap side: " + str(seq_id) + ":[" + str(flanking_upstream_region[0]) + "-" + str(flanking_upstream_region[1]) + "][" + str(gap_start) + "-" + str(gap_stop) + "][" + str(flanking_downstream_region[0]) + "-" + str(flanking_downstream_region[1]) + "]"
		if not gap_corr_start == "" :
			print >> sys.stderr, "###### Mate side: " + str(mate_id) + ":[" + str(upstream_corr_region_region[0]) + "-" + str(upstream_corr_region_region[1]) + "][" + str(gap_corr_start) + "-" + str(gap_corr_stop) + "][" + str(downstream_corr_region_region[0]) + "-" + str(downstream_corr_region_region[1]) + "]"

		region_info = {}
		region_info["seq_id"] = seq_id
		region_info["mate_id"] = mate_id
		region_info["gap_region"] = [ gap_start , gap_stop ]
		region_info["flanking_upstream_content"] = flanking_upstream_content
		region_info["flanking_upstream_region"] = flanking_upstream_region
		region_info["flanking_downstream_content"] = flanking_downstream_content
		region_info["flanking_downstream_region"] = flanking_downstream_region
		region_info["gap_corr_region"] = [ gap_corr_start , gap_corr_stop ]
		region_info["gap_corr_region_content"] = corr_region_content
		region_info["upstream_corr_region_content"] = upstream_corr_region_content
		region_info["upstream_corr_region_region"] = upstream_corr_region_region
		region_info["downstream_corr_region_content"] = downstream_corr_region_content
		region_info["downstream_corr_region_region"] = downstream_corr_region_region

		region_db = {}
		region_db["flanking_upstream"] = str( seq_fasta[seq_id][ int(flanking_upstream_region[0]) : int(flanking_upstream_region[1]) ] )
		region_db["flanking_downstream"] = str(seq_fasta[seq_id][ int(flanking_downstream_region[0]) : int(flanking_downstream_region[1]) ])
		if not gap_corr_start == "" :
			region_db["corr_seq"] = str(mate_fasta[mate_id][ int(gap_corr_start) : int(gap_corr_stop) ])
			region_db["corr_seq_upstream"] = str(mate_fasta[mate_id][ int(upstream_corr_region_region[0]) : int(upstream_corr_region_region[1]) ])
			region_db["corr_seq_downstream"] = str(mate_fasta[mate_id][ int(downstream_corr_region_region[0]) : int(downstream_corr_region_region[1]) ])
		else :
			region_db["corr_seq"] = "0"
			region_db["corr_seq_upstream"] = "0"
			region_db["corr_seq_downstream"] = "0"

		chars = {}
		chars["flanking_upstream_content"] = analyze_signal( flanking_upstream_content )
		chars["flanking_downstream_content"] = analyze_signal( flanking_downstream_content )
		chars["gap_corr_region_content"] = analyze_signal( corr_region_content )
		chars["upstream_corr_region_content"] = analyze_signal( upstream_corr_region_content )
		chars["downstream_corr_region_content"] = analyze_signal( downstream_corr_region_content )

		region_info["fasta_sequences"] = region_db
		region_info["sequences_characteristics"] = chars

		region_info["patching_strategy"] = status_to_strategy(region_info["sequences_characteristics"] )

		region_list.append(region_info)
	# Save regions to file
	json.dump( region_list, gzip.open(region_file, 'wb') , indent=4 )

	return region_file


def check_extracted_sequences( chunk_db ) :
	chr_list = []
	for seq_id in sorted(chunk_db["sequences"].keys()) :
		if not chunk_db["sequences"][seq_id]["hap"] == "U" :
			done_file = chunk_db["sequences"][seq_id]["region_of_interest_file"] + ".done"
			if not os.path.exists(chunk_db["sequences"][seq_id]["region_of_interest_file"] + ".done") :
				chr_list.append(seq_id)
	return chr_list


def analyze_signal( signal ):
	# signal codes:
	# "0" - no coverage
	# "1" - haploid coverage
	# "H" - Too high coverage
	# "2" - diploid coverage
	# "G" - Gap
	# "R" - Repeat

	signal_length = float(len(signal))
	# Haploid signal: no read mapping coverage or half than expected
	signal_h = [x for x in signal if x in (0 , "0" , 1 , "1") ]
	# Diploid coverage signal
	signal_2 = [x for x in signal if x in (2 , "2") ]
	# Gap Signal
	signal_g = [x for x in signal if x == "G" ]
	# Repetitive signal: Repeat elements or too high mapping coverage
	signal_r = [x for x in signal if x in ("H" , "R") ]

	# Usable values:
	# DIP: Sequence >75% with coverage within diploid level value
	# GAP: Sequence made >75% of gap
	# REP: Sequence made >75% of repeats
	# HAP: Sequence made >75% of haploid regions
	# OK: dip + rep >75%
	# NO: Unreliable

	if len(signal_g) > (0.75 * signal_length) :
		# the sequence is made mostly of gap
		summary = "GAP"
	elif len(signal_r) > (0.75 * signal_length) :
		# the sequence is made mostly of repeats
		summary = "REP"
	elif len(signal_h) > (0.75 * signal_length) :
		# the sequence is made mostly of haploid regions
		summary = "HAP"
	elif len(signal_2) > (0.75 * signal_length) :
		# the sequence is made mostly of repeats
		summary = "DIP"
	elif len(signal_2) + len(signal_r) > (0.75 * signal_length) :
		summary = "OK"
	else :
		summary = "NO"

	return summary


def collect_files_names( chunk_db , descriptor ) :
	file_list = []
	for seq_id in sorted(chunk_db["sequences"].keys()) :
		if descriptor in chunk_db["sequences"][seq_id] :
			file_list.append(chunk_db["sequences"][seq_id][descriptor])
	return file_list


def status_to_strategy( status_db ) :
	# Input content
	#	status_db["flanking_upstream_content"] = analyze_signal( flanking_upstream_content )
	#	status_db["flanking_downstream_content"] = analyze_signal( flanking_downstream_content )
	#	status_db["gap_corr_region_content"] = analyze_signal( corr_region_content )
	#	status_db["upstream_corr_region_content"] = analyze_signal( upstream_corr_region_content )
	#	status_db["downstream_corr_region_content"] = analyze_signal( downstream_corr_region_content )
	# Each of them may assume the following status:
	# 	DIP: Sequence >75% with coverage within diploid level value
	#	GAP: Sequence made >75% of gap
	#	REP: Sequence made >75% of repeats
	#	HAP: Sequence made >75% of haploid regions
	#	OK: dip + rep >75%
	#	NO: Unreliable
	status_ref = [ status_db["flanking_upstream_content"] , status_db["gap_corr_region_content"] , status_db["flanking_downstream_content"] ]
	status_alt =[ status_db["upstream_corr_region_content"] , status_db["gap_corr_region_content"] , status_db["downstream_corr_region_content"] ]
	status_alt = [ x if not ( x == "HAP" ) else "NO" for x in status_alt ]
	strategy = {}
	#print >> sys.stderr, status_ref
	#print >> sys.stderr, status_alt
	# output
	# Strategy is list of values:
	#	Strategy = {}
	#	#	Strategy["map_gap"] = "value"
	#	#	#	"NONE"  					|	unreliable	|		gap		|	unreliable	|
	#	#	#	"flanking_left"  			|	reliable	|		gap		|	unreliable	|
	#	#	#	"flanking_left:gap:rep" 	|	reliable	|		gap		|	repeat 		|
	#	#	#	"flanking_left:gap:right" 	|	reliable	|		gap		|	reliable	|
	#	#	#	"flanking_right" 			|	unreliable	|		gap		|	reliable 	|
	#	#	#	"flanking_rep:gap:right" 	|	repeat		|		gap		|	unreliable	|
	#	#	#	"flanking_rep:gap:rep" 		|	repeat		|		gap		|	repeat		|
	#	#	#	"hybrid_left:alt:right" 	|	reliable	|	[patch alt]	|	reliable	|
	#	#	#	"hybrid_left:alt" 			|	reliable	|	[patch alt]	|	unreliable	|
	#	#	#	"hybrid_left:alt:rep" 		|	reliable	|	[patch alt]	|	repeat		|
	#	#	#	"hybrid_alt:right" 			|	unreliable	|	[patch alt]	|	reliable	|
	#	#	#	"hybrid_rep:alt:right" 		|	repeat		|	[patch alt]	|	reliable	|
	#	#	#	"hybrid_rep:alt:rep" 		|	repeat		|	[patch alt]	|	repeat		|
	#	#	#	"hybrid_rep:alt" 			|	repeat		|	[patch alt]	|	unreliable	|
	#	#	#	"hybrid_alt:rep" 			|	unreliable	|	[patch alt]	|	repeat		|
	if status_ref == ["DIP","DIP","DIP"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["DIP","DIP","GAP"] : strategy["map_gap"] = "hybrid_left:alt"
	elif status_ref == ["DIP","DIP","REP"] : strategy["map_gap"] = "hybrid_left:alt:rep"
	elif status_ref == ["DIP","DIP","HAP"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["DIP","DIP","OK"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["DIP","DIP","NO"] : strategy["map_gap"] = "hybrid_left:alt"
	elif status_ref == ["DIP","GAP","DIP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["DIP","GAP","GAP"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["DIP","GAP","REP"] : strategy["map_gap"] = "flanking_left:gap:rep"
	elif status_ref == ["DIP","GAP","HAP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["DIP","GAP","OK"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["DIP","GAP","NO"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["DIP","REP","DIP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["DIP","REP","GAP"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["DIP","REP","REP"] : strategy["map_gap"] = "flanking_left:gap:rep"
	elif status_ref == ["DIP","REP","HAP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["DIP","REP","OK"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["DIP","REP","NO"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["DIP","HAP","DIP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["DIP","HAP","GAP"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["DIP","HAP","REP"] : strategy["map_gap"] = "flanking_left:gap:rep"
	elif status_ref == ["DIP","HAP","HAP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["DIP","HAP","OK"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["DIP","HAP","NO"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["DIP","OK","DIP"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["DIP","OK","GAP"] : strategy["map_gap"] = "hybrid_left:alt"
	elif status_ref == ["DIP","OK","REP"] : strategy["map_gap"] = "hybrid_left:alt:rep"
	elif status_ref == ["DIP","OK","HAP"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["DIP","OK","OK"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["DIP","OK","NO"] : strategy["map_gap"] = "hybrid_left:alt"
	elif status_ref == ["DIP","NO","DIP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["DIP","NO","GAP"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["DIP","NO","REP"] : strategy["map_gap"] = "flanking_left:gap:rep"
	elif status_ref == ["DIP","NO","HAP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["DIP","NO","OK"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["DIP","NO","NO"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["GAP","DIP","DIP"] : strategy["map_gap"] = "hybrid_alt:right"
	elif status_ref == ["GAP","DIP","REP"] : strategy["map_gap"] = "hybrid_alt:rep"
	elif status_ref == ["GAP","DIP","HAP"] : strategy["map_gap"] = "hybrid_alt:right"
	elif status_ref == ["GAP","DIP","OK"] : strategy["map_gap"] = "hybrid_alt:right"
	elif status_ref == ["GAP","GAP","DIP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["GAP","GAP","HAP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["GAP","GAP","OK"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["GAP","REP","DIP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["GAP","REP","HAP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["GAP","REP","OK"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["GAP","HAP","DIP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["GAP","HAP","HAP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["GAP","HAP","OK"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["GAP","OK","DIP"] : strategy["map_gap"] = "hybrid_alt:right"
	elif status_ref == ["GAP","OK","REP"] : strategy["map_gap"] = "hybrid_alt:rep"
	elif status_ref == ["GAP","OK","HAP"] : strategy["map_gap"] = "hybrid_alt:right"
	elif status_ref == ["GAP","OK","OK"] : strategy["map_gap"] = "hybrid_alt:right"
	elif status_ref == ["GAP","NO","DIP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["GAP","NO","HAP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["GAP","NO","OK"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["REP","DIP","DIP"] : strategy["map_gap"] = "hybrid_rep:alt:right"
	elif status_ref == ["REP","DIP","GAP"] : strategy["map_gap"] = "hybrid_rep:alt"
	elif status_ref == ["REP","DIP","REP"] : strategy["map_gap"] = "hybrid_rep:alt:rep"
	elif status_ref == ["REP","DIP","HAP"] : strategy["map_gap"] = "hybrid_rep:alt:right"
	elif status_ref == ["REP","DIP","OK"] : strategy["map_gap"] = "hybrid_rep:alt:right"
	elif status_ref == ["REP","DIP","NO"] : strategy["map_gap"] = "hybrid_rep:alt"
	elif status_ref == ["REP","GAP","DIP"] : strategy["map_gap"] = "flanking_rep:gap:right"
	elif status_ref == ["REP","GAP","REP"] : strategy["map_gap"] = "flanking_rep:gap:rep"
	elif status_ref == ["REP","GAP","HAP"] : strategy["map_gap"] = "flanking_rep:gap:right"
	elif status_ref == ["REP","GAP","OK"] : strategy["map_gap"] = "flanking_rep:gap:right"
	elif status_ref == ["REP","REP","DIP"] : strategy["map_gap"] = "flanking_rep:gap:right"
	elif status_ref == ["REP","REP","REP"] : strategy["map_gap"] = "flanking_rep:gap:rep"
	elif status_ref == ["REP","REP","HAP"] : strategy["map_gap"] = "flanking_rep:gap:right"
	elif status_ref == ["REP","REP","OK"] : strategy["map_gap"] = "flanking_rep:gap:right"
	elif status_ref == ["REP","HAP","DIP"] : strategy["map_gap"] = "flanking_rep:gap:right"
	elif status_ref == ["REP","HAP","REP"] : strategy["map_gap"] = "flanking_rep:gap:rep"
	elif status_ref == ["REP","HAP","HAP"] : strategy["map_gap"] = "flanking_rep:gap:right"
	elif status_ref == ["REP","HAP","OK"] : strategy["map_gap"] = "flanking_rep:gap:right"
	elif status_ref == ["REP","OK","DIP"] : strategy["map_gap"] = "hybrid_rep:alt:right"
	elif status_ref == ["REP","OK","GAP"] : strategy["map_gap"] = "hybrid_rep:alt"
	elif status_ref == ["REP","OK","REP"] : strategy["map_gap"] = "hybrid_rep:alt:rep"
	elif status_ref == ["REP","OK","HAP"] : strategy["map_gap"] = "hybrid_rep:alt:right"
	elif status_ref == ["REP","OK","OK"] : strategy["map_gap"] = "hybrid_rep:alt:right"
	elif status_ref == ["REP","OK","NO"] : strategy["map_gap"] = "hybrid_rep:alt"
	elif status_ref == ["REP","NO","DIP"] : strategy["map_gap"] = "flanking_rep:gap:right"
	elif status_ref == ["REP","NO","REP"] : strategy["map_gap"] = "flanking_rep:gap:rep"
	elif status_ref == ["REP","NO","HAP"] : strategy["map_gap"] = "flanking_rep:gap:right"
	elif status_ref == ["REP","NO","OK"] : strategy["map_gap"] = "flanking_rep:gap:right"
	elif status_ref == ["HAP","DIP","DIP"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["HAP","DIP","GAP"] : strategy["map_gap"] = "hybrid_left:alt"
	elif status_ref == ["HAP","DIP","REP"] : strategy["map_gap"] = "hybrid_left:alt:rep"
	elif status_ref == ["HAP","DIP","HAP"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["HAP","DIP","OK"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["HAP","DIP","NO"] : strategy["map_gap"] = "hybrid_left:alt"
	elif status_ref == ["HAP","GAP","DIP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["HAP","GAP","GAP"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["HAP","GAP","REP"] : strategy["map_gap"] = "flanking_left:gap:rep"
	elif status_ref == ["HAP","GAP","HAP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["HAP","GAP","OK"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["HAP","GAP","NO"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["HAP","REP","DIP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["HAP","REP","GAP"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["HAP","REP","REP"] : strategy["map_gap"] = "flanking_left:gap:rep"
	elif status_ref == ["HAP","REP","HAP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["HAP","REP","OK"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["HAP","REP","NO"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["HAP","HAP","DIP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["HAP","HAP","GAP"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["HAP","HAP","REP"] : strategy["map_gap"] = "flanking_left:gap:rep"
	elif status_ref == ["HAP","HAP","HAP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["HAP","HAP","OK"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["HAP","HAP","NO"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["HAP","OK","DIP"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["HAP","OK","GAP"] : strategy["map_gap"] = "hybrid_left:alt"
	elif status_ref == ["HAP","OK","REP"] : strategy["map_gap"] = "hybrid_left:alt:rep"
	elif status_ref == ["HAP","OK","HAP"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["HAP","OK","OK"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["HAP","OK","NO"] : strategy["map_gap"] = "hybrid_left:alt"
	elif status_ref == ["HAP","NO","DIP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["HAP","NO","GAP"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["HAP","NO","REP"] : strategy["map_gap"] = "flanking_left:gap:rep"
	elif status_ref == ["HAP","NO","HAP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["HAP","NO","OK"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["HAP","NO","NO"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["OK","DIP","DIP"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["OK","DIP","GAP"] : strategy["map_gap"] = "hybrid_left:alt"
	elif status_ref == ["OK","DIP","REP"] : strategy["map_gap"] = "hybrid_left:alt:rep"
	elif status_ref == ["OK","DIP","HAP"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["OK","DIP","OK"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["OK","DIP","NO"] : strategy["map_gap"] = "hybrid_left:alt"
	elif status_ref == ["OK","GAP","DIP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["OK","GAP","GAP"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["OK","GAP","REP"] : strategy["map_gap"] = "flanking_left:gap:rep"
	elif status_ref == ["OK","GAP","HAP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["OK","GAP","OK"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["OK","GAP","NO"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["OK","REP","DIP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["OK","REP","GAP"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["OK","REP","REP"] : strategy["map_gap"] = "flanking_left:gap:rep"
	elif status_ref == ["OK","REP","HAP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["OK","REP","OK"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["OK","REP","NO"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["OK","HAP","DIP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["OK","HAP","GAP"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["OK","HAP","REP"] : strategy["map_gap"] = "flanking_left:gap:rep"
	elif status_ref == ["OK","HAP","HAP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["OK","HAP","OK"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["OK","HAP","NO"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["OK","OK","DIP"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["OK","OK","GAP"] : strategy["map_gap"] = "hybrid_left:alt"
	elif status_ref == ["OK","OK","REP"] : strategy["map_gap"] = "hybrid_left:alt:rep"
	elif status_ref == ["OK","OK","HAP"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["OK","OK","OK"] : strategy["map_gap"] = "hybrid_left:alt:right"
	elif status_ref == ["OK","OK","NO"] : strategy["map_gap"] = "hybrid_left:alt"
	elif status_ref == ["OK","NO","DIP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["OK","NO","GAP"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["OK","NO","REP"] : strategy["map_gap"] = "flanking_left:gap:rep"
	elif status_ref == ["OK","NO","HAP"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["OK","NO","OK"] : strategy["map_gap"] = "flanking_left:gap:right"
	elif status_ref == ["OK","NO","NO"] : strategy["map_gap"] = "flanking_left"
	elif status_ref == ["NO","DIP","DIP"] : strategy["map_gap"] = "hybrid_alt:right"
	elif status_ref == ["NO","DIP","REP"] : strategy["map_gap"] = "hybrid_alt:rep"
	elif status_ref == ["NO","DIP","HAP"] : strategy["map_gap"] = "hybrid_alt:right"
	elif status_ref == ["NO","DIP","OK"] : strategy["map_gap"] = "hybrid_alt:right"
	elif status_ref == ["NO","GAP","DIP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["NO","GAP","HAP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["NO","GAP","OK"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["NO","REP","DIP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["NO","REP","HAP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["NO","REP","OK"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["NO","HAP","DIP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["NO","HAP","HAP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["NO","HAP","OK"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["NO","OK","DIP"] : strategy["map_gap"] = "hybrid_alt:right"
	elif status_ref == ["NO","OK","REP"] : strategy["map_gap"] = "hybrid_alt:rep"
	elif status_ref == ["NO","OK","HAP"] : strategy["map_gap"] = "hybrid_alt:right"
	elif status_ref == ["NO","OK","OK"] : strategy["map_gap"] = "hybrid_alt:right"
	elif status_ref == ["NO","NO","DIP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["NO","NO","HAP"] : strategy["map_gap"] = "flanking_right"
	elif status_ref == ["NO","NO","OK"] : strategy["map_gap"] = "flanking_right"
	else : strategy["map_gap"] = "NONE"


	#	#	Strategy["map_alt"] = "value"
	#	#	#	"NONE" :					|	unreliable	|	unreliable	|	unreliable	|
	#	#	#	"alt"						|	unreliable	|	reliable	|	unreliable	|
	#	#	#	"altF"						|	unreliable	|	fill		|	unreliable	|
	#	#	#	"altR"						|	unreliable	|	repeat		|	unreliable	|
	#	#	#	"alt:right"					|	unreliable	|	reliable	|	reliable	|
	#	#	#	"altF:right"				|	unreliable	|	fill		|	reliable	|
	#	#	#	"altR:right"				|	unreliable	|	repeat		|	reliable	|
	#	#	#	"rep:alt:right"				|	repeat		|	reliable	|	reliable	|
	#	#	#	"rep:altF:right"			|	repeat		|	fill		|	reliable	|
	#	#	#	"rep:altR:right"			|	repeat		|	repeat		|	reliable	|
	#	#	#	"left:alt"					|	reliable	|	reliable	|	unreliable	|
	#	#	#	"left:altF"					|	reliable	|	fill		|	unreliable	|
	#	#	#	"left:altR"					|	reliable	|	repeat		|	unreliable	|
	#	#	#	"left:alt:rep"				|	reliable	|	reliable	|	repeat		|
	#	#	#	"left:altF:rep"				|	reliable	|	fill		|	repeat		|
	#	#	#	"left:altR:rep"				|	reliable	|	repeat		|	repeat		|
	#	#	#	"left:alt:right"			|	reliable	|	reliable	|	reliable	|
	#	#	#	"left:altF:right"			|	reliable	|	fill		|	reliable	|
	#	#	#	"left:altR:right"			|	reliable	|	repeat		|	reliable	|
	#	#	#	"rep:alt:rep"				|	repeat		|	reliable	|	repeat		|
	#	#	#	"rep:altF:rep"				|	repeat		|	fill		|	repeat		|
	#	#	#	"rep:altR:rep"				|	repeat		|	repeat		|	repeat		|
	#	#	#	"left:gap:right"			|	reliable	|	gap			|	reliable	|
	#	#	#	"left:gap:rep"				|	reliable	|	gap			|	repeat		|
	#	#	#	"rep:gap:right"				|	repeat		|	gap			|	reliable	|
	#	#	#	"rep:alt"					|	repeat		|	reliable	|	unreliable	|
	#	#	#	"rep:altF"					|	repeat		|	fill		|	unreliable	|
	#	#	#	"rep:altR"					|	repeat		|	repeat		|	unreliable	|
	#	#	#	"alt:rep"					|	unreliable	|	reliable	|	repeat		|
	#	#	#	"altF:rep"					|	unreliable	|	fill		|	repeat		|
	#	#	#	"altR:rep"					|	unreliable	|	repeat		|	repeat		|
	#	#	#	"right"						|	unreliable	|	unreliable	|	reliable	|
	#	#	#	"left"						|	reliable	|	unreliable	|	unreliable	|

	if status_alt == ["DIP","DIP","DIP"] : strategy["map_alt"] = "left:altF:right"
	elif status_alt == ["DIP","DIP","GAP"] : strategy["map_alt"] = "left:altF"
	elif status_alt == ["DIP","DIP","REP"] : strategy["map_alt"] = "left:altF:rep"
	elif status_alt == ["DIP","DIP","OK"] : strategy["map_alt"] = "left:altF:right"
	elif status_alt == ["DIP","DIP","NO"] : strategy["map_alt"] = "left:altF"
	elif status_alt == ["DIP","GAP","DIP"] : strategy["map_alt"] = "left:gap:right"
	elif status_alt == ["DIP","GAP","GAP"] : strategy["map_alt"] = "left"
	elif status_alt == ["DIP","GAP","REP"] : strategy["map_alt"] = "left:gap:rep"
	elif status_alt == ["DIP","GAP","OK"] : strategy["map_alt"] = "left:gap:right"
	elif status_alt == ["DIP","GAP","NO"] : strategy["map_alt"] = "left"
	elif status_alt == ["DIP","REP","DIP"] : strategy["map_alt"] = "left:altR:right"
	elif status_alt == ["DIP","REP","GAP"] : strategy["map_alt"] = "left:altR"
	elif status_alt == ["DIP","REP","REP"] : strategy["map_alt"] = "left:altR:rep"
	elif status_alt == ["DIP","REP","OK"] : strategy["map_alt"] = "left:altR:right"
	elif status_alt == ["DIP","REP","NO"] : strategy["map_alt"] = "left:altR"
	elif status_alt == ["DIP","OK","DIP"] : strategy["map_alt"] = "left:alt:right"
	elif status_alt == ["DIP","OK","GAP"] : strategy["map_alt"] = "left:alt"
	elif status_alt == ["DIP","OK","REP"] : strategy["map_alt"] = "left:alt:rep"
	elif status_alt == ["DIP","OK","OK"] : strategy["map_alt"] = "left:alt:right"
	elif status_alt == ["DIP","OK","NO"] : strategy["map_alt"] = "left:alt"
	elif status_alt == ["DIP","NO","DIP"] : strategy["map_alt"] = "left:gap:right"
	elif status_alt == ["DIP","NO","GAP"] : strategy["map_alt"] = "left"
	elif status_alt == ["DIP","NO","REP"] : strategy["map_alt"] = "left:gap:rep"
	elif status_alt == ["DIP","NO","OK"] : strategy["map_alt"] = "left:gap:right"
	elif status_alt == ["DIP","NO","NO"] : strategy["map_alt"] = "left"
	elif status_alt == ["GAP","DIP","DIP"] : strategy["map_alt"] = "altF:right"
	elif status_alt == ["GAP","DIP","GAP"] : strategy["map_alt"] = "altF"
	elif status_alt == ["GAP","DIP","REP"] : strategy["map_alt"] = "altF:rep"
	elif status_alt == ["GAP","DIP","OK"] : strategy["map_alt"] = "altF:right"
	elif status_alt == ["GAP","DIP","NO"] : strategy["map_alt"] = "altF"
	elif status_alt == ["GAP","GAP","DIP"] : strategy["map_alt"] = "right"
	elif status_alt == ["GAP","GAP","OK"] : strategy["map_alt"] = "right"
	elif status_alt == ["GAP","REP","DIP"] : strategy["map_alt"] = "altR:right"
	elif status_alt == ["GAP","REP","GAP"] : strategy["map_alt"] = "altR"
	elif status_alt == ["GAP","REP","REP"] : strategy["map_alt"] = "altR:rep"
	elif status_alt == ["GAP","REP","OK"] : strategy["map_alt"] = "altR:right"
	elif status_alt == ["GAP","REP","NO"] : strategy["map_alt"] = "altR"
	elif status_alt == ["GAP","OK","DIP"] : strategy["map_alt"] = "alt:right"
	elif status_alt == ["GAP","OK","GAP"] : strategy["map_alt"] = "alt"
	elif status_alt == ["GAP","OK","REP"] : strategy["map_alt"] = "alt:rep"
	elif status_alt == ["GAP","OK","OK"] : strategy["map_alt"] = "alt:right"
	elif status_alt == ["GAP","OK","NO"] : strategy["map_alt"] = "alt"
	elif status_alt == ["GAP","NO","DIP"] : strategy["map_alt"] = "right"
	elif status_alt == ["GAP","NO","OK"] : strategy["map_alt"] = "right"
	elif status_alt == ["REP","DIP","DIP"] : strategy["map_alt"] = "rep:altF:right"
	elif status_alt == ["REP","DIP","GAP"] : strategy["map_alt"] = "rep:altF"
	elif status_alt == ["REP","DIP","REP"] : strategy["map_alt"] = "rep:altF:rep"
	elif status_alt == ["REP","DIP","OK"] : strategy["map_alt"] = "rep:altF:right"
	elif status_alt == ["REP","DIP","NO"] : strategy["map_alt"] = "rep:altF"
	elif status_alt == ["REP","GAP","DIP"] : strategy["map_alt"] = "rep:gap:right"
	elif status_alt == ["REP","GAP","OK"] : strategy["map_alt"] = "rep:gap:right"
	elif status_alt == ["REP","REP","DIP"] : strategy["map_alt"] = "rep:altR:right"
	elif status_alt == ["REP","REP","GAP"] : strategy["map_alt"] = "rep:altR"
	elif status_alt == ["REP","REP","REP"] : strategy["map_alt"] = "rep:altR:rep"
	elif status_alt == ["REP","REP","OK"] : strategy["map_alt"] = "rep:altR:right"
	elif status_alt == ["REP","REP","NO"] : strategy["map_alt"] = "rep:altR"
	elif status_alt == ["REP","OK","DIP"] : strategy["map_alt"] = "rep:alt:right"
	elif status_alt == ["REP","OK","GAP"] : strategy["map_alt"] = "rep:alt"
	elif status_alt == ["REP","OK","REP"] : strategy["map_alt"] = "rep:alt:rep"
	elif status_alt == ["REP","OK","OK"] : strategy["map_alt"] = "rep:alt:right"
	elif status_alt == ["REP","OK","NO"] : strategy["map_alt"] = "rep:alt"
	elif status_alt == ["REP","NO","DIP"] : strategy["map_alt"] = "rep:gap:right"
	elif status_alt == ["REP","NO","OK"] : strategy["map_alt"] = "rep:gap:right"
	elif status_alt == ["OK","DIP","DIP"] : strategy["map_alt"] = "left:altF:right"
	elif status_alt == ["OK","DIP","GAP"] : strategy["map_alt"] = "left:altF"
	elif status_alt == ["OK","DIP","REP"] : strategy["map_alt"] = "left:altF:rep"
	elif status_alt == ["OK","DIP","OK"] : strategy["map_alt"] = "left:altF:right"
	elif status_alt == ["OK","DIP","NO"] : strategy["map_alt"] = "left:altF"
	elif status_alt == ["OK","GAP","DIP"] : strategy["map_alt"] = "left:gap:right"
	elif status_alt == ["OK","GAP","GAP"] : strategy["map_alt"] = "left"
	elif status_alt == ["OK","GAP","REP"] : strategy["map_alt"] = "left:gap:rep"
	elif status_alt == ["OK","GAP","OK"] : strategy["map_alt"] = "left:gap:right"
	elif status_alt == ["OK","GAP","NO"] : strategy["map_alt"] = "left"
	elif status_alt == ["OK","REP","DIP"] : strategy["map_alt"] = "left:altR:right"
	elif status_alt == ["OK","REP","GAP"] : strategy["map_alt"] = "left:altR"
	elif status_alt == ["OK","REP","REP"] : strategy["map_alt"] = "left:altR:rep"
	elif status_alt == ["OK","REP","OK"] : strategy["map_alt"] = "left:altR:right"
	elif status_alt == ["OK","REP","NO"] : strategy["map_alt"] = "left:altR"
	elif status_alt == ["OK","OK","DIP"] : strategy["map_alt"] = "left:alt:right"
	elif status_alt == ["OK","OK","GAP"] : strategy["map_alt"] = "left:alt"
	elif status_alt == ["OK","OK","REP"] : strategy["map_alt"] = "left:alt:rep"
	elif status_alt == ["OK","OK","OK"] : strategy["map_alt"] = "left:alt:right"
	elif status_alt == ["OK","OK","NO"] : strategy["map_alt"] = "left:alt"
	elif status_alt == ["OK","NO","DIP"] : strategy["map_alt"] = "left:gap:right"
	elif status_alt == ["OK","NO","GAP"] : strategy["map_alt"] = "left"
	elif status_alt == ["OK","NO","REP"] : strategy["map_alt"] = "left:gap:rep"
	elif status_alt == ["OK","NO","OK"] : strategy["map_alt"] = "left:gap:right"
	elif status_alt == ["OK","NO","NO"] : strategy["map_alt"] = "left"
	elif status_alt == ["NO","DIP","DIP"] : strategy["map_alt"] = "altF:right"
	elif status_alt == ["NO","DIP","GAP"] : strategy["map_alt"] = "altF"
	elif status_alt == ["NO","DIP","REP"] : strategy["map_alt"] = "altF:rep"
	elif status_alt == ["NO","DIP","OK"] : strategy["map_alt"] = "altF:right"
	elif status_alt == ["NO","DIP","NO"] : strategy["map_alt"] = "altF"
	elif status_alt == ["NO","GAP","DIP"] : strategy["map_alt"] = "right"
	elif status_alt == ["NO","GAP","OK"] : strategy["map_alt"] = "right"
	elif status_alt == ["NO","REP","DIP"] : strategy["map_alt"] = "altR:right"
	elif status_alt == ["NO","REP","GAP"] : strategy["map_alt"] = "altR"
	elif status_alt == ["NO","REP","REP"] : strategy["map_alt"] = "altR:rep"
	elif status_alt == ["NO","REP","OK"] : strategy["map_alt"] = "altR:right"
	elif status_alt == ["NO","REP","NO"] : strategy["map_alt"] = "altR"
	elif status_alt == ["NO","OK","DIP"] : strategy["map_alt"] = "alt:right"
	elif status_alt == ["NO","OK","GAP"] : strategy["map_alt"] = "alt"
	elif status_alt == ["NO","OK","REP"] : strategy["map_alt"] = "alt:rep"
	elif status_alt == ["NO","OK","OK"] : strategy["map_alt"] = "alt:right"
	elif status_alt == ["NO","OK","NO"] : strategy["map_alt"] = "alt"
	elif status_alt == ["NO","NO","DIP"] : strategy["map_alt"] = "right"
	elif status_alt == ["NO","NO","OK"] : strategy["map_alt"] = "right"
	else : strategy["map_alt"] = "NONE"
	#print >> sys.stderr, strategy
	return strategy


def agp_and_fasta_from_target_info( seq_id , db , new_name , dir ):
	#json.dump( db , sys.stderr , indent=4)
	targets_db = dict(db)
	# Write files
	file_prefix = dir + "/" + new_name
	agp_file , region_db , fasta , fasta_file , target_signal = ranges_to_agp( seq_id , db["data"] , new_name , file_prefix , gaplen = 100 )

	# Pet up output
	targets_db["id"] = new_name
	targets_db["sequence"] = fasta
	targets_db["sequence_file"] = fasta_file
	targets_db["agp_file"] = agp_file
	targets_db["signal"] = target_signal
	targets_db["agp_regions"] = region_db
	# targets_db["agp_regions"][new_name]
	#	 targets_db["agp_regions"][new_name][partNum]
	#	 	{
	#	 		"coords" : [ new_name , int(Obj_start) , int(Obj_End) ] ,
	#			"region_type" : [ UP | SEQ | DOWN | GAP ]
	#	 		"Compnt_Type" : [ "N" , "W" , "gap_substitute" ]
	#	 		"desc" : [ seq_id , CompntStart , CompntEnd ,  Orientation ] / int(gap_len)
	#	 	}

	return targets_db


def get_fasta_lengths_from_config(chunk_db) :
	length_db = {}
	for chr in chunk_db["sequences"] :
		length_db[chr] = chunk_db["sequences"][chr]["length"]
	return length_db


def make_sequences_and_signals( info , numeric_id ,  chunk_db , workdir ) :

	patch_db = {}
	seq_id = info["seq_id"]
	alt_seq_id = info["mate_id"]
	patching_strategy = info["patching_strategy"]
	gap_start =  int( info["gap_region"][0] )
	gap_stop = int( info["gap_region"][1] )
	# Make mapping target AGPs and FASTAs
	## Target regions, sequences and signals
	flanking_up_start , flanking_up_stop = info["flanking_upstream_region"]
	flanking_down_start , flanking_down_stop = info["flanking_downstream_region"]
	flanking_upstream_seq = info["fasta_sequences"]["flanking_upstream"]
	flanking_upstream_signal = info["flanking_upstream_content"]
	flanking_downstream_seq = info["fasta_sequences"]["flanking_downstream"]
	flanking_downstream_signal = info["flanking_downstream_content"]
	print >> sys.stderr , "##### Gap flanking regions: " + str( [ [ flanking_up_start , flanking_up_stop ] , [ flanking_down_start , flanking_down_stop ]  ] )

	if not info["gap_corr_region"][0] == "" :
		alt_start = int(info["gap_corr_region"][0])
		alt_stop = int(info["gap_corr_region"][1])
		alt_flanking_up_start , alt_flanking_up_stop = info["upstream_corr_region_region"]
		alt_flanking_down_start , alt_flanking_down_stop = info["downstream_corr_region_region"]
		alt_seq = info["fasta_sequences"]["corr_seq"]
		alt_signal = info["gap_corr_region_content"]
		alt_flanking_upstream_seq = info["fasta_sequences"]["corr_seq_upstream"]
		alt_flanking_upstream_signal = info["upstream_corr_region_content"]
		alt_flanking_downstream_seq = info["fasta_sequences"]["corr_seq_downstream"]
		alt_flanking_downstream_signal = info["downstream_corr_region_content"]
		print >> sys.stderr , "##### Gap flanking regions on alternative allele: " + str ([ [ alt_flanking_up_start , alt_flanking_up_stop ] , [ alt_flanking_down_start , alt_flanking_down_stop ]  ])
	else :
		alt_start = ""
		alt_stop = ""
		alt_flanking_up_start = ""
		alt_flanking_up_stop = ""
		alt_flanking_down_start = ""
		alt_flanking_down_stop = ""
		alt_seq = ""
		alt_signal = "0"
		alt_flanking_upstream_seq = ""
		alt_flanking_upstream_signal = "0"
		alt_flanking_downstream_seq = ""
		alt_flanking_downstream_signal = "0"
		print >> sys.stderr , "##### Gap region has no corresponent on alternative allele"

	if patching_strategy["map_gap"] in ( "flanking_rep:gap:right" , "flanking_left:gap:rep" , "flanking_left:gap:right", "flanking_rep:gap:rep" ) :
		target_data = [
			[seq_id , flanking_up_start , flanking_up_stop , flanking_upstream_seq , flanking_upstream_signal ] ,
			[seq_id , gap_start , gap_stop] ,
			[seq_id , flanking_down_start ,  flanking_down_stop, flanking_downstream_seq , flanking_downstream_signal ]
					]
	elif patching_strategy["map_gap"] == "flanking_left" :
		target_data = [
			[seq_id , flanking_up_start , flanking_up_stop , flanking_upstream_seq , flanking_upstream_signal ] ,
			[seq_id , gap_start , gap_stop] ,
			[seq_id , flanking_down_start ,  flanking_down_stop ]
					]
	elif patching_strategy["map_gap"] == "flanking_right" :
		target_data = [
			[seq_id , flanking_up_start , flanking_up_stop ] ,
			[seq_id , gap_start , gap_stop] ,
			[seq_id , flanking_down_start ,  flanking_down_stop, flanking_downstream_seq , flanking_downstream_signal ]
					]
	### Hybrid
	elif patching_strategy["map_gap"] in ("hybrid_left:alt:right" , "hybrid_left:alt:rep" , "hybrid_rep:alt:right" , "hybrid_rep:alt:rep") :
		target_data = [
			[seq_id , flanking_up_start , flanking_up_stop , flanking_upstream_seq , flanking_upstream_signal ] ,
			[alt_seq_id, alt_start , alt_stop , alt_seq , alt_signal ] ,
			[seq_id , flanking_down_start ,  flanking_down_stop, flanking_downstream_seq , flanking_downstream_signal ]
					]
	elif patching_strategy["map_gap"] in ("hybrid_left:alt" , "hybrid_rep:alt" ) :
		target_data = [
			[seq_id , flanking_up_start , flanking_up_stop , flanking_upstream_seq , flanking_upstream_signal ] ,
			[alt_seq_id, alt_start , alt_stop , alt_seq , alt_signal ] ,
			[seq_id , flanking_down_start ,  flanking_down_stop ]
					]
	elif patching_strategy["map_gap"] in ("hybrid_alt:right", "hybrid_alt:rep" ) :
		target_data = [
			[seq_id , flanking_up_start , flanking_up_stop ] ,
			[alt_seq_id, alt_start , alt_stop , alt_seq , alt_signal ] ,
			[seq_id , flanking_down_start ,  flanking_down_stop, flanking_downstream_seq , flanking_downstream_signal ]
					]
	else :
		target_data = ""
	patch_db["target_1"] = {}
	patch_db["target_1"]["data"] = target_data

	## Target alternative
	if patching_strategy["map_alt"] in ("alt" , "altF" , "altR") :
		target_data = [
			[alt_seq_id, alt_flanking_up_start , alt_flanking_up_stop ] ,
			[alt_seq_id, alt_start , alt_stop , alt_seq , alt_signal ] ,
			[alt_seq_id, alt_flanking_down_start , alt_flanking_down_stop ]
					]
	elif patching_strategy["map_alt"] in ("alt:right" , "altF:right" , "altR:right", "alt:rep" , "altF:rep" , "altR:rep") :
		target_data = [
			[alt_seq_id, alt_flanking_up_start , alt_flanking_up_stop  ] ,
			[alt_seq_id, alt_start , alt_stop , alt_seq , alt_signal ] ,
			[alt_seq_id, alt_flanking_down_start , alt_flanking_down_stop , alt_flanking_downstream_seq , alt_flanking_downstream_signal ]
					]
	elif patching_strategy["map_alt"] in ("left:alt" , "left:altF" , "left:altR", "rep:alt" , "rep:altF" , "rep:altR") :
		target_data = [
			[alt_seq_id, alt_flanking_up_start , alt_flanking_up_stop , alt_flanking_upstream_seq , alt_flanking_upstream_signal ] ,
			[alt_seq_id, alt_start , alt_stop , alt_seq , alt_signal ] ,
			[alt_seq_id, alt_flanking_down_start , alt_flanking_down_stop ]
					]
	elif patching_strategy["map_alt"] in ("rep:alt:right" , "rep:altF:right" , "rep:altR:right" , "left:alt:rep" , "left:altF:rep" , "left:altR:rep" , "left:alt:right" , "left:altF:right" , "left:altR:right" , "rep:alt:rep" , "rep:altF:rep" , "rep:altR:rep") :
		target_data = [
			[alt_seq_id, alt_flanking_up_start , alt_flanking_up_stop , alt_flanking_upstream_seq , alt_flanking_upstream_signal ] ,
			[alt_seq_id, alt_start , alt_stop , alt_seq , alt_signal ] ,
			[alt_seq_id, alt_flanking_down_start , alt_flanking_down_stop , alt_flanking_downstream_seq , alt_flanking_downstream_signal ]
					]
	elif patching_strategy["map_alt"] in ("left:gap:right" , "left:gap:rep" , "rep:gap:right") :
				target_data = [
			[alt_seq_id, alt_flanking_up_start , alt_flanking_up_stop , alt_flanking_upstream_seq , alt_flanking_upstream_signal ] ,
			[alt_seq_id, alt_start , alt_stop ] ,
			[alt_seq_id, alt_flanking_down_start , alt_flanking_down_stop , alt_flanking_downstream_seq , alt_flanking_downstream_signal ]
					]
	elif patching_strategy["map_alt"] == "right" :
		target_data = [
			[alt_seq_id, alt_flanking_up_start , alt_flanking_up_stop ] ,
			[alt_seq_id, alt_start , alt_stop ] ,
			[alt_seq_id, alt_flanking_down_start , alt_flanking_down_stop , alt_flanking_downstream_seq , alt_flanking_downstream_signal ]
					]
	elif patching_strategy["map_alt"] == "left" :
		target_data = [
			[alt_seq_id, alt_flanking_up_start , alt_flanking_up_stop , alt_flanking_upstream_seq , alt_flanking_upstream_signal ] ,
			[alt_seq_id, alt_start , alt_stop ] ,
			[alt_seq_id, alt_flanking_down_start , alt_flanking_down_stop ]
					]
	else :
		target_data = ""
	patch_db["target_2"] = {}
	patch_db["target_2"]["data"] = target_data

	# Write AGPSs and FASTAs
	## Target
	if not patch_db["target_1"] == "" :
		patch_db["target_1"] = agp_and_fasta_from_target_info( seq_id , patch_db["target_1"] , "Gap_" + str(numeric_id) + "_target_1" , workdir)
		# patch_db["patching_strategy"] = patching_strategy
		# patch_db["target_1"]
		#	patch_db["target_1"]["id"] = "Gap_" + str(numeric_id) + "_target_1"
		#	patch_db["target_1"]["signal"] = signal
		#	patch_db["target_1"]["sequence"] = fasta
		#	patch_db["target_1"]["sequence_file"] = fasta_file
		#	patch_db["target_1"]["agp_file"] = agp_file
		#	patch_db["target_1"]["agp_regions"] = region_db
		#	patch_db["target_1"]["agp_regions"]["target_1"]
		#		patch_db["target_1"]["agp_regions"]["target_1"][partNum]
		#		 	{
		#		 		"coords" : [ "target_1" , int(Obj_start) , int(Obj_End) ] ,
		#				"region_type" : [ UP | SEQ | DOWN | GAP ]
		#		 		"Compnt_Type" : [ "N" , "W" , "gap_substitute" ]
		#		 		"desc" : [ seq_id , CompntStart , CompntEnd ,  Orientation ] / int(gap_len)
		#		 	}
	else :
		patch_db["target_1"] = ""
	## Target alternative
	if not patch_db["target_2"] == "" :
		patch_db["target_2"] = agp_and_fasta_from_target_info( alt_seq_id , patch_db["target_2"] , "Gap_" + str(numeric_id) + "_target_2" , workdir)
		# patch_db["patching_strategy"] = patching_strategy
		# patch_db["target_2"]
		#	patch_db["target_2"]["id"] = "Gap_" + str(numeric_id) + "_target_2"
		#	patch_db["target_2"]["signal"] = signal
		#	patch_db["target_2"]["sequence"] = fasta
		#	patch_db["target_2"]["sequence_file"] = fasta_file
		#	patch_db["target_2"]["agp_file"] = agp_file
		#	patch_db["target_2"]["agp_regions"] = region_db
		#	patch_db["target_2"]["agp_regions"]["target_2"]
		#		patch_db["target_2"]["agp_regions"]["target_2"][partNum]
		#		 	{
		#		 		"coords" : [ "target_2" , int(Obj_start) , int(Obj_End) ] ,
		#				"region_type" : [ UP | SEQ | DOWN | GAP ]
		#		 		"Compnt_Type" : [ "N" , "W" , "gap_substitute" ]
		#		 		"desc" : [ seq_id , CompntStart , CompntEnd ,  Orientation ] / int(gap_len)
		#		 	}
	else :
		patch_db["target_2"] = ""

	return patch_db


def map_on_gap( align_db , align_db_file ,  target_1 , target_2 , signal_1 , signal_2 , chunk_db , cores , workdir , gap_db , unwanted_pairs_db , known_grouped_db_seqid ) :
	print >> sys.stderr , '### Gap side regions'
	if "target_1" not in align_db:
		align_db["target_1"] = {"map_file" : "TODO" , "map_info": "TODO"}
		align_db["target_2"] = {"map_file" : "TODO" , "map_info": "TODO"}
	## Target
	### Map
	if ( not "map_file" in align_db["target_1"] ) or ( not os.path.exists(align_db["target_1"]["map_file"]+".done") ) :
		print >> sys.stderr , '#### Mapping'
		align_db["target_1"]["map_file"] , unplaced_len = map_nucmer_unplaced_on_target(target_1, "map_unplaced_on_flaking.coords" , int(cores) , chunk_db, workdir)
		touch(align_db["target_1"]["map_file"]+".done")
		json.dump(align_db , gzip.open(align_db_file , "wb") , indent=4)
	else :
		print >> sys.stderr , '#### Mapping results already present, loading'
		unplaced_len = json.load(open(workdir + "/tmp.unplaced.len.json"))

	print >> sys.stderr , '### Alternative allele regions'
	## Target alternative
	### Map
	if ( not "map_file" in align_db["target_2"] ) or ( not os.path.exists(align_db["target_2"]["map_file"]+".done") ) :
		print >> sys.stderr , '#### Mapping'
		align_db["target_2"]["map_file"] , unplaced_len = map_nucmer_unplaced_on_target(target_2, "raw_map_unplaced_on_alternative.coords" , int(cores) , chunk_db , workdir)
		touch(align_db["target_2"]["map_file"]+".done")
		json.dump(align_db , gzip.open(align_db_file , "wb") , indent=4)
	else :
		print >> sys.stderr , '#### Mapping results already present, loading'
		unplaced_len = json.load(open(workdir + "/tmp.unplaced.len.json"))
	### Uniquify

	print >> sys.stderr , '#### Identify longest unique alignment paths'
	align_db["target_1"]["map_info"] , align_db["target_2"]["map_info"] = analize_unplaced_hits( align_db, signal_1 , signal_2 , unplaced_len, workdir , gap_db , unwanted_pairs_db , known_grouped_db_seqid )
	json.dump(align_db , gzip.open(align_db_file , "wb") , indent=4)

	return align_db


#def patch_or_fill( info , patching_strategy , chunk_db , workdir , status):
#
#	# check if still working to classify!
#
#	patch_db = info["patch"]
#
#
#	#	patching_strategy["map_gap"] = "value"
#	#	#	"NONE"  					|	unreliable	|		gap		|	unreliable	|
#	#	#	"flanking_left"  			|	reliable	|		gap		|	unreliable	|
#	#	#	"flanking_left:gap:rep" 	|	reliable	|		gap		|	repeat 		|
#	#	#	"flanking_left:gap:right" 	|	reliable	|		gap		|	reliable	|
#	#	#	"flanking_right" 			|	unreliable	|		gap		|	reliable 	|
#	#	#	"flanking_rep:gap:right" 	|	repeat		|		gap		|	unreliable	|
#	#	#	"flanking_rep:gap:rep" 		|	repeat		|		gap		|	repeat		|
#	#	#	"hybrid_left:alt:right" 	|	reliable	|	[patch alt]	|	reliable	|
#	#	#	"hybrid_left:alt" 			|	reliable	|	[patch alt]	|	unreliable	|
#	#	#	"hybrid_left:alt:rep" 		|	reliable	|	[patch alt]	|	repeat		|
#	#	#	"hybrid_alt:right" 			|	unreliable	|	[patch alt]	|	reliable	|
#	#	#	"hybrid_rep:alt:right" 		|	repeat		|	[patch alt]	|	reliable	|
#	#	#	"hybrid_rep:alt:rep" 		|	repeat		|	[patch alt]	|	repeat		|
#	#	#	"hybrid_rep:alt" 			|	repeat		|	[patch alt]	|	unreliable	|
#	#	#	"hybrid_alt:rep" 			|	unreliable	|	[patch alt]	|	repeat		|
#	#	patching_strategy["map_alt"] = "value"
#	#	#	"NONE" :					|	unreliable	|	unreliable	|	unreliable	|
#	#	#	"alt"						|	unreliable	|	reliable	|	unreliable	|
#	#	#	"altF"						|	unreliable	|	fill		|	unreliable	|
#	#	#	"altR"						|	unreliable	|	repeat		|	unreliable	|
#	#	#	"alt:right"					|	unreliable	|	reliable	|	reliable	|
#	#	#	"altF:right"				|	unreliable	|	fill		|	reliable	|
#	#	#	"altR:right"				|	unreliable	|	repeat		|	reliable	|
#	#	#	"rep:alt:right"				|	repeat		|	reliable	|	reliable	|
#	#	#	"rep:altF:right"			|	repeat		|	fill		|	reliable	|
#	#	#	"rep:altR:right"			|	repeat		|	repeat		|	reliable	|
#	#	#	"left:alt"					|	reliable	|	reliable	|	unreliable	|
#	#	#	"left:altF"					|	reliable	|	fill		|	unreliable	|
#	#	#	"left:altR"					|	reliable	|	repeat		|	unreliable	|
#	#	#	"left:alt:rep"				|	reliable	|	reliable	|	repeat		|
#	#	#	"left:altF:rep"				|	reliable	|	fill		|	repeat		|
#	#	#	"left:altR:rep"				|	reliable	|	repeat		|	repeat		|
#	#	#	"left:alt:right"			|	reliable	|	reliable	|	reliable	|
#	#	#	"left:altF:right"			|	reliable	|	fill		|	reliable	|
#	#	#	"left:altR:right"			|	reliable	|	repeat		|	reliable	|
#	#	#	"rep:alt:rep"				|	repeat		|	reliable	|	repeat		|
#	#	#	"rep:altF:rep"				|	repeat		|	fill		|	repeat		|
#	#	#	"rep:altR:rep"				|	repeat		|	repeat		|	repeat		|
#	#	#	"left:gap:right"			|	reliable	|	gap			|	reliable	|
#	#	#	"left:gap:rep"				|	reliable	|	gap			|	repeat		|
#	#	#	"rep:gap:right"				|	repeat		|	gap			|	reliable	|
#	#	#	"rep:alt"					|	repeat		|	reliable	|	unreliable	|
#	#	#	"rep:altF"					|	repeat		|	fill		|	unreliable	|
#	#	#	"rep:altR"					|	repeat		|	repeat		|	unreliable	|
#	#	#	"alt:rep"					|	unreliable	|	reliable	|	repeat		|
#	#	#	"altF:rep"					|	unreliable	|	fill		|	repeat		|
#	#	#	"altR:rep"					|	unreliable	|	repeat		|	repeat		|
#	#	#	"right"						|	unreliable	|	unreliable	|	reliable	|
#	#	#	"left"						|	reliable	|	unreliable	|	unreliable	|
#
#
#	# 	patch_db["target_x"]["map_info"]
#	#	#	patch_db["target_x"]["map_info"]["signal_length"] = signal_length
#	#	#	patch_db["target_x"]["map_info"]["ex_signal_length"] = ex_signal_length
#	#	#	patch_db["target_x"]["map_info"]["min_signal_alignment"] = min_signal_alignment
#	#	#	patch_db["target_x"]["map_info"]["min_ext_signal_alignment"] = min_ext_signal_alignment
#	#	#	patch_db["target_x"]["map_info"]["min_alignment_size_ratio"] = threshold
#	#	#	patch_db["target_x"]["map_info"]["classified_hits"]
#	#	#	#	patch_db["target_x"]["map_info"]["classified_hits"][Qid]
#	#	#	#	#	patch_db["target_x"]["map_info"]["classified_hits"][Qid][hit_type] with hit_type in [ "non-rep-alt" , "non-rep-ext" , "alt" , "ext"]
#	#	#	#	#	#	patch_db["target_x"]["map_info"]["classified_hits"][Qid][hit_type]["path"] = good_path
#	#	#	#	#	#	patch_db["target_x"]["map_info"]["classified_hits"][Qid][hit_type]["path_length"] = good_path_length
#	#	#	patch_db["target_x"]["map_info"]["best_paths"]
#	#	#	#	patch_db["target_x"]["map_info"]["best_paths"][hit_type] with hit_type in [ "non-rep-alt" , "non-rep-ext" , "alt" , "ext"]
#	#	#	#	#	patch_db["target_x"]["map_info"]["best_paths"][hit_type][good_path_length]=[ .. , Qid , ... ]
#	#	#	patch_db["target_x"]["map_info"]["best_match"]
#	#	#	#	patch_db["target_x"]["map_info"]["best_match"]["id"] = Query_id
#	#	#	#	patch_db["target_x"]["map_info"]["best_match"][hit_type] with hit_type in [ "non-rep-alt" , "non-rep-ext" , "alt" , "ext"] :
#	#	#	#	#	patch_db["target_x"]["map_info"]["best_match"]["hit_type"]["path"] = longest_path
#	#	#	#	#	patch_db["target_x"]["map_info"]["best_match"]["hit_type"]["path_length"] = longest_path_length
#
#		# Validate
#		if patch_db["status"] == "mapped" :
#			patch_db["filler_region"] = "NONE"
#			## select best patching sequence, otherwise return alt allele as fill up sequence if possible
#			if not patch_db["target_2"] == "" :
#				if "map_info" in patch_db["target_2"] :
#					if "best_match" in patch_db["target_2"]["map_info"] :
#						if "id" in patch_db["target_2"]["map_info"]["best_match"] :
#							patch_id = patch_db["target_2"]["map_info"]["best_match"]["id"][:-2]
#							Orientation = patch_db["target_2"]["map_info"]["best_match"]["id"][-1]
#							start = 1
#							stop = unplaced_len[patch_db["target_2"]["map_info"]["best_match"]["id"]]
#							patch_db["filler_region"] = [ patch_id , int(start) , int(stop) , Orientation]
#
#			if (patch_db["filler_region"] == "NONE") and (not patch_db["target_1"] == "") :
#				# target_2 gave no result, check target_1
#				if "map_info" in patch_db["target_1"] :
#					if "best_match" in patch_db["target_1"]["map_info"] :
#						if "id" in patch_db["target_1"]["map_info"]["best_match"] :
#							patch_id = patch_db["target_1"]["map_info"]["best_match"]["id"][:-2]
#							Orientation = patch_db["target_1"]["map_info"]["best_match"]["id"][-1]
#							start = 1
#							stop = unplaced_len[patch_db["target_1"]["map_info"]["best_match"]["id"]]
#							patch_db["filler_region"] = [ patch_id , int(start) , int(stop) , Orientation]
#
#			if patch_db["filler_region"] == "NONE" :
#				### Try to patch with alt hap sequence as no map info was good enough
#				if info["sequences_characteristics"]["gap_corr_region_content"] in [ "DIP" , "OK" ] :
#					patch_db["filler_region"] = [ info["mate_id"] , int(info["gap_corr_region"][0]) , int(info["gap_corr_region"][1]) , "+"]
#
#			patch_db["status"] = "DONE"
#		else :
#			# Integrity Control
#			sys.exit(1)
#	return patch_db
#
#
#def filter_maps( paf_hits , feat = "matches" , threshold=0.1 ) :
#	###### PAF format ######
#	# https://github.com/lh3/miniasm/blob/master/PAF.md
#	#
#	# Tab separated format
#	#
#	#   0 - Query sequence name
#	#   1 - Query sequence length
#	#   2 - Query start (0-based)
#	#   3 - Query end (0-based)
#	#   4 -	Relative strand: "+" or "-"
#	#   5 - Target sequence name
#	#   6 - Target sequence length
#	#   7 - Target start on original strand (0-based)
#	#   8 - Target end on original strand (0-based)
#	#   9 - Number of residue matches
#	#  10 - Alignment block length
#	#  11 - Mapping quality (0-255; 255 for missing)
#	#
#	########################
#	good_hits = {}
#	json.dump( paf_hits )
#	for path_id in paf_hits :
#		path_matches = float(paf_hits[path_id][9])
#		path_cov = float(paf_hits[path_id][10])
#		qlen = int(paf_hits[path_id][1])
#		tlen = int(paf_hits[path_id][6])
#		min_len = threshold*min( qlen , tlen )
#		#if ( feat == "cov" and path_cov > min_len) or ( feat == "matches" and path_matches > min_len) :
#		if ( path_cov > min_len) :
#			good_hits[path_id] = paf_hits[path_id]
#
#	return good_hits


def analize_unplaced_hits( alignment_info_db , signal_1_db , signal_2_db , query_seq_len , workdir , gap_db , unwanted_pairs_db , known_grouped_db_seqid, threshold = 0.2) :
	# unwanted_pairs_db
	# unwanted_pairs_db[unplaced_id] = [ chr_id , ... ]
	# unwanted_pairs_db[unplaced_id|+] = [ chr_id , ... ]
	# unwanted_pairs_db[unplaced_id|-] = [ chr_id , ... ]
	# unwanted_pairs_db[unplaced_id|.] = [ chr_id , ... ]

	# known_grouped_db_seqid
	# known_grouped_db_seqid[unplaced_id] = [ chr_id , ... ]
	# known_grouped_db_seqid[unplaced_id|+] = [ chr_id , ... ]
	# known_grouped_db_seqid[unplaced_id|-] = [ chr_id , ... ]
	# known_grouped_db_seqid[unplaced_id|.] = [ chr_id , ... ]

	#gap_db[gap_id] = {
	#	"dir" : gap_instance_dir ,
	#	"file": gap_instance_file ,
	#	"sequence_id" : chr ,
	#	"coordinates" : el["gap_region"] ,
	#	"flanking_upstream_region" : el["flanking_upstream_region"] ,
	#	"flanking_downstream_region" : el["flanking_downstream_region"] ,
	#	"corresponding_sequence_id" : el["mate_id"]  ,
	#	"corresponding_coordinates" : el["gap_corr_region"] ,
	#	"upstream_corr_region_region" : el["upstream_corr_region_region"] ,
	#	"downstream_corr_region_region" : el["downstream_corr_region_region"]
	#	"target_1_id" : el["patch"]["target_1"]["id"]
	#	"target_1_strategy" :  el["patching_strategy"]["map_gap"] # "hybrid_left:alt:right"
	#	"target_2_id" : el["patch"]["target_2"]["id"]
	#	"target_2_strategy" : el["patching_strategy"]["map_alt"] # "left:alt:right"
	#	}

	# Extract the chromosome of origin
	target_id_to_chr = {}
	for gap_id in gap_db.keys() :
		target_1_id = gap_db[gap_id]["target_1_id"]
		target_2_id = gap_db[gap_id]["target_2_id"]
		chr_id = gap_db[gap_id]["sequence_id"]
		target_id_to_chr[target_1_id] = chr_id
		target_id_to_chr[target_2_id] = chr_id

	pass_threshold = 100
	all_threshold = 1000
	map_info_1 = {}
	map_info_2 = {}
	classified_pass_hits = {}
	signal_db = {}
	signal_db.update(signal_1_db)
	signal_db.update(signal_2_db)
	raw_map_hits = read_nucmer_coords( alignment_info_db["target_1"]["map_file"] )
	raw_map_hits.update( read_nucmer_coords( alignment_info_db["target_2"]["map_file"] ) )
	# format
	# raw_map_hits[(Tid,Qid)] = [ ... , [ id , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ] , ...]
	# Qid = unplaced_id|orientation
	# Tid = target_sequence

	# Filter rejected wrong chromosome matching
	map_hits = {}
	for pair in raw_map_hits.keys() :
		target_sequence_id , unplaced_id = pair
		target_chr = target_id_to_chr[target_sequence_id]
		if ( unplaced_id in unwanted_pairs_db ) and ( target_chr in unwanted_pairs_db[unplaced_id] ):
			# There is constrain on unplaced_id placement and it has not be in that chr
			continue
		else :
			if unplaced_id in known_grouped_db_seqid :
				# There is constrain on unplaced_id placement in a chr
				if target_chr in known_grouped_db_seqid[unplaced_id] :
					map_hits[pair] = raw_map_hits[pair]
				else :
					# Constrain not met, alignment on the wrong chr
					continue
			else :
				# no given constrain on unplaced_id placement in a chr
				map_hits[pair] = raw_map_hits[pair]

	# Filter no pass hits
	print >> sys.stderr, "#### Parsing and filtering " + str(len(map_hits.keys())) + " pairwise alignment"
	counter = 0
	for key in map_hits :
		counter += 1
		if counter % 1000 == 0 :
			print >> sys.stderr, "##### Processed " + str(counter) + " pairs"

		target_name = key[0]
		query_id = key[1]

		#print >> sys.stderr, "##### " + str(key[0]) + " Vs. " + str(key[1])
		for hit in map_hits[key] :
			classified_hit = classify_hit( hit , signal_db[target_name] )
			# format:
			# classified_hit {
			#	 "filter" : "PASS"|"NOPASS"
			#	 "hit" : [ id , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ]
			#	}
			if classified_hit["filter"]	== "PASS" :
				if key not in classified_pass_hits :
					classified_pass_hits[key] = []
				classified_pass_hits[key].append(classified_hit["hit"])

	# Make graphs, extract longest path between each target ans query
	all_paths = {}
	pass_paths = []
	print >> sys.stderr, "#### Generating tiling paths"
	for entry in sorted(classified_pass_hits.keys()) :
		Tid , Qid = entry
		map_pass_graph = make_map_graph( classified_pass_hits[entry], 1000000 )

		longest_pass_path = nx.dag_longest_path(map_pass_graph, weight='align')
		longest_pass_graph = get_subgraph_from_path_tuples( map_pass_graph , longest_pass_path )
		longest_pass_graph_matches = longest_pass_graph.size(weight='match')

		#except :
		#	longest_pass_path = ""
		#	longest_pass_graph = ""
		#	longest_pass_graph_matches = 0

		if not longest_pass_graph_matches == 0 :
			map_all_graph = make_map_graph( map_hits[entry], 1000000 )
			longest_all_path = nx.dag_longest_path(map_all_graph, weight='align')
			longest_all_graph = get_subgraph_from_path_tuples( map_all_graph , longest_all_path )
			longest_all_graph_matches = longest_all_graph.size(weight='match')
			if entry not in all_paths :
				all_paths[entry] = {}
			all_paths[entry]["longest_path"] = longest_all_path
			all_paths[entry]["longest_graph"] = longest_all_graph
			all_paths[entry]["longest_graph_matches"] = longest_all_graph_matches
			pass_paths.append( [longest_pass_graph_matches , longest_all_graph_matches , entry[0] , entry[1] , longest_pass_path , longest_pass_graph] )

	# Sort paths my match length, assign the longest
	print >> sys.stderr, "#### Selecting longest tiling paths"
	pass_paths = sorted( pass_paths , key = lambda x: (x[0], x[1]) , reverse=True)
	pickle.dump(pass_paths , open( workdir + "/paths.pass.pkl" , 'w+') , pickle.HIGHEST_PROTOCOL )
	pickle.dump(all_paths  , open( workdir + "/paths.all.pkl"  , 'w+') , pickle.HIGHEST_PROTOCOL )
	while len(pass_paths) > 0 :
		longest_pass = pass_paths.pop(0)
		try :
			longest_pass_graph_matches , longest_all_graph_matches , Tid , Qid , longest_pass_path ,  longest_pass_graph  = longest_pass
		except :
			print >> sys.stderr, "[ERROR] line 3159 failed (longest_pass_graph_matches , longest_all_graph_matches , Tid , Qid , longest_pass_path ,  longest_pass_graph  = longest_pass)"
			print >> sys.stderr, longest_pass
			exit(1)
		else :
			longest_all_path = all_paths[(Tid , Qid)]["longest_path"]
			# Assign best hit
			if (int(longest_pass_graph_matches) > pass_threshold) and (int(longest_all_graph_matches) > all_threshold) :
				target_strategy = ""
				for id in gap_db :
					if gap_db[id]["target_1_id"] == Tid :
						target_strategy = gap_db[id]["target_1_strategy"]
						map_info_1[Tid] = {}
						map_info_1[Tid]["strategy"] = target_strategy
						map_info_1[Tid]["best_match"] = Qid
						map_info_1[Tid]["best_match_sequence_len"] = query_seq_len[Qid]
						map_info_1[Tid]["longest_pass_path"] = longest_pass_path
						map_info_1[Tid]["longest_pass_path_matches"] = longest_pass_graph_matches
						map_info_1[Tid]["longest_all_path"] = longest_all_path
						map_info_1[Tid]["longest_all_path_matches"] = longest_all_graph_matches
						pass_paths = remove_ids_from_pass(pass_paths , Tid , Qid )
					elif gap_db[id]["target_2_id"] == Tid :
						target_strategy = gap_db[id]["target_2_strategy"]
						map_info_2[Tid] = {}
						map_info_2[Tid]["strategy"] = target_strategy
						map_info_2[Tid]["best_match"] = Qid
						map_info_2[Tid]["best_match_sequence_len"] = query_seq_len[Qid]
						map_info_2[Tid]["longest_pass_path"] = longest_pass_path
						map_info_2[Tid]["longest_pass_path_matches"] = longest_pass_graph_matches
						map_info_2[Tid]["longest_all_path"] = longest_all_path
						map_info_2[Tid]["longest_all_path_matches"] = longest_all_graph_matches
						pass_paths = remove_ids_from_pass(pass_paths , Tid , Qid )

	return map_info_1 , map_info_2


def remove_ids_from_pass( path_list , target_id , query_id ) :
	new_list = []
	print >> sys.stderr, "target_id: " + target_id + " | query_id: " + query_id
	blacklist = [ target_id , query_id , target_id[:-1]+"+" , target_id[:-1]+"-" , query_id[:-1]+"+" , query_id[:-1]+"-" ]
	for element in path_list :
		if ( not element[2] in blacklist ) and (not element[3] in blacklist ) :
			#print >> sys.stderr , "Keep: " + element[2] + " -> " + element[3]
			new_list.append(element)
		#else :
		#	print >> sys.stderr , "Remove: " + element[2] + " -> " + element[3]

	return sorted( new_list , key = lambda x: (x[0], x[1]) , reverse=True)


def classify_hit( hit , target_signal ) :
	# format:
	# hits = [ id , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ]
	#
	# agp_regions[partNum] = {
	#	"coords" : [ target_id , int(Obj_start) , int(Obj_End) ] ,
	#	"region_type" : [ UP | SEQ | DOWN | GAP ]
	#	"Compnt_Type" : [ "N" , "W" , "gap_substitute" ]
	#	"desc" : [ seq_id , CompntStart , CompntEnd ,  Orientation ] / int(gap_len)
	#	}

	id , Tstart , Tstop , Qstart , Qstop , matches , hitLen = hit
	# Classify by position
	classified_element = {}
	classified_element["hit"] = hit

	hit_region_signal_type = analyze_signal( target_signal[ int(Tstart) : int(Tstop) ] )

	# Possible values:
	# DIP: Sequence >75% with coverage within diploid level value
	# GAP: Sequence made >75% of gap
	# REP: Sequence made >75% of repeats
	# HAP: Sequence made >75% of haploid regions
	# OK: dip + rep >75%
	# NO: Unreliable
	if hit_region_signal_type in [ "GAP" , "NO" , "REP" ] :
		classified_element["filter"] = "NOPASS"
	else :
		classified_element["filter"] = "PASS"

	# format:
	# classified_element = {
	#	 "filter" : "PASS"|"NOPASS"
	#	 "hit" : [ id , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ]
	#	}
	return classified_element


def check_patch_status( gap_dir ) :
	patch_descriptor_json_file_gz = gap_dir + "/gap_instance.json.gz"
	if not os.path.exists(patch_descriptor_json_file_gz) :
		return "TODO"
	else :
		gap_info = json.load( gzip.open(patch_descriptor_json_file_gz) )
		if not "patch" in gap_info :
			return "TODO"
		else :
			if not "status" in gap_info["patch"] :
				return "TOCHECK"
			else :
				return gap_info["patch"]["status"]


def whole_genome_dotplot( tIds , qIds , out_file_name_prefix, workdir, plot_folder, coords , reuse_dotplots = False , skip_dotplots = False, minimum_size = 3000 , minimum_identity = 85) :
	# tIds and qIds = string(chr01_id,chr02_id,...)
	# coords_file , "tID\ttLen\ttStart\ttStop\tqID\tqLen\tqStart\tqStop\tidentity\tmatch"
	output_dir = workdir + "/" + plot_folder
	output_file = output_dir + "/" + out_file_name_prefix + ".dotplot"
	output_file_relative = plot_folder + "/" + out_file_name_prefix + ".dotplot"
	mkdir(output_dir)

	## all inter-chromosomal plots
	all_dotplots = {}
	if not skip_dotplots :
		for tid in tIds.split(",") :
			all_dotplots[tid] = {}
			for qid in qIds.split(",") :
				file_prefix	= qid + ".on." + tid + ".dotplot"
				if not reuse_dotplots :
					file_full_path = output_dir + "/" + file_prefix
					SimpleDotplot_script = scriptDirectory + "/SimpleDotplot.Rmd"

					log_connection = open( output_dir + "/." + file_prefix + ".log" , 'w')
					err_connection = open( output_dir + "/." + file_prefix + ".err" , 'w')
					command = "Rscript -e 'library(rmarkdown) ; rmarkdown::render(\"" + SimpleDotplot_script + "\" , knit_root_dir = \"" + os.path.realpath(".") + "\" , output_file = \"" + file_prefix + ".html\" , output_dir = \"" + os.path.realpath(output_dir) + "\" , params=list(coords = \"" + coords + "\" , filename = \"" + file_full_path + "\" , queryID = \"" + qid + "\" , refID = \"" + tid + "\" , identity = \"" + str(minimum_identity) + "\" , match = \"" + str(minimum_size) + "\"))'"
					print >> sys.stderr, "#### Running command: " + command
					reportProcess = subprocess.Popen( command , shell=True , stdout=log_connection , stderr=err_connection )
					output, error = reportProcess.communicate()
					log_connection.close()
					err_connection.close()

				all_dotplots[tid][qid] = {}
				all_dotplots[tid][qid]["html"] = file_prefix + ".html"
				all_dotplots[tid][qid]["pdf"] = file_prefix + ".pdf"
				all_dotplots[tid][qid]["png"] = file_prefix + ".png"

	if not reuse_dotplots :
		original_script = scriptDirectory + "/WholeGenomeDotplot.Rmd"
		script = output_file + ".Rmd"
		script_file = open(script, 'a+')
		shutil.copyfileobj(open(original_script), script_file)
		# Add lines to the script to build an index of all chr-vs-chr dotplots
		print >> script_file , ""
		print >> script_file , "# Chromosome vs Chromosome Dotplot index"
		print >> script_file , ""
		print >> script_file , "Dotplots of by chromosome sequence"
		print >> script_file , ""
		print >> script_file , "Target sequence:"
		print >> script_file , ""
		print >> script_file , "## Column {.tabset}"
		for tid in sorted(all_dotplots.keys()):
			box_height = 42 + ( 22 * len(all_dotplots[tid].keys()) )
			print >> script_file , "### " + tid + " {data-height=" + str(box_height) + "}"
			print >> script_file , ""
			for query_id in sorted(all_dotplots[tid].keys()) :
				html_file_path = all_dotplots[tid][query_id]["html"]
				pdf_file_path = all_dotplots[tid][query_id]["pdf"]
				png_file_path = all_dotplots[tid][query_id]["png"]
				print >> script_file , "* Query: " + str(query_id) + " >>> [html](" + html_file_path + ") | [pdf](" + pdf_file_path + ") | [png](" + png_file_path + ")"
				print >> script_file , ""
		script_file.close()

		log_connection = open( output_dir + "/." + out_file_name_prefix + ".dotplot" + ".log" , 'w')
		err_connection = open( output_dir + "/." + out_file_name_prefix + ".dotplot" + ".err" , 'w')
		command = "Rscript -e 'library(rmarkdown) ; rmarkdown::render(\"" + os.path.realpath(script) + "\" , knit_root_dir = \"" + os.path.realpath(".") + "\" , output_file = \"" + out_file_name_prefix + ".dotplot.html\" , output_dir = \"" + os.path.realpath(output_dir) + "\" , params=list(coords = \"" + coords + "\" , filename = \"" + output_file + "\" , identity = \"" + str(minimum_identity) + "\" , match = \"" + str(minimum_size) + "\"))'"
		print >> sys.stderr, "#### Running command: " + command
		reportProcess = subprocess.Popen( command , shell=True , stdout=log_connection , stderr=err_connection )
		output, error = reportProcess.communicate()
		log_connection.close()
		err_connection.close()
	outfiles = {}
	outfiles["html"] = output_file_relative + ".html"
	outfiles["pdf"] = output_file_relative + ".pdf"
	outfiles["png"] = output_file_relative + ".png"


	return outfiles , all_dotplots


def make_coords_table( delta_prefix , workdir , command_line) :
	coords_table_file_name = workdir + "/" + delta_prefix + ".coords.txt"
	input_delta = workdir + "/" + delta_prefix + ".delta"
	extract_coords_process = command_line + " -Tl " + input_delta + " | tail -n +5 "
	temp_file = open( coords_table_file_name + ".tmp" , 'w' )
	coordsProcess = subprocess.Popen(extract_coords_process , shell=True, stdout=temp_file)
	output, error = coordsProcess.communicate()
	coords_file = open( coords_table_file_name , 'w' )
	#awk -F \"\t\" 'BEGIN {OFS=\"\t\" ; print \"tID\",\"tLen\",\"tStart\",\"tStop\",\"qID\",\"qLen\",\"qStart\",\"qStop\",\"identity\",\"match\"} {print $10,$8,$1,$2,$11,$9,$3,$4,$7,$5}' "
	print >> coords_file , "tID\ttLen\ttStart\ttStop\tqID\tqLen\tqStart\tqStop\tidentity\tmatch"
	for line in open(coords_table_file_name + ".tmp" , 'r') :
		coords_columns = line.rstrip().split("\t")
		print >> coords_file , "\t".join([ coords_columns[9] ,coords_columns[7] ,coords_columns[0] ,coords_columns[1] ,coords_columns[10] ,coords_columns[8] ,coords_columns[2] ,coords_columns[3] ,coords_columns[6] ,coords_columns[4] ] )
	coords_file.close()
	os.remove(coords_table_file_name + ".tmp")
	return coords_table_file_name


def make_chunk_coords_table_from_self_paf(paf_file, workdir, filename) :
	coords_file = open( filename , 'w' )
	print >> coords_file , "tID\ttLen\ttStart\ttStop\tqID\tqLen\tqStart\tqStop\tidentity\tmatch"
	for line in open(paf_file) :
		coords_columns = line.rstrip().split("\t")[0:11]
		# coords_columns = [Qid , Qlen , Qstart , Qstop , strand , Tid , Tlen , Tstart , Tstop , matches , hitLen]
		#					0     1      2        3       4        5     6      7        8       9         10
		# Extract for Qid_name and correct by orientation
		qID , qStart , qStop = coords_columns[0].split("#")
		pos_start = int(qStart)
		pos_top = int(qStop)
		q_fragment_len = int(coords_columns[1])
		q_fragment_start = int(coords_columns[2])
		q_fragment_stop = int(coords_columns[3])
		delta_start = q_fragment_start
		delta_stop = q_fragment_len - q_fragment_stop
		q_start = pos_start + max( 0, delta_start)
		q_stop = pos_top - max( 0, delta_stop)
		if coords_columns[4] == "+" :
			qStart = min( q_start , q_stop)
			qStop = max( q_start , q_stop)
		elif coords_columns[4] == "-" :
			qStart = max( q_start , q_stop)
			qStop = min( q_start , q_stop)
		else :
			print >> sys.stdout, "[ERROR] Unexpected file format for paf file " + paf_file
			print >> sys.stderr, "[ERROR] Unexpected file format for paf file " + paf_file
			sys.exit(3)
		identity = float(coords_columns[9]) / float(coords_columns[10])
		#                                  tID                tLen               tStart             tStop              qID  qLen               qStart       qStop        identity        match
		print >> coords_file , "\t".join([ coords_columns[5] ,coords_columns[6] ,coords_columns[7] ,coords_columns[8] ,qID ,coords_columns[6] ,str(qStart) ,str(qStop) , str(identity) , coords_columns[10] ] )
	coords_file.close()
	return filename


def make_pair_html_report(coords, coords_self, workdir, output_dir, queryID, refID, structure , legacy , markers , dup_markers , hap1ID , hap2ID , hap1Len , hap2Len, counts_hap1 ="diploid_gene_count_trace.hap1.txt", counts_hap2 ="diploid_gene_count_trace.hap2.txt", min_align ="3000", similarity ="90", ratio="0.33") :
	# 	All scripts accept the same input, make use only of partial info
	if hap1ID == "" :
		hap1ID = refID
		hap2ID = queryID

	if not structure == "":
		if not legacy == "":
			if not markers == "":
				if not coords_self == "" :
					script=scriptDirectory + "/ChrBoard.html.Rmd"
				else :
					script=scriptDirectory + "/ChrBoard.html.noself.Rmd"
			else :
				if not coords_self == "" :
					script=scriptDirectory + "/ChrBoard.html.nomarkers.Rmd"
				else :
					script=scriptDirectory + "/ChrBoard.html.nomarkers.noself.Rmd"
		else :
			if not markers == "":
				if not coords_self == "" :
					script=scriptDirectory + "/ChrBoard.html.nolegacy.Rmd"
				else :
					script=scriptDirectory + "/ChrBoard.html.nolegacy.noself.Rmd"
			else :
				if not coords_self == "" :
					script=scriptDirectory + "/ChrBoard.html.nolegacy.nomarkers.Rmd"
				else :
					script=scriptDirectory + "/ChrBoard.html.nolegacy.nomarkers.noself.Rmd"
	else :
		if not legacy == "":
			if not markers == "":
				if not coords_self == "" :
					script=scriptDirectory + "/ChrBoard.html.nostructure.Rmd"
				else :
					script=scriptDirectory + "/ChrBoard.html.nostructure.noself.Rmd"
			else :
				if not coords_self == "" :
					script=scriptDirectory + "/ChrBoard.html.nostructure.nomarkers.Rmd"
				else :
					script=scriptDirectory + "/ChrBoard.html.nostructure.nomarkers.noself.Rmd"
		else :
			if not markers == "":
				if not coords_self == "" :
					script=scriptDirectory + "/ChrBoard.html.nostructure.nolegacy.Rmd"
				else :
					script=scriptDirectory + "/ChrBoard.html.nostructure.nolegacy.noself.Rmd"
			else :
				if not coords_self == "" :
					script=scriptDirectory + "/ChrBoard.html.nostructure.nolegacy.nomarkers.Rmd"
				else :
					script=scriptDirectory + "/ChrBoard.html.nostructure.nolegacy.nomarkers.noself.Rmd"

	report_file = queryID + ".on." + refID + ".report.html"
	log_connection = open( output_dir + "/." + report_file + ".log" , 'w')
	err_connection = open( output_dir + "/." + report_file + ".err", 'w')
	command = "Rscript -e 'library(rmarkdown) ; rmarkdown::render(\"" + script + "\" , knit_root_dir = \"" + workdir + "\" , output_file = \"" + report_file + "\" , output_dir = \"" + output_dir + "\" , params=list(coords = \"" + coords + "\" , coords_self = \"" + coords_self + "\" , counts_hap1 = \"" + counts_hap1 + "\" , counts_hap2 = \"" + counts_hap2 + "\" , min_align = \"" + str(min_align) + "\" , similarity = \"" + str(similarity) + "\" , queryID = \"" + queryID + "\" , refID = \"" + refID + "\" , hap1ID = \"" + hap1ID + "\" , hap2ID = \"" + hap2ID + "\" , hap1Len = \"" + str(hap1Len) + "\" , hap2Len = \"" + str(hap2Len) + "\" , ratio= \"" + str(ratio) + "\" , structure = \"" + structure + "\" , legacy = \"" + legacy + "\" , markers = \"" + markers + "\" , dup_markers = \"" + dup_markers + "\" ))'"
	print >> sys.stderr , "#### Running command: " + command
	reportProcess = subprocess.Popen( command , shell=True , stdout=log_connection , stderr=err_connection )
	output, error = reportProcess.communicate()
	log_connection.close()
	err_connection.close()
	relative_dir = os.path.basename(output_dir)
	return relative_dir + "/" + report_file


def make_pair_pdf_report(coords, coords_self, workdir, output_dir, queryID, refID, structure , legacy , markers , dup_markers , hap1ID , hap2ID , counts_hap1 , min_align ="3000", similarity ="90", ratio="0.33") :
	# Rscript --vanilla support_scripts/ChrBoard.pdf.R -c coords.txt -s self_map.txt -g diploid_gene_count_trace.hap1.txt -m 3000 -i 90 -q "VITMroTrayshed_v2.0.hap2.chr01" -t "VITMroTrayshed_v2.0.hap1.chr01" -r 0.33 -d ~/Desktop -o test.pdf -a "structure.txt" -b "markers.bed" -e "dup_markers.bed"
	# 	All scripts accept the same input, make use only of partial info
	report_file = queryID + ".on." + refID + ".report.pdf"
	log_connection = open( output_dir + "/." + queryID + ".on." + refID + ".report.pdf.log" , 'w')
	err_connection = open( output_dir + "/." + queryID + ".on." + refID + ".report.pdf.err" , 'w' )
	if not structure == "" :
		a = workdir + "/" + structure
		if not legacy == "":
			l = workdir + "/" + legacy
			if not markers == "":
				b = workdir + "/" + markers
				e = workdir + "/" + dup_markers
				if not coords_self == "" :
					s = workdir + "/" + coords_self
					script=scriptDirectory + "/ChrBoard.pdf.R"
				else :
					s = "0"
					script=scriptDirectory + "/ChrBoard.pdf.noself.R"
			else :
				b = "0"
				e = "0"
				if not coords_self == "" :
					s = workdir + "/" + coords_self
					script=scriptDirectory + "/ChrBoard.pdf.nomarkers.R"
				else :
					s = "0"
					script=scriptDirectory + "/ChrBoard.pdf.nomarkers.noself.R"
		else :
			l = "0"
			if not markers == "":
				b = workdir + "/" + markers
				e = workdir + "/" + dup_markers
				if not coords_self == "" :
					s = workdir + "/" + coords_self
					script=scriptDirectory + "/ChrBoard.pdf.nolegacy.R"
				else :
					s = "0"
					script=scriptDirectory + "/ChrBoard.pdf.nolegacy.noself.R"
			else :
				b = "0"
				e = "0"
				if not coords_self == "" :
					s = workdir + "/" + coords_self
					script=scriptDirectory + "/ChrBoard.pdf.nolegacy.nomarkers.R"
				else :
					s = "0"
					script=scriptDirectory + "/ChrBoard.pdf.nolegacy.nomarkers.noself.R"
	else :
		a = "0"
		if not legacy == "":
			l = workdir + "/" + legacy
			if not markers == "":
				b = workdir + "/" + markers
				e = workdir + "/" + dup_markers
				if not coords_self == "" :
					s = workdir + "/" + coords_self
					script=scriptDirectory + "/ChrBoard.pdf.nostructure.R"
				else :
					s = "0"
					script=scriptDirectory + "/ChrBoard.pdf.nostructure.noself.R"
			else :
				b = "0"
				e = "0"
				if not coords_self == "" :
					s = workdir + "/" + coords_self
					script=scriptDirectory + "/ChrBoard.pdf.nostructure.nomarkers.R"
				else :
					s = "0"
					script=scriptDirectory + "/ChrBoard.pdf.nostructure.nomarkers.noself.R"
		else :
			l = "0"
			if not markers == "":
				b = workdir + "/" + markers
				e = workdir + "/" + dup_markers
				if not coords_self == "" :
					s = workdir + "/" + coords_self
					script=scriptDirectory + "/ChrBoard.pdf.nostructure.nolegacy.R"
				else :
					s = "0"
					script=scriptDirectory + "/ChrBoard.pdf.nostructure.nolegacy.noself.R"
			else :
				b = "0"
				e = "0"
				if not coords_self == "" :
					s = workdir + "/" + coords_self
					script=scriptDirectory + "/ChrBoard.pdf.nostructure.nolegacy.nomarkers.R"
				else :
					s = "0"
					script=scriptDirectory + "/ChrBoard.pdf.nostructure.nolegacy.nomarkers.noself.R"

	command = "Rscript --vanilla " + script + " -d " + output_dir + " -o " + report_file + " -c " + workdir + "/" + coords + " -s " + s + " -g " + workdir + "/" + counts_hap1 + " -m " + str(min_align) + " -i " + str(similarity) + " -q " + queryID + " -t " + refID + " -r " + str(ratio) + " -a " + a + " -l " + l + " -b " + b + " -e " + e
	print >> sys.stderr , "#### Running command: " + command
	reportProcess = subprocess.Popen( command , shell=True , stdout=log_connection , stderr=err_connection )
	output, error = reportProcess.communicate()
	log_connection.close()
	err_connection.close()
	relative_dir = os.path.basename(output_dir)
	return relative_dir + "/" + report_file


def make_no_genes_html_report(coords, coords_self, workdir, output_dir, queryID, refID, structure = "" , legacy = "" , markers = "" , dup_markers = "" , min_align ="3000", similarity ="90") :
	if not structure == "":
		if not legacy == "":
			if not markers == "":
				script=scriptDirectory + "/ChrReport_nogene.html.Rmd"
			else :
				script=scriptDirectory + "/ChrReport_nogene.html.nomarkers.Rmd"
		else :
			if not markers == "":
				script=scriptDirectory + "/ChrReport_nogene.html.nolegacy.Rmd"
			else :
				script=scriptDirectory + "/ChrReport_nogene.html.nolegacy.nomarkers.Rmd"
	else :
		if not legacy == "":
			if not markers == "":
				script=scriptDirectory + "/ChrReport_nogene.html.nostructure.Rmd"
			else :
				script=scriptDirectory + "/ChrReport_nogene.html.nostructure.nomarkers.Rmd"
		else :
			if not markers == "":
				script=scriptDirectory + "/ChrReport_nogene.html.nostructure.nolegacy.Rmd"
			else :
				script=scriptDirectory + "/ChrReport_nogene.html.nostructure.nolegacy.nomarkers.Rmd"

	report_file = queryID + ".on." + refID + ".report.html"
	log_connection = open( output_dir + "/." + report_file + ".log" , 'w')
	err_connection = open( output_dir + "/." + report_file + ".err", 'w')
	command = "Rscript -e 'library(rmarkdown) ; rmarkdown::render(\"" + script + "\" , knit_root_dir = \"" + workdir + "\" , output_file = \"" + report_file + "\" , output_dir = \"" + output_dir + "\" , params=list(coords = \"" + coords + "\" , coords_self = \"" + coords_self + "\" , min_align = \"" + str(min_align) + "\" , similarity = \"" + str(similarity) + "\" , queryID = \"" + queryID + "\" , refID = \"" + refID + "\" , structure = \"" + structure + "\" , legacy = \"" + legacy + "\" , markers = \"" + markers + "\" , dup_markers = \"" + dup_markers + "\" ))'"
	print >> sys.stderr , "#### Running command: " + command
	reportProcess = subprocess.Popen( command , shell=True , stdout=log_connection , stderr=err_connection )
	output, error = reportProcess.communicate()
	log_connection.close()
	err_connection.close()
	relative_dir = os.path.basename(output_dir)
	return relative_dir + "/" + report_file


def make_no_genes_pdf_report(coords, coords_self, workdir, output_dir, queryID, refID, structure = "" , legacy = "" , markers = "" , dup_markers = "" , min_align ="3000", similarity ="90") :
	# Rscript --vanilla support_scripts/ChrReport_nogene.pdf.R -c coords.txt -s self_map.txt -m 3000 -i 90 -q "VITMroTrayshed_v2.0.hap2.chr01" -t "VITMroTrayshed_v2.0.hap1.chr01" -r 0.33 -d ~/Desktop -o test.pdf
	if not structure =="" :
		a = workdir + "/" + structure
		if not legacy == "":
			l = workdir + "/" + legacy
			if not markers == "":
				b = workdir + "/" + markers
				e = workdir + "/" + dup_markers
				script=scriptDirectory + "/ChrReport_nogene.pdf.R"
			else :
				b = "0"
				e = "0"
				script=scriptDirectory + "/ChrReport_nogene.pdf.nomarkers.R"
		else :
			l = "0"
			if not markers == "":
				b = workdir + "/" + markers
				e = workdir + "/" + dup_markers
				script=scriptDirectory + "/ChrReport_nogene.pdf.nolegacy.R"
			else :
				b = "0"
				e = "0"
				script=scriptDirectory + "/ChrReport_nogene.pdf.nolegacy.nomarkers.R"
	else :
		a = 0
		if not legacy == "":
			l = workdir + "/" + legacy
			if not markers == "":
				b = workdir + "/" + markers
				e = workdir + "/" + dup_markers
				script=scriptDirectory + "/ChrReport_nogene.pdf.nostructure.R"
			else :
				b = "0"
				e = "0"
				script=scriptDirectory + "/ChrReport_nogene.pdf.nostructure.nomarkers.R"
		else :
			l = 0
			if not markers == "":
				b = workdir + "/" + markers
				e = workdir + "/" + dup_markers
				script=scriptDirectory + "/ChrReport_nogene.pdf.nostructure.nolegacy.R"
			else :
				b = "0"
				e = "0"
				script=scriptDirectory + "/ChrReport_nogene.pdf.nostructure.nolegacy.nomarkers.R"

	report_file = queryID + ".on." + refID + ".report.pdf"
	log_connection = open(output_dir + "/." + report_file + ".log" , 'w')
	err_connection = open(output_dir + "/." + report_file + ".err" , 'w' )
	command = "Rscript --vanilla " + script + " -d " + output_dir + " -o " + report_file + "-c " + coords + " -s " + coords_self + " -m " + str(min_align) + " -i " + str(similarity) + " -q " + queryID + " -t " + refID + " -a " + a + " -l " + l + " -b " + b + " -e " + e
	print >> sys.stderr , "#### Running command: " + command
	reportProcess = subprocess.Popen( command , shell=True , stdout=log_connection , stderr=err_connection )
	output, error = reportProcess.communicate()
	log_connection.close()
	err_connection.close()
	relative_dir = os.path.basename(output_dir)
	return relative_dir + "/" + report_file


def read_block( file_name ) :
	# block file format:
	#	#	>OutSeqID_1
	#	#	SeqID_1	Start	Stop	Strand
	#	#	SeqID_2	Start	Stop	Strand
	#	#	...
	#	#	>OutSeqID_2
	#	#	SeqID_N	Start	Stop	Strand
	#	#	SeqID_M	Start	Stop	Strand
	#	#	...

	block_db = {}
	sequence_db = {}
	sequence_id = ""
	block_id = 0
	for line in open(file_name) :
		if line.rstrip() == "" or line[0] == "#" :
			continue
		elif line[0] == ">" :
			# New sequence
			if not sequence_id == "" :
				block_db[sequence_id] = sequence_db
			sequence_db = {}
			sequence_id = str(line.rstrip()[1:].split()[0])
			block_id = 0
		else :
			block_id += 1
			try :
				seqID , start , stop , strand = line.rstrip().split()
			except :
				print >> sys.stdout , "[Error] Parsing " + file_name + ": block line not in the proper format. See standard error for more details"
				print >> sys.stderr , "[Error] Parsing " + file_name + ": block line not in the proper format"
				print >> sys.stderr , "[Error] Expected 4 fields: [ SeqID , Start , Stop , Strand ] , " + str(len(line.rstrip().split())) + " found"
				print >> sys.stderr , "[Error] Conflicting line: "
				print >> sys.stderr , line
				print >> sys.stderr , "[Error] Splitting result: "
				print >> sys.stderr , str( line.rstrip().split() )
				sys.exit(1)
			sequence_db[block_id] = [seqID , int(start) , int(stop) , strand]
	block_db[sequence_id] = sequence_db
	return block_db


def dodge_overlaps(block_dict , fasta_db , fasta_len_db , gap_size , spacer , cores, mapper , annotation , tempdir , paths , add_unplaced ) :
	new_agp = {}
	gene_position = feature_ranges( annotation , "gene" )
	harmed_loci = []
	for chr in sorted(block_dict.keys()) :
		print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Sequence: ' + str(chr)
		print >> sys.stderr , '## Sequence: ' + str(chr)
		left_block = {}
		new_agp[chr] = {}
		# block_dict[chr][block_id] = [seqID , int(start) , int(stop) , strand]
		for block_id in sorted(block_dict[chr].keys()) :
			right_block = { "region_given" : block_dict[chr][block_id] }
			right_block["region_corrected"] , right_block["seq"] , right_block["annot_on_fasta"] = get_block_extremities(fasta_db, fasta_len_db, right_block["region_given"][0], int(right_block["region_given"][1]), int(right_block["region_given"][2]) , right_block["region_given"][3] , gene_position)
			if not left_block == {} :
				# Compare left block with right block
				# map with blat or nucmer
				print >> sys.stderr , '### Mapping ' + str(left_block["region_corrected"]) + " on " + str(right_block["region_corrected"])
				mappings = map_regions( left_block["seq"] , right_block["seq"] , mapper , cores , tempdir , paths )
				# Returns a coords dict:
				# mappings[(left_block,right_block)].append([ id , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ]
				# Parse mappings
				print >> sys.stderr , '### Searching extremities overlap boundaries'
				overlap_region = find_extremity_overlap( mappings , left_block , right_block)
				# Correct coordinates to respect genes
				print >> sys.stderr , '### Update regions upon overlap boundaries'
				left_block , right_block , left_harmed_loci , refining_status = refine_regions( left_block , right_block , overlap_region )
				# Add left block to the
				new_agp[chr] = add_block_to_agp(new_agp , chr , left_block , gap_size , spacer , refining_status )
				harmed_loci += left_harmed_loci

			if "region_trimmed" in right_block:
				# update right_block region_corrected to reflect trimming if needed
				print >> sys.stderr , "#### Update region " + str(right_block["region_given"]) + " >>> trimmed " + str(right_block["region_corrected"])
				right_block["region_corrected"] , right_block["seq"] , right_block["annot_on_fasta"] = get_block_extremities(fasta_db, fasta_len_db, right_block["region_trimmed"][0], int(right_block["region_trimmed"][1]), int(right_block["region_trimmed"][2]) , right_block["region_trimmed"][3] , gene_position)
				del(right_block["region_trimmed"])

			left_block = dict(right_block)

		new_agp[chr] = add_block_to_agp(new_agp , chr , left_block , gap_size , spacer , "last" )

	if add_unplaced :
		print >> sys.stderr , '### Adding unplaced sequences to output'
		used_sequences = {}
		# new_agp[seq_id][int(start)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation ]
		for seq_id in new_agp.keys() :
			for start in new_agp[seq_id].keys() :
				if new_agp[seq_id][int(start)][4] == "W" :
					used_sequences[ new_agp[seq_id][int(start)][5] ] = "1"
		print >> sys.stderr , '#### Used sequences: ' + str( len(used_sequences.keys()) ) + "/" + str( len(fasta_len_db.keys()) )

		for seq_id in sorted(fasta_len_db.keys()) :
			if seq_id not in used_sequences :
				seq_length = fasta_len_db(seq_id)
				new_agp[seq_id][1] = [ seq_id , 1 , seq_length , 1 , "W" , seq_id , 1 , seq_length ,  "+" ]

	return new_agp , harmed_loci


def refine_regions( left_db , right_db , overlapping_regions_pair ) :
	#print >> sys.stderr , "## Refining overlapping regions: [ Left = " + str(left_db["region_given"]) + " ] | [ Right region = " + str(right_db["region_given"]) + " ]"
	new_left_db = left_db
	new_right_db = right_db
	left_matched_loci_names = []
	right_matched_loci = []

	if overlapping_regions_pair == {} :
		new_left_db["region_trimmed"] = left_db["region_corrected"]
		new_right_db["region_trimmed"] = right_db["region_corrected"]
		refining_status = "disjointed"
	else :
		left_start , left_stop , right_start , right_stop = overlapping_regions_pair
		# Try trim the right sequence up to right_stop if no gene overlaps the break
		# Otherwise >> edit the trim regions on the left to keep the gene present on the right intact

		# Check if right stop is on top of a gene or left start
		for gene in right_db["annot_on_fasta"] :
			#gene = [ gene_start , gene_stop , gene_id ]
			if int(gene[0]) < int(right_stop) < int(gene[1]) :
				print >> sys.stderr, "#### Gene in the right sequence in the region overlap: " + str(gene)
				right_matched_loci.append(gene)

		if right_matched_loci == [] :
			# no genes is harmed , trim the right region at the right stop overlap
			shift = 0
		else :
			print >> sys.stderr, "###### Shifting region selection to accommodate gene(s) in the right sequence"
			# use the leftmost gene start to define the shift
			min_start = 10000000000000
			for gene in right_matched_loci :
				if int(gene[0]) < min_start :
					min_start = int(gene[0])

			# add some space upstream, 1 codon or 0
			min_start = max( 0, min_start - 3)
			shift = int(right_stop) - min_start
			print >> sys.stderr, "###### Shifting: -" + str(shift) + "bp"
			# Check if new breakpoints harms any gene on the left, in that case issue a warning
			new_left_stop = left_stop - shift

			for gene in left_db["annot_on_fasta"] :
				#gene = [ gene_start , gene_stop , gene_id ]
				if int(gene[0]) < int(new_left_stop) < int(gene[1]) :
					left_matched_loci_names.append(gene[2])
					#print >> sys.stderr, "Gene in the left sequence in the region overlap: " + str(gene)

			if not left_matched_loci_names == [] :
				righ_matched_loci_names = [ x[2] for x in right_matched_loci ]
				print >> sys.stdout , "[WARNING] Overlap between left region " + str(left_db["region_corrected"]) + " and right region " + str(right_db["region_corrected"]) + " leads to harm at least one annotated locus with any selection of a junction position. See standard error for more details."
				print >> sys.stderr , "[WARNING] # Overlap between left region " + str(left_db["region_corrected"]) + " and right region " + str(right_db["region_corrected"]) + " leads to harm at least one annotated locus with any selection of a junction position."
				print >> sys.stderr , "[WARNING] # This issue may arise from the original sequence fragmentation that may have lead to partial annotation"
				print >> sys.stderr , "[WARNING] # Right region loci will be preserved in the final annotation, manual curation of the missing loci is suggested as assembled sequence may still support their structure"
				print >> sys.stderr , "[WARNING] ## Left region loci involved: " + str(left_matched_loci_names)
				print >> sys.stderr , "[WARNING] ## Right region loci involved: " + str(righ_matched_loci_names)

		l_id , l_start , l_end , l_orientation = left_db["region_corrected"]
		if l_orientation == "+" :
			l_new_end = int(l_end) - shift
			new_left_db["region_trimmed"] = [ l_id , int(l_start) , l_new_end , l_orientation ]
		else :
			l_new_start = int(l_start) + shift
			new_left_db["region_trimmed"] = [ l_id , l_new_start , int(l_end) , l_orientation ]

		r_id , r_start , r_end , r_orientation = right_db["region_corrected"]
		if r_orientation == "+" :
			r_new_start = int(r_start) + int(right_stop) - shift
			new_right_db["region_trimmed"] = [ r_id , r_new_start , int(r_end) , r_orientation ]
		else :
			r_new_end = int(r_end) - int(right_stop) + shift
			new_right_db["region_trimmed"] = [ r_id , int(r_start) , r_new_end , r_orientation ]

		refining_status = "overlapping"

	return 	new_left_db , new_right_db , left_matched_loci_names , refining_status


def export_from_agp(out_prefix, no_print_fasta, agp_db, sequences, mode, sequences_len, annotation_gff3 ="") :

	out_fasta = agp2fasta( agp_db , sequences , mode )
	#### Generate FASTA
	if not no_print_fasta :
		print >> sys.stdout, "=== Writing update sequence FASTA"
		out_fasta_file = write_fasta_from_db( out_fasta , out_prefix + ".fasta")

	agp_ranges = agp2range( agp_db , "new" )
	agp_ranges.update( agp2range( agp_db , "old" ) )

	#### - GFF3 IF GIVEN - convert annotation
	if not annotation_gff3 == "" :
		print >> sys.stdout, "=== Converting coordinates"
		if mode == "old_to_new" :
			### OLD -> NEW
			translation_db = translate_from_AGP_whole_genome(agp_db)
			### Convert coordinates
			new_gff3 = translate_gff3(annotation_gff3 , translation_db , out_prefix + ".dropped_loci.txt" , out_prefix + ".multiple_copy_loci.txt" )

		elif mode == "new_to_old" :
			translation_db = translate_from_AGP_whole_genome_reverse(agp_db)
			### Convert coordinates
			new_gff3 = translate_gff3(annotation_gff3 , translation_db , out_prefix + ".dropped_loci.txt" , out_prefix + ".multiple_copy_loci.txt" )

		else :
			new_gff3 = ""

		### Write GFF3
		if not new_gff3 == "" :
			print >> sys.stdout, "=== Writing update GFF3"
			write_gff3( new_gff3 , out_prefix + ".annotation.gff3" , get_length_from_fasta_db( out_fasta ) )


def add_counts_to_dict( feat_dict ) :
	# feat_dict:
	# feat_dict[feat_name][feat_description] = []
	new_dict = {}

	description_counts = {}
	for feat_name in sorted(feat_dict.keys())	:
		for feat_description in sorted(feat_dict[feat_name].keys()) :
			if not feat_description in description_counts :
				description_counts[feat_description] = 0
			description_counts[feat_description] += 1

	for feat_name in sorted(feat_dict.keys())	:
		new_dict[feat_name] = {}
		for feat_description in sorted(feat_dict[feat_name].keys()) :
			new_dict[feat_name][feat_description] = description_counts[feat_description]

	return new_dict


def check_marker_copies( marker_hits_by_id_db , copies=2 ) :
	# marker_hits_by_id_db[marker_id]=[ hit1 , hit2, ... ]
	# with hitN == [ seq_id , int(start) , int(stop) ]
	unique_markers = {}
	multi_copy_list = []
	print >> sys.stderr, "## Checking marker copy number"
	for marker_id in sorted( marker_hits_by_id_db.keys() ) :
		marker_copies = len(marker_hits_by_id_db[marker_id])
		if marker_copies > int(copies):
			print >> sys.stderr, "### " + marker_id + " present in " + str(marker_copies) + " copies. Discarded."
			multi_copy_list.append(marker_id)
		else :
			unique_markers[marker_id] = marker_hits_by_id_db[marker_id]

	return multi_copy_list , unique_markers


def check_in_sequence_duplications(marker_hits_by_seq_db, unique_marker_hits_by_id , multi_copy_marker_list , sequence_fasta_db , annotation_db ,  sequence_agp_db , out_prefix , cores , paths , skip_qc=False , copies=2 ) :
	# sequence_agp_db can be empty
	chimeric_sequences_db = {}
	markers_itradup = []
	unique_distinct_marker_hits_by_id = {}
	#if copies > 2 :
	#	print >> sys.stderr, "[ERROR] Unsupported ploidy level"
	#	sys.exit(1)

	for marker_id in sorted(unique_marker_hits_by_id.keys()) :
		# marker_hits_by_id_db[marker_id] ==[ [ seq_id_1 , int(start_1) , int(stop_1) ] , [ seq_id_2 , int(start_2) , int(stop_2) ] ]
		if copies <= 2 :
			if len(unique_marker_hits_by_id[marker_id]) < 2 :
				unique_distinct_marker_hits_by_id[marker_id] = unique_marker_hits_by_id[marker_id]
			elif len(unique_marker_hits_by_id[marker_id]) == 2:
				match_seq_1 = unique_marker_hits_by_id[marker_id][0][0]
				match_seq_2 = unique_marker_hits_by_id[marker_id][1][0]
				if match_seq_1 == match_seq_2 :
					if match_seq_1 not in chimeric_sequences_db :
						chimeric_sequences_db[match_seq_1] = []
					chimeric_sequences_db[match_seq_1].append(marker_id)
					markers_itradup.append(marker_id)
				else :
					unique_distinct_marker_hits_by_id[marker_id] = unique_marker_hits_by_id[marker_id]
			else:
				if len(unique_marker_hits_by_id[marker_id]) > copies :
					print >> sys.stderr, "[WARNING] Marker " + marker_id + " present in >" + str(copies) + " copies (" + str(len(unique_marker_hits_by_id[marker_id])) + ") after filtering. Check code"
		else :
			match_seq = []
			for marker_hit in unique_marker_hits_by_id[marker_id] :
				match_seq_id = marker_hit[0]
				if match_seq_id in match_seq :
					if match_seq_id not in chimeric_sequences_db :
						chimeric_sequences_db[match_seq_id] = []
					chimeric_sequences_db[match_seq_id].append(marker_id)
					markers_itradup.append(marker_id)
				else :
					match_seq.append(match_seq_id)
			if not marker_id in markers_itradup :
				unique_distinct_marker_hits_by_id[marker_id] = unique_marker_hits_by_id[marker_id]


	# if there is any chimeric sequence, make reports and generate index
	if not chimeric_sequences_db == {} and not skip_qc:
		report_db = { "Input_Sequences" : { "Reports" : {} } }
		out_dir = out_prefix + ".input_sequence_QC"
		mkdir(out_dir)
		for sequence_id in chimeric_sequences_db.keys() :
			print >> sys.stderr , "### Producing intra-sequence duplication report for " + sequence_id
			report_db["Input_Sequences"]["Reports"][sequence_id] = sequence_duplication_report( sequence_id , sequence_fasta_db , annotation_db , sequence_agp_db , chimeric_sequences_db[sequence_id] , marker_hits_by_seq_db , out_dir , cores , paths )
		index_name= "index.sequence_duplication_QC.html"
		make_index_from_report_db(index_name , out_dir , out_dir , report_db)

	return chimeric_sequences_db, markers_itradup , unique_distinct_marker_hits_by_id


def make_index_from_report_db(index_file_name , workdir , out_dir , report_db) :
	# report_db for output sequences:
	# report_db[comparison]["Whole"]["html"|"pdf"|"png"] = "comparison/file.html"
	# report_db[comparison]["Reports"][queryID]["html"|"pdf"] = "comparison/queryID.on.targetID.html"
	# with comparison any of this (Hap2 may not be present):
	#	"Hap1_vs_Reference"
	#	"Reference_vs_Hap1"
	#	"Hap1_vs_Hap1" (no report)
	#	"Hap2_vs_Reference"
	#	"Reference_vs_Hap2"
	#	"Hap2_vs_Hap1"
	#	"Hap1_vs_Hap2"
	#	"Hap2_vs_Hap2" (no report)
	# All files names are reported relative to output_dir
	# report_db for input sequences:
	# report_db["Input_Sequences"]["Reports"][sequence_id]

	# For each comparison ( report_db.keys() ) make a new tab
	# For each set of plots ( report_db[comparison].keys() ) make a box
	# for each element in each set add a listing element with html and pdf

	# if plots_set == "Whole" -> add png as image to the box

	index_file_rmd = out_dir + "/" + index_file_name + ".Rmd"
	index_rmd = open( index_file_rmd , 'w')

	# Generate Rmarkdown with the file list
	print >> index_rmd , "---"
	print >> index_rmd , "title: Marker duplication QC index"
	print >> index_rmd , "output:"
	print >> index_rmd , "  flexdashboard::flex_dashboard:"
	print >> index_rmd , "    storyboard: true"
	print >> index_rmd , "    vertical_layout: scroll"
	print >> index_rmd , "    smooth_scroll: true"
	print >> index_rmd , "    runtime: shiny"
	print >> index_rmd , "--- "
	print >> index_rmd , ""
	print >> index_rmd , " <style type=\"text/css\"> "
	print >> index_rmd , ""
	print >> index_rmd , " .chart-title {"
	print >> index_rmd , "    font-size: 24px;"
	print >> index_rmd , " }"
	print >> index_rmd , " "
	print >> index_rmd , " .chart-stage {"
	print >> index_rmd , "    font-size: 16px;"
	print >> index_rmd , " }"
	print >> index_rmd , " "
	print >> index_rmd , " </style>"
	print >> index_rmd , ""
	print >> index_rmd , ""
	print >> index_rmd , '<div style="margin-top: 50px;"></div>'
	print >> index_rmd , ""
	print >> index_rmd , ""

	for comparison in sorted( report_db.keys() ) :
		comparison_title = comparison.replace("_"," ")
		print >> index_rmd , "# " + comparison_title
		print >> index_rmd , ""
		plots_set = report_db[comparison].keys()
		if "Whole" in plots_set :
			print >> index_rmd , "##"
			print >> index_rmd , ""
			print >> index_rmd , "### Whole genome dotplot {data-height=66}"
			html_file_path = report_db[comparison]["Whole"]["html"]
			png_file_path = report_db[comparison]["Whole"]["png"]
			print >> index_rmd , "[Click to investigate dotplots](" + html_file_path + ") [![Interactive html]("  + png_file_path + "){ width=96% }](" + html_file_path + ")"
			print >> index_rmd , ""
			print >> index_rmd , ""
			print >> index_rmd , "##"
			print >> index_rmd , ""
			box_height = 42 + (22*len(report_db[comparison]["Reports"].keys()))
			print >> index_rmd , "### Report by chromosome {data-height=" + str(box_height) + "}"
			for query_id in sorted(report_db[comparison]["Reports"].keys()) :
				html_file_path = report_db[comparison]["Reports"][query_id]["html"]
				pdf_file_path = report_db[comparison]["Reports"][query_id]["pdf"]
				print >> index_rmd , "* " + str(query_id) + " [html](" + html_file_path + ") | [pdf](" + pdf_file_path + ")"
			print >> index_rmd , ""

		elif comparison=="Rejected":
			print >> index_rmd , ""
			print >> index_rmd , "Comparative analysis of sequence structure between pseudomolecules and unplaced sequences to them associated but remained unplaced."
			print >> index_rmd , ""
			print >> index_rmd , "Please, select a chromosome to start from the menu."
			print >> index_rmd , ""
			#report_db["Rejected"][chr_id][seq_id]
			for plot_list_id in sorted(plots_set) :
				plot_list_title = plot_list_id.replace("_"," ")
				print >> index_rmd , ""
				print >> index_rmd , "# " + plot_list_title
				print >> index_rmd , ""
				print >> index_rmd , "Associated unplaced sequences: "
				print >> index_rmd , ""
				for seq_id in sorted(report_db[comparison][plot_list_id].keys()) :
					html_file_path = report_db[comparison][plot_list_id][seq_id]["html"]
					pdf_file_path = report_db[comparison][plot_list_id][seq_id]["pdf"]
					png_file_path = report_db[comparison][plot_list_id][seq_id]["png"]
					sequence_size = report_db[comparison][plot_list_id][seq_id]["size"]
					print >> index_rmd , "* " + str(seq_id) + ", (Length " + '{:,}'.format(int(sequence_size)) + "bp) - Plots: [html](" + html_file_path + ") | [pdf](" + pdf_file_path + ") | [png](" + png_file_path + ")"
				print >> index_rmd , ""

		else :
			for plot_list_id in sorted(plots_set) :
				plot_list_title = plot_list_id.replace("_"," ")
				print >> index_rmd , "##"
				print >> index_rmd , ""
				box_height = 42 + (22*len(report_db[comparison][plot_list_id].keys()))
				print >> index_rmd , "### " + plot_list_title + " {data-height=" + str(box_height) + "}"
				for element_id in sorted(report_db[comparison][plot_list_id].keys()) :
					html_file_path = report_db[comparison][plot_list_id][element_id]["html"]
					pdf_file_path = report_db[comparison][plot_list_id][element_id]["pdf"]
					if "png" in report_db[comparison][plot_list_id][element_id] :
						png_file_path = report_db[comparison][plot_list_id][element_id]["png"]
						print >> index_rmd , "* " + str(element_id) + " [html](" + html_file_path + ") | [pdf](" + pdf_file_path + ") | [png](" + png_file_path + ")"
					else :
						print >> index_rmd , "* " + str(element_id) + " [html](" + html_file_path + ") | [pdf](" + pdf_file_path + ")"
				print >> index_rmd , ""

	index_rmd.close()
	# Generate html index from Rmarkdown file
	log_connection = open(out_dir + "/." + index_file_name + ".Rmd" + ".conversion.log" , 'w')
	err_connection = open(out_dir + "/." + index_file_name + ".Rmd" + ".conversion.err" , 'w')
	command = "Rscript -e 'library(rmarkdown) ; rmarkdown::render( \"" + index_file_rmd + "\" , knit_root_dir = \"" + os.path.realpath(workdir) + "\" ,  output_dir = \"" + os.path.realpath(out_dir) + "\" , output_file = \"" + index_file_name + "\")'"
	print >> sys.stderr, "#### Running command: " + command
	reportProcess = subprocess.Popen( command , shell=True , stdout=log_connection , stderr=err_connection )
	output, error = reportProcess.communicate()
	log_connection.close()
	err_connection.close()
	return index_file_name


def sequence_duplication_report( seq_id , fasta_db , annotation_db , agp_db , dup_markers_list , markers_db , outdir , cores , paths) :
	# markers_db[seq_id] = [ ... , [ int(start) , int(stop) , marker_id , marker_chr , int(marker_pos) ] , ... ]

	# Prepare track files
	## Make chinks of fasta
	seq_fasta = fasta_db[seq_id]
	seq_fasta_len = {seq_id : len(fasta_db[seq_id])}
	seq_fasta_file_name = outdir + "/" + seq_id + ".fasta"
	seq_fasta_file_name = write_fasta_from_db( {seq_id : seq_fasta} , seq_fasta_file_name , False )

	seq_chunks_fasta = make_chunks_from_fasta( seq_id , seq_fasta , 1000)
	seq_chunks_fasta_file_name = outdir + "/" + seq_id + ".chunks.fasta"
	seq_chunks_fasta_file_name = write_fasta_from_db( seq_chunks_fasta , seq_chunks_fasta_file_name , False)

	## Make coords file for dotplot
	# (seq_id , seq_fasta_file_name , seq_id , seq_fasta_file_name , outdir , cores , paths , True )
	coords_file_name = seq_id + ".coords.txt"
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
		print >> sys.stderr , "[ERROR] Wrong or no path to minimap2"
		sys.exit(1)

	else :
		map_file_name = seq_id + ".selfmap_chunks.paf"
		mappingCommand = command_line + " --cs -x asm20 -r 1000 -t " + str(cores) + " "

		map_file = open( outdir + "/" + map_file_name , "w")
		mapProcess = subprocess.Popen(mappingCommand + seq_fasta_file_name + " " + seq_chunks_fasta_file_name + " ", shell=True, stdout=map_file)
		output, error = mapProcess.communicate()
		map_file.close()
		coord_table_file = make_chunk_coords_table_from_self_paf(outdir + "/" + map_file_name, outdir, outdir + "/" + seq_id + ".selfmap_chunks.txt")
		# coord_table_file => table [tID , tLen , tStart , tStop , qID , qLen , qStart , qStop , identity , match]
		filtered_coord_table = remove_self_hits(coord_table_file)
		# filtered_coord_table == [ ... , [tID , tLen , tStart , tStop , qID , qLen , qStart , qStop , identity , match] , ... ]
		write_table( filtered_coord_table , outdir + "/" + coords_file_name )

	## Write markers files
	all_markers = markers_db[seq_id]
	duplicated_markers = []
	for marker in all_markers :
		if marker[2] in dup_markers_list :
			duplicated_markers.append(marker)
	seq_all_markers_file_name = seq_id + ".all_markers.bed"
	seq_all_markers_file = open(outdir + "/" + seq_all_markers_file_name , 'w')
	for marker in sorted(all_markers) :
		print >> seq_all_markers_file, seq_id + "\t" + str(marker[0]) + "\t" + str(marker[1]) + "\t" + str(marker[2])
	seq_all_markers_file.close()
	seq_duplicated_markers_file_name = seq_id + ".duplicated_markers.bed"
	seq_duplicated_markers_file = open(outdir + "/" + seq_duplicated_markers_file_name , 'w')
	for marker in sorted(duplicated_markers) :
		print >> seq_duplicated_markers_file, seq_id + "\t" + str(marker[0]) + "\t" + str(marker[1]) + "\t" + str(marker[2])
	seq_duplicated_markers_file.close()

	## Structure from agp
	agp_table = []
	agp_table_file_name = seq_id + ".composition.txt"
	if seq_id in agp_db :
		seq_agp = agp_db[seq_id]
		for start in sorted(seq_agp.keys()) :
			Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = seq_agp[start]
			if Compnt_Type == "W" :
				agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation ])
	agp_table_file = write_table(agp_table, outdir + "/" + agp_table_file_name)

	## Run gene count analysis
	genes_hit_counts_file_name = seq_id + ".gene_count_trace.txt"
	seq_gff3 = {}
	for gene_id in annotation_db.keys() :
		gff_line , chr , start , mRNA_dict = annotation_db[gene_id]
		if chr == seq_id :
			seq_gff3[gene_id] = annotation_db[gene_id]
			#json.dump(seq_gff3[gene_id] , sys.stderr , indent = 4)
			#print >> sys.stderr , ''
	seq_gff3_file_name = outdir + "/" + seq_id + ".annotation.gff3"
	seq_gff3_file_name = write_gff3(seq_gff3 , seq_gff3_file_name , seq_fasta_len )
	## Run gmap
	CDS_file = get_sequence( seq_gff3 , {seq_id : seq_fasta} ,  outdir + "/" + seq_id , "CDS")
	### Index genome sequences
	index_dir = outdir + "/gmap_index"
	indexing_out_file = open( outdir + "/" + seq_id + ".gmap_index.log" ,"w" )
	indexing_err_file = open( outdir + "/" + seq_id + ".gmap_index.err" ,"w" )
	mkdir(index_dir)
	indexing_command = "gmap_build -D " + index_dir + " -d " + seq_id + ".fasta" + " " + seq_fasta_file_name
	indexProcess = subprocess.Popen(indexing_command, shell=True, stdout=indexing_out_file , stderr=indexing_err_file)
	output, error = indexProcess.communicate()
	indexing_out_file.close()
	indexing_err_file.close()
	### map CDS on results
	gmap_results_file_name = seq_id + ".self.gmap.gff3"
	gmap_gff3 = open( outdir + "/" + gmap_results_file_name , "w" )
	gmap_err = open( outdir + "/" + gmap_results_file_name + ".err" , "w" )
	gmapCommand = "gmap -D " + index_dir + " -d " + seq_id + ".fasta" + " -f 2 -n 500 -t " + str(cores) + " " + CDS_file
	gmapProcess = subprocess.Popen(gmapCommand, shell=True, stdout=gmap_gff3 , stderr=gmap_err)
	output, error = gmapProcess.communicate()
	gmap_err.close()
	gmap_gff3.close()
	genes_hit_counts_table = read_gmap_results(seq_gff3, outdir + "/" + gmap_results_file_name, 95, 95 , "gene" )
	genes_hit_counts = write_table(genes_hit_counts_table, outdir + "/" + genes_hit_counts_file_name )

	# Make report
	reports = {}
	report_file_html = seq_id + ".duplication_report.html"
	report_file_pdf_prefix = seq_id + ".duplication_report"
	reports["html"] = make_single_html_report( os.path.realpath(outdir) , os.path.realpath(outdir) , report_file_html , seq_id , agp_table_file_name , coords_file_name , genes_hit_counts_file_name , seq_all_markers_file_name , seq_duplicated_markers_file_name )
	reports["pdf"] , reports["png"] = make_single_pdf_report( os.path.realpath(outdir) , os.path.realpath(outdir) , report_file_pdf_prefix , seq_id , agp_table_file_name , coords_file_name , genes_hit_counts_file_name , seq_all_markers_file_name , seq_duplicated_markers_file_name )
	return reports


def remove_self_hits( hit_table_file ) :
	# hit_table_file line = [tID , tLen , tStart , tStop , qID , qLen , qStart , qStop , identity , match]
	filtered_table = []
	for line in open(hit_table_file) :
		if line.rstrip() == "" or line[0] == "#" or line[0:3] == "tID":
			continue
		else :
			tID , tLen , tStart , tStop , qID , qLen , qStart , qStop , identity , match = line.rstrip().split("\t")
			if ( int(tStart) == int(qStart) ) and ( int(tStop) == int(qStop) ) :
				continue
			else :
				filtered_table.append([tID , tLen , tStart , tStop , qID , qLen , qStart , qStop , identity , match])

	return filtered_table


def make_single_html_report( workdir, output_dir , output_file, queryID , structure , coords, gene_counts , markers_all , markers_dup ) :
	script=scriptDirectory + "/SelfReport.html.Rmd"
	log_connection = open( output_dir + "/." + queryID + ".report.html.log" , 'w')
	err_connection = open( output_dir + "/." + queryID + ".report.html.err" , 'w')
	command = "Rscript -e 'library(rmarkdown) ; rmarkdown::render(\"" + script + "\" , knit_root_dir = \"" + workdir + "\" , output_file = \"" + output_file + "\" , output_dir = \"" + output_dir + "\" , params=list(coords = \"" + coords + "\" , structure = \"" + structure + "\" , gene_counts = \"" + gene_counts + "\" , markers_all = \"" + markers_all + "\" , markers_dup = \"" + str(markers_dup) + "\" , queryID = \"" + queryID + "\"))'"
	print >> sys.stderr, "#### Running command: " + command
	reportProcess = subprocess.Popen( command , shell=True , stdout=log_connection , stderr=err_connection )
	output, error = reportProcess.communicate()
	log_connection.close()
	err_connection.close()
	return output_file


def make_single_pdf_report( workdir, output_dir , output_file , queryID , structure , coords, gene_counts , markers_all , markers_dup ) :
	# Rscript --vanilla support_scripts/SelfReport.pdf.R -d ~/Desktop -o test.pdf -c coords.txt -g gene_counts.txt -m markers_all.bed -n markers_dup.bed -q "VITMroTrayshed_v2.0.hap2.chr01" -s structure.txt
	log_connection = open( output_dir + "/." + queryID + ".report.pdf.log" , 'w')
	err_connection = open( output_dir + "/." + queryID + ".report.pdf.err" , 'w')
	script=scriptDirectory + "/SelfReport.pdf.R"
	command = "Rscript --vanilla " + script + " -d " + output_dir + " -o " + output_file + " -c " + workdir + "/" + coords + " -g " +  workdir + "/" + gene_counts + " -m " +  workdir + "/" + markers_all + " -n " +  workdir + "/" + markers_dup + " -q " + queryID + " -s " +  workdir + "/" + structure
	print >> sys.stderr, "#### Running command: " + command
	reportProcess = subprocess.Popen( command , shell=True , stdout=log_connection , stderr=err_connection )
	output, error = reportProcess.communicate()
	log_connection.close()
	err_connection.close()
	return output_file + ".pdf" , output_file + ".png"


def write_table( row_list , file_name , compress=False) :
	if compress :
		file_out = gzip.open(file_name , 'wb')
	else :
		file_out = open(file_name , 'w')
	for row in sorted(row_list) :
		print >> file_out, "\t".join([str(x) for x in row])
	file_out.close()
	return file_name


def read_table( file_name ) :
	list_of_lists = []
	try:
		for line in gzip.open( file_name , 'rb') :
			el = line.rstrip().split("\t")
			list_of_lists.append( el )
	except IOError:
		for line in open( file_name ) :
			el = line.rstrip().split("\t")
			list_of_lists.append( el )
	return list_of_lists


def clean_markers( unique_distinct_marker_hits_by_id , chimera_id_db , marker_map_order_seq , marker_map_order_id , output_prefix , fasta_db , nucmer_param , cores , remove_chimera , extend , forced_1 , forced_2 , use_forced_1 , use_forced_2) :
	# unique_distinct_marker_hits_by_id[marker_id] ==[ [ seq_id_1 , int(start_1) , int(stop_1) ] , [ seq_id_2 , int(start_2) , int(stop_2) ] ]
	filtered_hits_by_seq = {}
	temp_dir_name = output_prefix + ".marker_cleaning_tmp"
	mkdir(temp_dir_name)
	log_file_name = output_prefix + ".marker_blocks.txt"
	forced_log_file_name = output_prefix + ".marker_blocks.forced.txt"
	log_file = open(log_file_name , 'w')
	forced_log_file = open(forced_log_file_name , 'w')

	forced_chr_direction = {}
	if use_forced_1 :
		for chr_id in sorted(forced_1.keys()) :
			for seq_id_dir in forced_1[chr_id] :
				seq_id = seq_id_dir[:-2]
				direction = seq_id_dir[-1]
				forced_chr_direction[seq_id] = [ chr_id , direction]
	if use_forced_2 :
		for chr_id in sorted(forced_2.keys()) :
			for seq_id_dir in forced_2[chr_id] :
				seq_id = seq_id_dir[:-2]
				direction = seq_id_dir[-1]
				forced_chr_direction[seq_id] = [ chr_id , direction]

	# 1 - remove sequences with ambiguity, if necessary
	for marker_id in sorted(unique_distinct_marker_hits_by_id.keys()):
		for element in unique_distinct_marker_hits_by_id[marker_id] :
			seq_id = element[0]
			start = int(element[1])
			stop = int(element[2])
			pos = int(marker_map_order_id[marker_id][1])
			chr_id = marker_map_order_id[marker_id][0]
			if (not remove_chimera) or (remove_chimera and (seq_id not in chimera_id_db )) :
				if seq_id not in filtered_hits_by_seq :
					filtered_hits_by_seq[seq_id] = []
				filtered_hits_by_seq[seq_id].append( [ seq_id, start , stop , chr_id , pos , marker_id ] )
				#print >> sys.stderr , unique_distinct_marker_hits_by_id[marker_id]
			else :
				print >> sys.stderr , "[WARNING] Marker " + marker_id + " removed as located on a chimeric sequence"

	# 	2 - Assign match direction, filtering inconsistent markers
	best_marker_set = {}
	unreliable_list = {}
	for seq_id in sorted(filtered_hits_by_seq.keys()) :
		# Markers can be associated to different chromosomes and different directions -> assign to the most abundant
		sorted_hits = sorted(filtered_hits_by_seq[seq_id])
		# Split by chr, find the marker congruent by strand
		markers_by_chr = {}
		for hit in sorted_hits :
			seq_id, start , stop , chr_id , chr_pos , marker_id = hit
			if chr_id not in markers_by_chr :
				markers_by_chr[chr_id] = []
			markers_by_chr[chr_id].append([ chr_id , chr_pos , marker_id , seq_id, start , stop ])
		#json.dump(markers_by_chr , sys.stderr)

		best_paths_by_chr_by_strand = {}
		# key = tuple (chr, strand)

		for chr_id in sorted(markers_by_chr.keys()) :
			best_paths_by_chr_by_strand[ ( chr_id , "+" ) ] = find_best_marker_set(markers_by_chr[chr_id], marker_map_order_seq[chr_id], "+")
			print >> log_file , seq_id + "\t" + chr_id + "\t+\t" + ";".join([ ",".join([ str(y) for y in x]) for x in best_paths_by_chr_by_strand[ ( chr_id , "+" ) ] ])
			best_paths_by_chr_by_strand[ ( chr_id , "-" ) ] = find_best_marker_set(markers_by_chr[chr_id], marker_map_order_seq[chr_id], "-")
			print >> log_file , seq_id + "\t" + chr_id + "\t-\t" + ";".join([ ",".join([ str(y) for y in x]) for x in best_paths_by_chr_by_strand[ ( chr_id , "-" ) ] ])
		# Take best marker set
		chr_strand_markers_count = []
		for pair in sorted(best_paths_by_chr_by_strand.keys()) :
			chr_strand_markers_count.append( [ len(best_paths_by_chr_by_strand[pair]) , pair[0] , pair[1] ] )
		# Sort descending by count of markers, ascending by chr and strand
		# hits on different direction of the same chromosome are reported one after the other
		chr_strand_markers_count = sorted( chr_strand_markers_count , key=lambda k: (-k[0], k[1] , k[2]) )
		#json.dump(chr_strand_markers_count , sys.stderr, indent=4)
		#print >> sys.stderr , ""


		if seq_id in forced_chr_direction :
			chr_id , orientation = forced_chr_direction[seq_id]
			try :
				marker_list_plus = best_paths_by_chr_by_strand[ ( chr_id , "+" ) ]
				marker_pos_plus = [ int(x[1]) for x in marker_list_plus ]
			except :
				print >> sys.stdout , "[ERROR] Sequence has no markers on the requested chromosome"
				print >> sys.stderr , "[ERROR] Sequence " + seq_id + " has no markers on chromosome " + chr_id + "(+)"
				exit(1)
			try :
				marker_list_minus = best_paths_by_chr_by_strand[ ( chr_id , "-" ) ]
				marker_pos_minus = [ int(x[1]) for x in marker_list_minus ]
			except :
				print >> sys.stdout , "[ERROR] Sequence has no markers on the requested chromosome"
				print >> sys.stderr , "[ERROR] Sequence " + seq_id + " has no markers on chromosome " + chr_id + "(-)"
				exit(1)
			extended_list = marker_pos_plus + marker_pos_minus
			extended_range = [ min(extended_list) , max(extended_list) ]
			if orientation == "+" :
				range = [ min(marker_pos_plus) , max(marker_pos_plus) ]
				print >> forced_log_file , seq_id + "\t" + chr_id + "\t+\t" + ";".join([ ",".join([ str(y) for y in x]) for x in marker_list_plus ])
			else :
				range = [ min(marker_pos_minus) , max(marker_pos_minus) ]
				print >> forced_log_file , seq_id + "\t" + chr_id + "\t-\t" + ";".join([ ",".join([ str(y) for y in x]) for x in marker_list_minus ])
			best_marker_set[seq_id] = { "chr" : [ chr_id , orientation ] , "markers" : best_paths_by_chr_by_strand[ chr_id , orientation ] , "range" : range , "extended_range" : extended_range}

		else :

			first_hit = chr_strand_markers_count[0]
			second_hit = chr_strand_markers_count[1]
			# first_hit = [ highest_marker_count_block , chr_id , chr_orientation ]

			if first_hit[0] > second_hit[0] :
				# First hit counts more markers than any other -> safe assignation
				chr_id = first_hit[1]
				orientation = first_hit[2]
				marker_list_plus = best_paths_by_chr_by_strand[ chr_id , "+" ]
				marker_pos_plus = [ int(x[1]) for x in marker_list_plus ]
				marker_list_minus = best_paths_by_chr_by_strand[ chr_id , "-" ]
				marker_pos_minus = [ int(x[1]) for x in marker_list_minus ]
				extended_list = marker_pos_plus + marker_pos_minus
				extended_range = [ min(extended_list) , max(extended_list) ]

				if orientation == "+" :
					range = [ min(marker_pos_plus) , max(marker_pos_plus) ]
				else :
					range = [ min(marker_pos_minus) , max(marker_pos_minus) ]

				best_marker_set[seq_id] = { "chr" : [ chr_id , orientation ] , "markers" : best_paths_by_chr_by_strand[ chr_id , orientation ] , "range" : range , "extended_range" : extended_range}
			else :
				# 1st and 2nd hit has the same number of markers, check if unreliable
				if not first_hit[1] == second_hit[1] :
					# Different chromosomes -> unreliable positioning
					# Situation unlikely to happen given the sorting
					unreliable_list[seq_id] = []
					unreliable_list[seq_id].append( [ first_hit[1] , first_hit[2] , best_paths_by_chr_by_strand[ ( first_hit[1] , first_hit[2]) ] ] )
					unreliable_list[seq_id].append( [ second_hit[1] , second_hit[2] , best_paths_by_chr_by_strand[ ( second_hit[1] , second_hit[2]) ] ] )
				else :
					# Same chromosome -> opposite strands have the same probability, check third to see if ti has the same number of markers -> unreliable
					if len(chr_strand_markers_count) > 2 :
						third_hit = chr_strand_markers_count[2]
						if first_hit[0] > third_hit[0] :
							# 3rd hit on another another chr but with less markers -> unknown direction
							chr_id = first_hit[1]
							marker_list_plus = best_paths_by_chr_by_strand[ chr_id , "+" ]
							marker_pos_plus = [ int(x[1]) for x in marker_list_plus ]
							marker_list_minus = best_paths_by_chr_by_strand[ chr_id , "-" ]
							marker_pos_minus = [ int(x[1]) for x in marker_list_minus ]
							extended_list = marker_pos_plus + marker_pos_minus
							extended_range = [ min(extended_list) , max(extended_list) ]
							range = {}
							range["+"] = [ min(marker_pos_plus) , max(marker_pos_plus) ]
							range["-"] = [ min(marker_pos_minus) , max(marker_pos_minus) ]
							best_marker_set[seq_id] = { "chr" : [ first_hit[1] , "." ] , "markers" : {} , "range" : range , "extended_range" : extended_range}
							best_marker_set[seq_id]["markers"][first_hit[2]] = best_paths_by_chr_by_strand[ (first_hit[1] , first_hit[2]) ]
							best_marker_set[seq_id]["markers"][second_hit[2]] = best_paths_by_chr_by_strand[ (second_hit[1] , second_hit[2]) ]
						elif first_hit[0] == third_hit[0] :
							# 3rd hit must be on another chr, has the same number of markers -> unreliable
							unreliable_list[seq_id] = []
							unreliable_list[seq_id].append( [ first_hit[1] , first_hit[2] , best_paths_by_chr_by_strand[ ( first_hit[1] , first_hit[2]) ] ] )
							unreliable_list[seq_id].append( [ second_hit[1] , second_hit[2] , best_paths_by_chr_by_strand[ ( second_hit[1] , second_hit[2]) ] ] )
							unreliable_list[seq_id].append( [ third_hit[1] , third_hit[2] , best_paths_by_chr_by_strand[ ( third_hit[1] , third_hit[2]) ] ] )
					else :
						# Only 2 hits on the same chromosome
						chr_id = first_hit[1]
						marker_list_plus = best_paths_by_chr_by_strand[ chr_id , "+" ]
						marker_pos_plus = [ int(x[1]) for x in marker_list_plus ]
						marker_list_minus = best_paths_by_chr_by_strand[ chr_id , "-" ]
						marker_pos_minus = [ int(x[1]) for x in marker_list_minus ]
						extended_list = marker_pos_plus + marker_pos_minus
						extended_range = [ min(extended_list) , max(extended_list) ]
						range = {}
						range["+"] = [ min(marker_pos_plus) , max(marker_pos_plus) ]
						range["-"] = [ min(marker_pos_minus) , max(marker_pos_minus) ]
						best_marker_set[seq_id] = { "chr" : [ first_hit[1] , "." ] , "markers" : {} , "range" : range , "extended_range" : extended_range}
						best_marker_set[seq_id]["markers"][first_hit[2]] = best_paths_by_chr_by_strand[ (first_hit[1] , first_hit[2]) ]
						best_marker_set[seq_id]["markers"][second_hit[2]] = best_paths_by_chr_by_strand[ (second_hit[1] , second_hit[2]) ]
	log_file.close()
	forced_log_file.close()

	# Print all unreliable
	if not unreliable_list == {} :
		unreliable_list_file_name = output_prefix + ".unreliable_sequences.txt"
		unreliable_list_file = open(unreliable_list_file_name , 'w')
		counter = 0
		for seq_id in sorted(unreliable_list.keys()) :
			counter +=1
			#print >> unreliable_list_file , ">" + seq_id
			for element in unreliable_list[seq_id] :
				print >> unreliable_list_file, element
		print >> sys.stderr , "## Unreliable sequences because of ambiguous maker usage: " + str(counter) + ". See " + unreliable_list_file_name + " file for more details."

	# Get marker used by sequences
	marker_mapping_sequences = {}
	# Each marker should hit max 2 distinct sequences
	for seq_id in sorted(best_marker_set.keys()) :
		chr_id , direction = best_marker_set[seq_id]["chr"]
		if not direction == "." :
			# Has direction, only one list of markers is present
			marker_list = best_marker_set[seq_id]["markers"]
		else :
			# concatenate all the hits from both directions
			marker_list = []
			if "+" in  best_marker_set[seq_id]["markers"] :
				marker_list += best_marker_set[seq_id]["markers"]["+"]
			if "-" in  best_marker_set[seq_id]["markers"] :
				marker_list += best_marker_set[seq_id]["markers"]["-"]
		# marker_list = [ [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... , ]
		for hit in marker_list :
			chr_id , chr_pos , marker_id , seq_id, start , stop = hit
			if marker_id not in marker_mapping_sequences :
				marker_mapping_sequences[marker_id] = []
			if seq_id not in marker_mapping_sequences[marker_id] :
				marker_mapping_sequences[marker_id].append(seq_id)

	# 	3 - generate a range of covered markers for each sequence and recover direction where possible
	clean_hits_by_seq = {}
	print >> sys.stderr , "# Comparing unoriented sequences to alternative allele to recover orientation"
	for seq_id in sorted(best_marker_set.keys()) :
		chr_id , direction = best_marker_set[seq_id]["chr"]
		marker_list = best_marker_set[seq_id]["markers"]
		# marker_list = [ [chr_id , chr_pos , marker_id , seq_id, start , stop ] , ... , ]
		# or
		# marker_list["+"] = [ [chr_id , chr_pos , marker_id , seq_id, start , stop ] , ... , ]
		# marker_list["-"] = [ [chr_id , chr_pos , marker_id , seq_id, start , stop ] , ... , ]
		if not direction == "." :
			if extend :
				marker_range = best_marker_set[seq_id]["extended_range"]
			else :
				marker_range = best_marker_set[seq_id]["range"]
			clean_hits_by_seq[seq_id] = { "id" : seq_id , "chr" : best_marker_set[seq_id]["chr"] ,  "markers" : best_marker_set[seq_id]["markers"] , "range" : marker_range , "orientation" : direction }
		else :
			print >> sys.stderr , "## " + seq_id + " missing orientation"
			# Recover direction of unsorted sequences
			marker_ids = []
			if "+" in marker_list :
				marker_ids += [ x[2] for x in marker_list["+"] ]
			if "-" in marker_list :
				marker_ids += [ x[2] for x in marker_list["-"] ]
			# find all marker mate sequence
			sequence_marker_mates = []
			for marker_id in marker_ids :
				sequence_marker_mates += marker_mapping_sequences[marker_id]
			sequence_marker_mates = list(set(sequence_marker_mates))
			if len(sequence_marker_mates) == 1 :
				# there is no other sequence matching the same marker(s), no recover possible
				print >> sys.stderr , "### No alternative allele sequence found, orientation could not be recovered "
				marker_range = best_marker_set[seq_id]["extended_range"]
				clean_hits_by_seq[seq_id] = { "id" : seq_id , "chr" : [ best_marker_set[seq_id]["chr"][0] , direction ], "markers" : marker_list , "range" : marker_range , "orientation" : direction }
			elif len(sequence_marker_mates) > 2 :
				# markers matching multiple sequences, no recover possible
				print >> sys.stderr , "### Markers match multiple sequences on the alternative allele, orientation could not be recovered "
				marker_range = best_marker_set[seq_id]["extended_range"]
				clean_hits_by_seq[seq_id] = { "id" : seq_id , "chr" : [ best_marker_set[seq_id]["chr"][0] , direction ] , "markers" : marker_list , "range" : marker_range , "orientation" : direction }
			elif len(sequence_marker_mates) == 2 :
				# markers are shared with only one other sequence, try to map and find direction
				if sequence_marker_mates[0] == seq_id :
					mate_id = sequence_marker_mates[1]
				else :
					mate_id = sequence_marker_mates[0]
				print >> sys.stderr , "### Comparing to " + mate_id
				new_direction = update_direction( seq_id , mate_id , best_marker_set[mate_id] , fasta_db , nucmer_param , cores , temp_dir_name)
				if new_direction == "." :
					# direction cannot be updated using the mate sequence
					marker_range = best_marker_set[seq_id]["extended_range"]
					clean_hits_by_seq[seq_id] = { "id" : seq_id , "chr" : [ best_marker_set[seq_id]["chr"][0] , new_direction ] , "markers" : marker_list , "range" : marker_range , "orientation" : new_direction }
				else :
					print >> sys.stderr , "### orientation updated --> " + new_direction
					if new_direction in marker_list :
						new_marker_list = marker_list[new_direction]
					else :
						only_key = marker_list.keys()[0]
						new_marker_list = marker_list[only_key]
					if extend :
						marker_range = best_marker_set[seq_id]["extended_range"]
					else :
						marker_range = best_marker_set[seq_id]["range"][new_direction]
					clean_hits_by_seq[seq_id] = { "id" : seq_id , "chr" : [ best_marker_set[seq_id]["chr"][0] , new_direction ] , "markers" : new_marker_list , "range" : marker_range , "orientation" : new_direction }

	#os.rmdir(temp_dir_name)

	return clean_hits_by_seq


def update_direction( seq_id , mate_id , mate_marker_set_db , fasta_db , nucmer_param , cores , temp_dir) :
	mate_orientation = mate_marker_set_db["chr"][1]
	if mate_orientation == "." :
		# mate has no rientation -> cannot mutuate it
		new_orientation = "."
		print >> sys.stderr , "### Alternative allele sequence lacking orientation information, orientation could not be recovered "
	else :
		query_db = {seq_id : fasta_db[seq_id].upper()}
		query_fasta = temp_dir + "/" + str(seq_id) + ".fasta"
		write_fasta_from_db( query_db , query_fasta )
		target_db = {}
		target_db[mate_id + "|+"] = fasta_db[mate_id].upper()
		target_db[mate_id + "|-"] = str(Seq(fasta_db[mate_id]).reverse_complement()).upper()
		target_fasta = temp_dir + "/" + str(mate_id) + ".fasta"
		write_fasta_from_db( target_db , target_fasta )
		query_on_target_prefix = temp_dir + "/" + str(seq_id) + ".on." + mate_id + ".coords"
		query_on_target_coords = map_nucmer( target_fasta , query_fasta ,  cores ,  query_on_target_prefix , nucmer_param["nucmer"] , nucmer_param["show-coords"] , " --forward " , " -l -r -T -H " )
		fasta_len_dict = get_length_from_fasta_db(target_db)
		fasta_len_dict.update(get_length_from_fasta_db(query_db))
		map_results = read_nucmer_coords( query_on_target_coords )
		uniq_alignments = hits_best_tiling_path(map_results, fasta_len_dict)
		# uniq_alignments[( Tid , Qid )] = [ .. , [ Tid , int(Tstart) , int(Tstop) , Qid, int(Qstart) , int(Qstop) , int(align_length) ,  int(match_length) ]
		if ( mate_id + "|+" , seq_id ) in uniq_alignments :
			map_match_plus = sum([ int(x[7]) for x in uniq_alignments[( mate_id + "|+" ,seq_id )] ] )
		else :
			map_match_plus = 0
		if ( mate_id + "|-" , seq_id ) in uniq_alignments :
			map_match_minus = sum([ int(x[7]) for x in uniq_alignments[( mate_id + "|-" ,seq_id )] ] )
		else :
			map_match_minus = 0

		if map_match_plus > map_match_minus :
			new_orientation = mate_orientation
		elif map_match_plus < map_match_minus :
			if mate_orientation == "+" :
				new_orientation = "-"
			else :
				new_orientation = "+"
		else :
			new_orientation = "."
			print >> sys.stderr , "### Alignment on alternative allele sequence inconclusive, orientation could not be recovered"
		# Clean up
		#os.remove(query_on_target_prefix)
		#os.remove(query_on_target_prefix.split(".coords")[0] + ".delta")
		#os.remove(query_on_target_prefix.split(".coords")[0] + ".err")
		#os.remove(query_fasta)
		#os.remove(target_fasta)

	return new_orientation


def find_best_marker_set(list_matches, marker_db_by_seq, strand="+") :
	# list_matches = [ ... , [ chr_id , chr_pos , marker_id , seq_id, start , stop ] , ... ]
	# Already split by chr
	# make use of networkx
	# marker -> node
	# sorting order -> edge (weight 1)
	# 0 -> marker == weight 1
	# marker -> last == weight 1

	hit_graph = nx.DiGraph()
	hit_number = len(list_matches)
	if strand == "+" :
		# Sort by chr in chr_pos ascending order
		hits = sorted( list_matches[:] ,  key=lambda x: int(x[1]) )
	else :
		# Sort by chr in chr_pos descending order
		hits = sorted( list_matches[:] , key=lambda x: int(x[1]) , reverse=True)
	#print >> sys.stderr , "Strand: " + strand
	#print >> sys.stderr , hits

	# generate graph
	# Add start and stop nodes == first and last marker position
	start = min( [x[0] for x in marker_db_by_seq ])
	stop = max( [x[0] for x in marker_db_by_seq ])
	hit_graph.add_node("start")
	hit_graph.add_node("stop")
	# add matching markers as nodes and connect all of the to start and stop -> weight =1
	for hit in hits :
		# hit == [ chr_id , chr_pos , marker_id , seq_id, start , stop ]
		marker_id = hit[2]
		hit_graph.add_node(marker_id)
		if strand=="+" :
			hit_graph.add_edge( "start" , marker_id , join=1 )
			hit_graph.add_edge( marker_id , "stop" , join=1 )
		else:
			hit_graph.add_edge( "stop" , marker_id , join=1 )
			hit_graph.add_edge( marker_id , "start" , join=1 )

	# connect markers with edges
	# markers hits are sorted according to the position in the map and the orientation of the chromosome
	# Two consecutive markers can be placed in the same block (make an edge) if their hits are in order in the sequence
	for i in range(0,hit_number-1) :
		for j in range(i + 1,hit_number) :
			pos_left = hits[i][1]
			marker_id_left = hits[i][2]
			stop_left = hits[i][5]
			pos_right = hits[j][1]
			marker_id_right = hits[j][2]
			start_right = hits[j][4]
			if start_right >= stop_left :
				hit_graph.add_edge( marker_id_left , marker_id_right , join=1 )

	best_marker_block = nx.dag_longest_path(hit_graph, weight='join')
	best_marker_set = []
	for hit in hits :
		if hit[2] in best_marker_block :
			best_marker_set.append(hit)
			# hit = [ chr_id , chr_pos , marker_id , seq_id, start , stop ]
	return best_marker_set


def add_orientation_to_ID( markers_list ) :
	new_markers_list = []
	# markers_list = [ ... ,
	# 	{
	# 		"id" : seq_id
	# 		"chr" : chr_id
	# 		"markers" : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ]
	# 		"range" : [marker_pos_min , marker_pos_max] ,
	# 		"orientation" : ["+" or "-" or "."] }
	# 	}
	for element in markers_list :
		new_element = dict(element)
		try :
			seq_id = element["id"]
		except :
			json.dump( element , sys.stderr , indent= 4)
			sys.exit(1)
		new_element["id"] = seq_id + "|" + element["orientation"]
		new_markers_list.append(new_element)

	return new_markers_list


def validate_marker_set( markers_list , forced_list_1 , forced_list_2 , black_list_1 , black_list_2 ,  resolution_method=2 ) :
	# resolution_method
	#	1: use marker orientation, remove from forced
	#	2: update forced orientation with calculated
	#	3: move to blacklist
	#	4: reject marker information and force the requested direction
	# Update undirected hits with forced orientation

	# markers_list = [ ... ,
	# 	{
	# 		"id" : seq_id|orientation
	# 		"chr" : chr_id
	# 		"markers" : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ] / "markers" : { "+" : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ] , "-" : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ] }
	# 		"range" : [marker_pos_min , marker_pos_max] ,
	# 		"orientation" : ["+" or "-" or "."] }
	# 	}

	new_markers_list = []
	new_forced_list_1 = []
	new_forced_list_2 = []
	new_black_list_1 = black_list_1[:]
	new_black_list_2 = black_list_2[:]
	forced_list_1_undirected = [ x[:-2] for x in forced_list_1 ]
	forced_list_2_undirected = [ x[:-2] for x in forced_list_2 ]
	#print >> sys.stderr, forced_list_1_undirected
	for component in markers_list :
		seq_id = component["id"]
		seq_name = seq_id[:-2]
		seq_orientation = seq_id[-1]
		#print >> sys.stderr , [seq_id , seq_name , seq_orientation]
		if seq_name in forced_list_1_undirected :
			# sequence is in forced_1 list, check direction
			if seq_id in forced_list_1 :
				# Same direction, no issue, copy both in new lists
				new_markers_list.append(component)
				new_forced_list_1.append(seq_id)
			else:
				# Check if the sequence has no orientation in the marker list
				if seq_orientation == "." :
					# Update marker information, notify the user and copy both in new lists
					# Find the orientation in forced_list_1
					forced_seq_id = ""
					for id in forced_list_1 :
						if id[:-2] == seq_name :
							forced_seq_id=id
					new_forced_list_1.append(forced_seq_id)
					component["id"] = forced_seq_id
					component["orientation"] = forced_seq_id[-1]
					component["markers"] = component["markers"][component["orientation"]]
					new_markers_list.append(component)
					print >> sys.stderr , "#### Sequence " + seq_name + " was missing orientation. Updated using forced list orientation --> [" +  component["orientation"] + "]"
				else :
					# Conflicting directions, use resolution_method
					try :
						a = int(resolution_method)
					except :
						print >> sys.stdout , "[ERROR] Unknown conflict resolution mode " + str(resolution_method)
						sys.exit(1)
					if int(resolution_method) == 1 :
						new_markers_list.append(component)
					elif int(resolution_method) == 2 :
						new_markers_list.append(component)
						new_forced_list_1.append(seq_id)
					elif int(resolution_method) == 3 :
						new_markers_list.append(component)
						new_black_list_1 += three_orientation_list( seq_id , True )
					elif int(resolution_method) == 4 :
						forced_seq_id = ""
						for id in forced_list_1 :
							if id[:-2] == seq_name :
								forced_seq_id=id
						new_forced_list_1.append(forced_seq_id)
						new_black_list_1 += three_orientation_list(seq_id)
					else :
						print >> sys.stdout , "[ERROR] Unknown conflict resolution mode " + str(resolution_method)
						sys.exit(1)

		elif seq_name in forced_list_2_undirected :
			# sequence is in forced_1 list, check direction
			if seq_id in forced_list_2 :
				# Same direction, no issue, copy both in new lists
				new_markers_list.append(component)
				new_forced_list_2.append(seq_id)
			else:
				# Check if the sequence has no orientation in the marker list
				if seq_orientation == "." :
					# Update marker information, notify the user and copy both in new lists
					# Find the orientation in forced_list_2
					forced_seq_id = ""
					for id in forced_list_2 :
						if id[:-2] == seq_name :
							forced_seq_id=id
					new_forced_list_2.append(forced_seq_id)
					component["id"] = forced_seq_id
					component["orientation"] = forced_seq_id[-1]
					component["markers"] = component["markers"][component["orientation"]]
					new_markers_list.append(component)
					print >> sys.stderr , "#### Sequence " + seq_name + " was missing orientation. Updated using forced list orientation --> [" +  component["orientation"] + "]"
				else :
					# Conflicting directions, use resolution_method
					try :
						a = int(resolution_method)
					except :
						print >> sys.stdout , "[ERROR] Unknown conflict resolution mode " + str(resolution_method)
						sys.exit(1)
					if int(resolution_method) == 1 :
						new_markers_list.append(component)
					elif int(resolution_method) == 2 :
						new_markers_list.append(component)
						new_forced_list_2.append(seq_id)
					elif int(resolution_method) == 3 :
						new_markers_list.append(component)
						new_black_list_2 += three_orientation_list( seq_id , True )
					elif int(resolution_method) == 4 :
						forced_seq_id = ""
						for id in forced_list_2 :
							if id[:-2] == seq_name :
								forced_seq_id=id
						new_forced_list_2.append(forced_seq_id)
						new_black_list_2 += three_orientation_list(seq_id)
					else :
						print >> sys.stdout , "[ERROR] Unknown conflict resolution mode " + str(resolution_method)
						sys.exit(1)
		else :
			# Hit is not in a forced list, just add it to new_markers_list
			if seq_orientation == "." :
				print >> sys.stderr , "#### Sequence " + seq_name + " was missing orientation. Updated to --> '+' "
				component["markers"] = component["markers"]["+"]
				new_markers_list.append(component)
			else :
				new_markers_list.append(component)

		new_black_list_1 = list(set(new_black_list_1))
		new_black_list_2 = list(set(new_black_list_2))

	return new_markers_list , new_forced_list_1 , new_forced_list_2 , new_black_list_1 , new_black_list_2


def markers_to_network( markers_list , chr_id , max_size, forced_list , blacklist , outprefix ) :
	# Assumes:
	#	if a sequence is in forced_list and in markers_list, orientation should be consistent -> QC done ahead
	#	Not all forced sequences must be in makers list
	error_list = []
	for element in forced_list :
		if element in blacklist :
			error_list.append(element)
	if not error_list == [] :
		incompatibility_file_name = outprefix + ".forced_list_incompatibility.txt"
		incompatibility_file = open(incompatibility_file_name , 'w')
		for element in error_list :
			print >> incompatibility_file , element
		incompatibility_file.close()
		print >> sys.stdout , "[ERROR] Incompatibility between list of sequences requested to be used and blacklisted ones."
		print >> sys.stderr , "[ERROR] Incompatibility between list of sequences requested to be used and blacklisted ones. " + str(len(error_list)) + " sequences present in both lists. See " + outprefix + ".forced_list_incompatibility.txt file for more details."
		sys.exit(1)

	# Extract forced sequences and filter the others from from blacklist sequences hits
	# markers_list = [ ... ,
	# 	{
	# 		"id" : seq_id|+
	# 		"chr" : chr_id
	# 		"markers" : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ]
	# 		"range" : [marker_pos_min , marker_pos_max] ,
	# 		"orientation" : ["+" or "-" or "."] }
	# 	}
	markers_list_forced = []
	markers_list_clean = []
	print >> sys.stderr, "#### Elements to use in the graph"
	for element in markers_list :
		seq_id = element["id"]
		marker_num = str(len(element["markers"]))
		marker_range = "[" + str(element["range"][0]) + ":" +  str(element["range"][1]) + "]"
		if seq_id in forced_list :
			markers_list_forced.append(element)
			print >> sys.stderr, "##### " + seq_id + ": forced - " + marker_num + " marker(s) - range " + marker_range + " - makers " + str(element["markers"])
		else :
			if not seq_id in blacklist :
				markers_list_clean.append(element)
				print >> sys.stderr, "##### " + seq_id + ": usable - " + marker_num + " marker(s) - range " + marker_range + " - makers " + str(element["markers"])

	forced_length = len(markers_list_forced)
	clean_length = len(markers_list_clean)
	# Generate nodes
	marker_grph = nx.DiGraph()
	marker_grph.add_node("ChrStart")
	marker_grph.add_node("ChrStop")

	# Generate forced block edges and connect all forced blocks to "ChrStart" == [0 , 0]  and  "ChrStop" == [stop_pos , stop_pos ]
	if forced_length > 0 :
		covered_ranges = []

		# Add forced markers first and connect them to start and stop
		for i in range(0, forced_length) :
			marker = markers_list_forced[i]
			marker_num = len(marker["markers"])
			seq_id = marker["id"]
			marker_pos_min , marker_pos_max = marker["range"]
			if marker_pos_min == marker_pos_max :
				marker_pos_min = float(marker_pos_max)
				marker_pos_max = marker_pos_min + 0.1
			else :
				marker_pos_min = float(marker_pos_min)
				marker_pos_max = float(marker_pos_max) + 0.1
			marker_grph.add_node(marker_pos_min)
			marker_grph.add_node(marker_pos_max)
			marker_grph.add_edge( "ChrStart" , marker_pos_min, name="" , start="ChrStart" , stop=seq_id , marker_num=0)
			marker_grph.add_edge( marker_pos_max , "ChrStop" , name="" , start=seq_id , stop="ChrStop" , marker_num=0)
			marker_grph.add_edge( marker_pos_min , marker_pos_max , name=seq_id , start=seq_id , stop=seq_id , marker_num=marker_num )
			covered_ranges.append([marker_pos_min , marker_pos_max])
		# Connect forced markers
		for i in range(0, forced_length) :
			for j in range(0, forced_length) :
				if i == j :
					continue
				else :
					# make link if i(marker_pos_min , marker_pos_max) comes strictly before j(marker_pos_min , marker_pos_max) are compatible
					marker_1 = markers_list_forced[i]
					marker_1_seq_id = marker_1["id"]
					marker_1_pos_min , marker_1_pos_max = marker_1["range"]
					if marker_1_pos_min == marker_1_pos_max :
						marker_1_pos_min = float(marker_1_pos_max)
						marker_1_pos_max = marker_1_pos_min + 0.1
					else :
						marker_1_pos_min = float(marker_1_pos_min)
						marker_1_pos_max = float(marker_1_pos_max) + 0.1

					marker_2 = markers_list_forced[j]
					marker_2_seq_id = marker_2["id"]
					marker_2_pos_min , marker_2_pos_max = marker_2["range"]
					if marker_2_pos_min == marker_2_pos_max :
						marker_2_pos_min = float(marker_2_pos_max)
						marker_2_pos_max = marker_2_pos_min + 0.1
					else :
						marker_2_pos_min = float(marker_2_pos_min)
						marker_2_pos_max = float(marker_2_pos_max) + 0.1

					if marker_1_pos_max < marker_2_pos_min :
						# Check if the link is already present
						if marker_grph.has_edge(marker_1_pos_max , marker_2_pos_min) :
							# if the link is not another update removing info on start and stop
							if not marker_grph[marker_1_pos_max][marker_2_pos_min]["name"] == "" :
								marker_grph[marker_1_pos_max][marker_2_pos_min]["start"] = ""
								marker_grph[marker_1_pos_max][marker_2_pos_min]["stop"] = ""
						else :
							marker_grph.add_edge( marker_1_pos_max , marker_2_pos_min , name= "" , start=marker_1_seq_id , stop=marker_2_seq_id , marker_num=0)

		# Update with markers_list_clean that do not overlap covered ranges
		for k in range(0, clean_length) :
			marker_c = markers_list_clean[k]
			marker_c_seq_id = marker_c["id"]
			marker_c_num = len(marker_c["markers"])
			marker_c_pos_min , marker_c_pos_max = marker_c["range"]
			if marker_c_pos_min == marker_c_pos_max :
				marker_c_pos_min = float(marker_c_pos_max)
				marker_c_pos_max = marker_c_pos_min + 0.1
			else :
				marker_c_pos_min = float(marker_c_pos_min)
				marker_c_pos_max = float(marker_c_pos_max) + 0.1
			#print >> sys.stderr , "- Checking: " +  marker_c_seq_id + " (" + str(marker_c_num) + " markers, corrected range: [" + str(marker_c_pos_min) + "," + str(marker_c_pos_max) + "] )"

			if range_overlap_any( [marker_c_pos_min , marker_c_pos_max] , covered_ranges ) :
				#print >> sys.stderr , "-- Range overlaps forced sequences ranges"
				continue
			else:
				# Add marker_c sequence to the network
				#print >> sys.stderr , "-- Adding sequence to the tiling network"
				marker_grph.add_node(marker_c_pos_min)
				marker_grph.add_node(marker_c_pos_max)
				marker_grph.add_edge( "ChrStart" , marker_c_pos_min, name="" , start="ChrStart" , stop=marker_c_seq_id , marker_num=0)
				marker_grph.add_edge( marker_c_pos_max , "ChrStop" , name="" , start=marker_c_seq_id , stop="ChrStop" , marker_num=0)
				marker_grph.add_edge( marker_c_pos_min , marker_c_pos_max , name=marker_c_seq_id , start=marker_c_seq_id , stop=marker_c_seq_id , marker_num=marker_c_num )

				# link to all compatible forced sequences
				for i in range(0, forced_length) :
					marker_f = markers_list_forced[i]
					marker_f_seq_id = marker_f["id"]
					marker_f_pos_min , marker_f_pos_max = marker_f["range"]
					if marker_f_pos_min == marker_f_pos_max :
						marker_f_pos_min = float(marker_f_pos_max)
						marker_f_pos_max = marker_f_pos_min + 0.1
					else :
						marker_f_pos_min = float(marker_f_pos_min)
						marker_f_pos_max = float(marker_f_pos_max) + 0.1
					#print >> sys.stderr , "--- Checking connectivity with forced sequence: " + marker_f_seq_id + "(range: " + str(marker_f_pos_min) + "-" +  str(marker_f_pos_max) + ")"
					if marker_f_pos_max < marker_c_pos_min :
						# marker_c follows marker_f
						# Check if the link marker_f_pos_max -> marker_c_pos_min is already present
						#print >> sys.stderr , "---- Linking downstream of forced sequence"
						if marker_grph.has_edge(marker_f_pos_max , marker_c_pos_min) :
							##print >> sys.stderr , "----- Link was already present"
							# if the link is not another update removing info on start and stop
							if not marker_grph[marker_f_pos_max][marker_c_pos_min]["name"] == "" :
								marker_grph[marker_f_pos_max][marker_c_pos_min]["start"] = ""
								marker_grph[marker_f_pos_max][marker_c_pos_min]["stop"] = ""
						else :
							#print >> sys.stderr , "----- Link created"
							#print >> sys.stderr , "----- Link: " + marker_f_seq_id + " --> " + marker_c_seq_id
							marker_grph.add_edge( marker_f_pos_max , marker_c_pos_min , name= "" , start=marker_f_seq_id , stop=marker_c_seq_id , marker_num=0)
					elif marker_c_pos_max < marker_f_pos_min :
						#print >> sys.stderr , "---- Linking upstream of forced sequence"
						# marker_f follows marker_c
						# Check if the link marker_c_pos_max -> marker_f_pos_min is already present
						if marker_grph.has_edge(marker_c_pos_max , marker_f_pos_min) :
							#print >> sys.stderr , "----- Link updated"
							# if the link is not another update removing info on start and stop
							if not marker_grph[marker_c_pos_max][marker_f_pos_min]["name"] == "" :
								marker_grph[marker_c_pos_max][marker_f_pos_min]["start"] = ""
								marker_grph[marker_c_pos_max][marker_f_pos_min]["stop"] = ""
						else :
							#print >> sys.stderr , "----- Link created"
							#print >> sys.stderr , "----- Link: " + marker_c_seq_id + " --> " + marker_f_seq_id
							marker_grph.add_edge( marker_c_pos_max , marker_f_pos_min , name= "" , start=marker_c_seq_id , stop=marker_f_seq_id , marker_num=0)

		# Update with links between compatible clean markers
		for i in range(0, clean_length) :
			for j in range(0, clean_length) :
				if i == j :
					continue
				else :
					# make link if i(marker_pos_min , marker_pos_max) comes strictly before j(marker_pos_min , marker_pos_max)
					marker_1 = markers_list_clean[i]
					marker_1_seq_id = marker_1["id"]
					marker_1_pos_min , marker_1_pos_max = marker_1["range"]
					if marker_1_pos_min == marker_1_pos_max :
						marker_1_pos_min = float(marker_1_pos_max)
						marker_1_pos_max = marker_1_pos_min + 0.1
					else :
						marker_1_pos_min = float(marker_1_pos_min)
						marker_1_pos_max = float(marker_1_pos_max) + 0.1

					marker_2 = markers_list_clean[j]
					marker_2_seq_id = marker_2["id"]
					marker_2_pos_min , marker_2_pos_max = marker_2["range"]
					if marker_2_pos_min == marker_2_pos_max :
						marker_2_pos_min = float(marker_2_pos_max)
						marker_2_pos_max = marker_2_pos_min + 0.1
					else :
						marker_2_pos_min = float(marker_2_pos_min)
						marker_2_pos_max = float(marker_2_pos_max) + 0.1

					#print >> sys.stderr , "---- Analysing : " + marker_1_seq_id + "(" + str(marker_1_pos_min) + "-" + str( marker_1_pos_max) + ") --> " + marker_2_seq_id + "(" + str(marker_2_pos_min) + "-" + str( marker_2_pos_max) + ")"
					if marker_1_pos_max < marker_2_pos_min :
						# Check if the link is already present
						if marker_grph.has_edge(marker_1_pos_max , marker_2_pos_min) :
							#print >> sys.stderr , "----- Link updated"
							# if the link is not another update removing info on start and stop
							if not marker_grph[marker_1_pos_max][marker_2_pos_min]["name"] == "" :
								marker_grph[marker_1_pos_max][marker_2_pos_min]["start"] = ""
								marker_grph[marker_1_pos_max][marker_2_pos_min]["stop"] = ""
						else :
							#print >> sys.stderr , "----- Link created"
							#print >> sys.stderr , "----- Link: " + marker_1_seq_id + " --> " + marker_2_seq_id
							marker_grph.add_edge( marker_1_pos_max , marker_2_pos_min , name= "" , start=marker_1_seq_id , stop=marker_2_seq_id , marker_num=0)

	else :
		# no forced, build over all clean sequences
		# Make nodes and connect all sequence to chr start and chr stop
		for k in range(0, clean_length) :
			marker_c = markers_list_clean[k]
			marker_c_seq_id = marker_c["id"]
			marker_c_num = len(marker_c["markers"])
			marker_c_pos_min , marker_c_pos_max = marker_c["range"]

			if marker_c_pos_min == marker_c_pos_max :
				marker_c_pos_min = float(marker_c_pos_max)
				marker_c_pos_max = marker_c_pos_min + 0.1
			else :
				marker_c_pos_min = float(marker_c_pos_min)
				marker_c_pos_max = float(marker_c_pos_max) + 0.1

			# Add marker_c sequence to the network
			marker_grph.add_node(marker_c_pos_min)
			marker_grph.add_node(marker_c_pos_max)
			marker_grph.add_edge( "ChrStart" , marker_c_pos_min, name="" , start="ChrStart" , stop=marker_c_seq_id , marker_num=0)
			marker_grph.add_edge( marker_c_pos_max , "ChrStop" , name="" , start=marker_c_seq_id , stop="ChrStop" , marker_num=0)
			marker_grph.add_edge( marker_c_pos_min , marker_c_pos_max , name=marker_c_seq_id , start=marker_c_seq_id , stop=marker_c_seq_id , marker_num=marker_c_num )

		# link to all compatible sequences
		for i in range(0, clean_length) :
			for j in range(0, clean_length) :
				if i == j :
					continue
				else :
					# make link if i(marker_pos_min , marker_pos_max) comes strictly before j(marker_pos_min , marker_pos_max) are compatible
					marker_1 = markers_list_clean[i]
					marker_1_seq_id = marker_1["id"]
					marker_1_pos_min , marker_1_pos_max = marker_1["range"]
					if marker_1_pos_min == marker_1_pos_max :
						marker_1_pos_min = float(marker_1_pos_max)
						marker_1_pos_max = marker_1_pos_min + 0.1
					else :
						marker_1_pos_min = float(marker_1_pos_min)
						marker_1_pos_max = float(marker_1_pos_max) + 0.1

					marker_2 = markers_list_clean[j]
					marker_2_seq_id = marker_2["id"]
					marker_2_pos_min , marker_2_pos_max = marker_2["range"]
					if marker_2_pos_min == marker_2_pos_max :
						marker_2_pos_min = float(marker_2_pos_max)
						marker_2_pos_max = marker_2_pos_min + 0.1
					else :
						marker_2_pos_min = float(marker_2_pos_min)
						marker_2_pos_max = float(marker_2_pos_max)

					if marker_1_pos_max < marker_2_pos_min :
						# Check if the link is already present
						if marker_grph.has_edge(marker_1_pos_max , marker_2_pos_min) :
							# if the link is not another update removing info on start and stop
							if not marker_grph[marker_1_pos_max][marker_2_pos_min]["name"] == "" :
								marker_grph[marker_1_pos_max][marker_2_pos_min]["start"] = ""
								marker_grph[marker_1_pos_max][marker_2_pos_min]["stop"] = ""
						else :
							#print >> sys.stderr , "----- Link: " + marker_1_seq_id + " --> " + marker_2_seq_id
							marker_grph.add_edge( marker_1_pos_max , marker_2_pos_min , name= "" , start=marker_1_seq_id , stop=marker_2_seq_id , marker_num=0)

	return marker_grph


def range_overlap_any( range , list_of_ranges) :
	# range = [start , stop]
	# list_of_ranges = [ ... , [start , stop] , ... ]
	# Test if the given range overlaps any of the ranges in list_of_ranges
	match=False
	for test_range in list_of_ranges :
		if int(test_range[0]) <= int(range[0]) <= int(test_range[1]) :
			match = True
		if int(test_range[0]) <= int(range[1]) <= int(test_range[1]) :
			match = True
		if int(range[0]) <= int(test_range[0]) <= int(range[1]) :
			match = True
		if int(range[0]) <= int(test_range[1]) <= int(range[1]) :
			match = True

	return match


def three_orientation_list( id , has_orientation=True ) :
	if has_orientation :
		clean_name = id[:-1]
	else :
		clean_name = id + "|"
	return [ clean_name + "+" , clean_name + "-" , clean_name + "." ]


def search_unwanted( edges_names , unwanted_pairs , query_len , forced_list ) :
	# Run search of the first unwanted pair to fix,
	# If a pair is found in the path:
	# 		If one is forced in -> return the other
	# 		else -> return the shortest
	output_list = ""
	for unwanted_id in sorted(unwanted_pairs.keys()) :
		unwanted_mates = unwanted_pairs[unwanted_id]
		for unwanted_mate in unwanted_mates :
			if (unwanted_id in edges_names) and (unwanted_mate in edges_names) :
				# both sequence in pair present in the path
				output_list = {}
				if unwanted_id in forced_list :
					output_list["keep"] = unwanted_id
					output_list["remove"] = three_orientation_list(unwanted_mate)
				elif unwanted_mate in forced_list :
					output_list["keep"] = unwanted_mate
					output_list["remove"] = three_orientation_list(unwanted_id)
				else :
					unwanted_id_len = query_len[unwanted_id]
					unwanted_mate_len = query_len[unwanted_mate]
					if unwanted_id_len < unwanted_mate_len :
						output_list["keep"] = unwanted_mate
						output_list["remove"] = three_orientation_list(unwanted_id)
					else :
						output_list["keep"] = unwanted_id
						output_list["remove"] = three_orientation_list(unwanted_mate)
			if not output_list == "" :
				break
		if not output_list == "" :
			# Found unwanted pair in the tapth, return it
			break

	# returns "" if no issue is found, a db with the seq_id to keep and the seq_id to remove otherwise
	# output_list["keep"] = keep_seq_id|+
	# output_list["remove"] = [ reject_seq_id|+ , reject_seq_id|- , reject_seq_id|. ]
	return output_list


def get_marker_subgraph_from_path( original_graph , selected_path ) :
	# selected_path = sorted list of nodes in the path
	edges = []
	edges_names=[]

	for i in range( 0 , len(selected_path) -1 ):
		from_node = selected_path[i]
		to_node = selected_path[i+1]
		element = original_graph[from_node][to_node]
		if not element["name"] == "" :
			edges_names.append(element["name"])
			edges.append([ element["name"] , "(" + str(from_node) + ":" + str(to_node) + ")" ,  element["marker_num"] ])

	return edges , edges_names


def make_list_from_marker_path( querylist ) :
	### Query list element format
	# Object:	[b40-14.HS_iter4_seq3733|+	,	(61252:6463804)	,	5 				]
	#			[ID|strand					,	(T_start:Tstop)	,	# of markers	]
	id_list = []
	seq_list = []
	for CompntId in querylist :
		Id, T_range , markers = CompntId
		Id_clean = Id[:-1]
		seq_list.append(Id)
		id_list.append(Id_clean+"+")
		id_list.append(Id_clean+"-")
		id_list.append(Id_clean+".")
	return id_list , seq_list


def fill_orientation( querylist , id_list , report_file_name , new_orientation = "+" ) :
	# querylist[chr_id]	= [ ... ,	[b40-14.HS_iter4_seq3733|+	,	(61252:6463804)	,	5 			]	, ... ]
	#								[ID|strand					,	(T_start:Tstop)	,	# of markers]
	# id_list[chr_id] = [ ... , "ID|strand" , ...  ]

	# Copy all info on the query sequences
	new_querylist = dict(querylist)
	new_id_list = dict(id_list)
	report_file = open(report_file_name ,"w")
	# Update query orientation information if it was missing
	for chr_id in sorted(new_querylist.keys()) :
		if (new_orientation == "+") or (new_orientation == "-") :
			# Fixed update
			item_number = len(new_querylist[chr_id])
			for i in range(0,item_number) :
				seq_id = new_querylist[chr_id][i][0]
				seq_id_clean = seq_id[:-2]
				strand = seq_id[-1]
				if strand == "." :
					new_querylist[chr_id][i][0] = seq_id_clean + "|" + new_orientation
					new_id_list[chr_id][i] = seq_id_clean + "|" + new_orientation
					#print >> report_file , seq_id_clean + " strand -> " + new_orientation
					#print >> sys.stderr , "### " + seq_id_clean + " strand -> " + new_orientation
		else :
			# Update from db
			# new_orientation[seq_id] = orientation
			item_number = len(new_querylist[chr_id])
			for i in range(0,item_number) :
				seq_id = new_querylist[chr_id][i][0]
				seq_id_clean = seq_id[:-2]
				strand = seq_id[-1]
				if strand == "." :
					if seq_id_clean in new_orientation :
						orientation = new_orientation[seq_id_clean]
					else :
						print >> sys.stderr, "[WARNING] Sequence " + seq_id_clean + " could not a proper orientation. Orientation will be forced to '+' "
						orientation = "+"
					new_querylist[chr_id][i][0] = seq_id_clean + "|" + orientation
					new_id_list[chr_id][i] = seq_id_clean + "|" + orientation
					#print >> report_file , seq_id_clean + " strand -> " + orientation
					#print >> sys.stderr , "### " + seq_id_clean + " strand -> " + orientation

	report_file.close()
	return new_querylist , new_id_list


def remove_sequence_from_graph( unwanted_ids_list , chr_id , max_size, markers_list , forced_list , blacklist ) :
	new_blacklist = blacklist[:]
	#print >> sys.stderr, "blacklist:"
	#print >> sys.stderr, blacklist
	#print >> sys.stderr, new_blacklist
	for element in unwanted_ids_list :
		new_blacklist.append(element)
	print >> sys.stderr, "##### used blacklist: " + str(new_blacklist)
	print >> sys.stderr, "##### used forced list: " + str(forced_list)
	new_graph = markers_to_network( markers_list , chr_id , max_size, forced_list , new_blacklist , "incompatibility" )
	return new_graph


def get_no_orientation( list_1 ,  list_2 ) :
	seq_id_db = {}
	for chr_id in list_1 :
		for seq_id in list_1[chr_id] :
			strand = seq_id[-1]
			if strand == "." :
				seq_id_clean = seq_id[:-2]
				seq_id_db[seq_id_clean] = chr_id
	for chr_id in list_2 :
		for seq_id in list_2[chr_id] :
			strand = seq_id[-1]
			if strand == "." :
				seq_id_clean = seq_id[:-2]
				seq_id_db[seq_id_clean] = chr_id
	# Return db seq_id -> chr
	return seq_id_db


def best_orientation( alignments ) :
	# alignments[Tid , Qid] =[ ... , [ Tid , int(Tstart) , int(Tstop) , Qid, int(Qstart) , int(Qstop) , int(align_length) ,  int(match_length) ] , ... ]

	orientation = "."
	if len(alignments.keys()) > 2 :
		print >> sys.stdout , "[ERROR] Expected one alignment for each orientation, " + str(len(alignments.keys())) + " found"
		print >> sys.stderr , "[ERROR] Expected one alignment for each orientation, " + str(len(alignments.keys())) + " found"
		print >> sys.stderr , "[ERROR] Reported matches:"
		for pair in alignments.keys() :
			print >> sys.stderr , pair
		sys.exit(1)
	elif len(alignments.keys()) == 0 :
		print >> sys.stdout , "[WARNING] No mapping found for sequence on chromosome, orientation couldn't be fixed"
	elif len(alignments.keys()) == 1 :
		match = alignments.keys()[0]
		best_query = match[1]
		orientation = best_query[-1]
	else :
		# 2 hits, find the one with highest identity
		match_1_id , match_2_id = alignments.keys()
		match_1_matches = sum([int(x[7]) for x in alignments[match_1_id]] )
		match_2_matches = sum([int(x[7]) for x in alignments[match_2_id]] )
		if match_1_matches > match_2_matches :
			orientation = match_1_id[1][-1]
		elif match_1_matches < match_2_matches :
			orientation = match_2_id[1][-1]
		else :
			# match_1_matches == match_2_matches
			print >> sys.stdout , "[WARNING] Mapping with the same identity on both strands, orientation couldn't be fixed"
	return orientation


def generate_fasta_from_path( paths_edges_db , tmp_dir , prefix , fasta_db , gap=1000) :
	gap_seq = "N"*int(gap)
	info_db ={}
	for chr_id in sorted(paths_edges_db.keys()) :
		info_db[chr_id] = {}
		new_id = prefix + "_" + chr_id
		print >> sys.stderr , "## " + chr_id + " | " +  new_id
		out_fasta_db = {new_id : ""}
		info_db[chr_id]["id"] = new_id
		info_db[chr_id]["structure"] = []
		for seq_oriented_id in paths_edges_db[chr_id] :
			#print >> sys.stderr, seq_oriented_id
			seq_id = seq_oriented_id[:-2]
			#print >> sys.stderr, seq_id
			orientation = seq_oriented_id[-1]
			#print >> sys.stderr, orientation
			if seq_id in fasta_db :
				seq_fasta = fasta_db[seq_id]
				if orientation == "+" :
					seq_fasta = fasta_db[seq_id]
				elif orientation == "-" :
					seq_fasta = str(Seq(fasta_db[seq_id]).reverse_complement())
				else :
					print >> sys.stderr , "[WARNING] Sequence requested without orientation (" + seq_id + "), using forward (+) orientation for it"
					seq_fasta = fasta_db[seq_id]
			else :
				print >> sys.stdout , "[ERROR] Unknown sequence requested"
				print >> sys.stderr , "[ERROR] Unknown sequence requested (" + seq_id + ")"
				print >> sys.stderr , fasta_db.keys()
				sys.exit(1)

			if not len(out_fasta_db[new_id]) == 0 :
				#I have previous sequences in the output fasta, add a gap
				out_fasta_db[new_id] += gap_seq

			out_start = len(out_fasta_db[new_id])
			out_fasta_db[new_id] += seq_fasta
			out_stop = len(out_fasta_db[new_id])
			info_db[chr_id]["structure"].append([ out_start , out_stop , seq_oriented_id ])
			print >> sys.stderr , "### " + seq_oriented_id + " (len: " + str(len(seq_fasta)) +  ") -> " + new_id + ":" + str(out_start) + "-" + str(out_stop)
		info_db[chr_id]["fasta_len"] = len(out_fasta_db[new_id])
		fasta_file = tmp_dir + "/" + new_id
		info_db[chr_id]["fasta_file"] = write_fasta_from_db( out_fasta_db , fasta_file , False)
	# info_db[chr_id]
	# info_db[chr_id]["id"]
	# info_db[chr_id]["fasta_file"]
	# info_db[chr_id]["fasta_len"]
	# info_db[chr_id]["structure"] = [ ... , [concat_start , concat_stop , component_seq_id ] , ...  ]
	return info_db


def get_component_alignment( map_db ) :
	# map_db[chr_id]
	# map_db[chr_id]["id"] = intermediate_Qid
	# map_db[chr_id]["fasta_file"]
	# map_db[chr_id]["fasta_len"]
	# map_db[chr_id]["structure"] = [ ... , [concat_start , concat_stop , component_seq_id ] , ...  ]
	# map_db[chr_id]["best_alignment"] = file >> [ ... , [ chr_id , int(Tstart) , int(Tstop) , intermediate_Qid , int(Qstart) , int(Qstop) , int(align_length) ,  int(match_length) ] , ... ]
	mapped_sequences = {}
	unmapped_sequences = {}
	splitted_regions = {}
	checkpoint_file = open(".get_component_alignment.checkpoint.txt" , 'w')
	for chr_id in sorted(map_db.keys()):
		splitted_regions[chr_id] = []
		mapped_sequences[chr_id] = []
		print >> sys.stderr, "### Splitting " + chr_id + " intermediate alignments across query sequences"
		region_list = sorted(map_db[chr_id]["structure"])
		# print >> sys.stderr, "### Regions: " + str(len(region_list))
		# region_list is sorted by intermediate_Qid start and stop
		hits = read_table(map_db[chr_id]["best_alignment"])
		print >> sys.stderr, "#### Hits: " + str(len(hits))
		# sort hits by intermediate_Qid start and stop
		#hits.sort(key = operator.itemgetter( 4, 5) )
		# Find mapping position for each region
		for region in region_list:
			#print >> sys.stderr, "## Region analysed : " + str(region)
			region_start , region_stop,  region_id = region
			# Find all hits relative to the ROI
			hits_on_ROI = []
			for hit in hits :
				chr_id , Tstart , Tstop , intermediate_Qid , Qstart , Qstop , align_length ,  match_length = hit
				if int(Qstart) > int(region_stop) or int(Qstop) < int(region_start) :
					# hit outside fo ROI
					#print >> sys.stderr, "#### Outside ROI"
					continue
				else :
					# TODO: SOMETHING MAY BE WRONG WITH COORDINATES TRANSLATION
					# Cap hit within the sequence region
					#print >> sys.stderr, "#### Hit within ROI range: " + str(hit)
					hit_portion_start = max( int(Qstart) , int(region_start) )
					hit_portion_stop = min( int(Qstop) , int(region_stop) )
					on_seq_start = hit_portion_start - region_start
					on_seq_stop = hit_portion_stop - region_start

					delta_start = max( 0 , int(region_start) - int(Qstart))
					delta_stop = max( 0 , int(Qstop) - int(region_stop))
					T_projection_start = int(Tstart) + delta_start
					T_projection_stop = int(Tstop) - delta_stop

					hit_len = int(Qstop) - int(Qstart)
					hit_portion_len = hit_portion_stop - hit_portion_start
					try :
						fraction = float(hit_portion_len) / float(hit_len)
					except :
						print >> sys.stderr , "[ERROR] in translation of hit coordinates"
						print >> sys.stderr , "[ERROR] Region: region_start , region_stop,  region_id "
						print >> sys.stderr , region
						print >> sys.stderr , "[ERROR] Hit: chr_id , Tstart , Tstop , intermediate_Qid , Qstart , Qstop , align_length ,  match_length "
						print >> sys.stderr , hit
						exit(4)

					hit_portion_aligned = int(float(align_length)*fraction)
					hit_portion_matched = int(float(match_length)*fraction)

					hits_on_ROI.append( [ chr_id , int(T_projection_start) , int(T_projection_stop) , region_id , int(on_seq_start) , int(on_seq_stop) ,  hit_portion_aligned , hit_portion_matched ] )
					#print >> sys.stderr, "#### Translated hit: " + str( [ chr_id , int(T_projection_start) , int(T_projection_stop) , region_id , int(on_seq_start) , int(on_seq_stop) ,  hit_portion_aligned , hit_portion_matched ] )
			if hits_on_ROI == [] :
				# No hit were reported for the given sequence -> issue a warning
				print >> sys.stderr, "#### [WARNING] Sequence " + region_id + " was placed with markers but couldn't find a proper alignment on the guide genome. The sequence will not appear in the results"
				if chr_id not in unmapped_sequences :
					unmapped_sequences[chr_id] = []
				unmapped_sequences[chr_id].append(region_id)
			else :
				#print >> sys.stderr, "#### Hits in the ROI to translate in mapping range: " + str(len(hits_on_ROI))
				mapped_sequences[chr_id].append(region_id)
				# Uniquify
				Tstart = min( [ int(x[1]) for x in hits_on_ROI ])
				Tstop = max( [ int(x[2]) for x in hits_on_ROI ])
				Qstart = min( [ int(x[4]) for x in hits_on_ROI ] )
				Qstop = max( [ int(x[5]) for x in hits_on_ROI ])
				matches = sum( [ int(x[6]) for x in hits_on_ROI ] )
				hitLen = sum( [ int(x[7]) for x in hits_on_ROI ] )
				splitted_regions[chr_id].append( [region_id , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ])
				print >> checkpoint_file , "\t".join([region_id , str(Tstart) , str(Tstop) , str(Qstart) , str(Qstop) , str(matches) , str(hitLen) ])
	# splitted_regions[Chr] =
	# 		# 	[ 	... ,
	# 		# 		[ sed_id|orintation , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ] ,
	# 	 	#		... ]
	# mapped_sequences[chr] = [ , sed_id|orintation , ]
	checkpoint_file.close()
	return splitted_regions , mapped_sequences , unmapped_sequences


def merge_sorted_lists( paths_edges , forced_list ) :
	new_list = {}
	for chr_id in sorted(set(paths_edges.keys() + forced_list.keys())) :
		new_list[chr_id] = []

		forced_id_graph = nx.DiGraph()
		if chr_id in forced_list :
			# forced_id_graph has paths between ids only if there is any ids to connect
			forced_id_graph.add_node("BEGIN")
			forced_id_graph.add_node("END")
			prev_id = "BEGIN"
			# Add all forced nodes and links
			for forced_id in forced_list[chr_id] :
				forced_id_graph.add_node(forced_id)
				forced_id_graph.add_edge(prev_id , forced_id , weight=1)
				prev_id = forced_id
			forced_id_graph.add_edge(prev_id , "END" , weight=1)

		sequence_id_graph = nx.DiGraph()
		sequence_id_graph.add_node("BEGIN")
		# Superimpose paths_edges links
		if chr_id in paths_edges :
			paths_edges_list = paths_edges[chr_id]
			paths_edges_list.append("END")
		else :
			paths_edges_list = ["END"]
			# If not path was calculated for the given chromosome, all forced ids are used
		prev_id = "BEGIN"
		for seq_id in paths_edges[chr_id] :
			if not forced_id_graph.has_node(seq_id) :
				# The sequence was not present in the forced list, add edge
				sequence_id_graph.add_node(seq_id)
				sequence_id_graph.add_edge(prev_id , seq_id , weight=1)
				prev_id = seq_id
			else :
				if not nx.has_path(forced_id_graph , prev_id , seq_id) :
					sequence_id_graph.add_node(seq_id)
					sequence_id_graph.add_edge(prev_id , seq_id , weight=1)
					prev_id = seq_id
				else :
					forced_path = nx.shortest_path(forced_id_graph , prev_id , seq_id)
					# is the list of nodes [prev_id , ... , seq_id]
					prev_id= forced_path[0]
					for node in range(1, len(forced_path)):
						sequence_id_graph.add_node(forced_path[node])
						sequence_id_graph.add_edge(prev_id , forced_path[node] , weight=1)
						prev_id = forced_path[node]

		new_list[chr_id] = nx.shortest_path(sequence_id_graph , "BEGIN" , "END")[1:-1]

	# Report warnings for all forced sequences that are not compatible with the marker tiling path
	for chr_id in sorted(forced_list.keys()) :
		unused_seq_ids = []
		for seq_id in forced_list[chr_id] :
			if seq_id not in new_list[chr_id] :
				unused_seq_ids.append(seq_id)
		if not unused_seq_ids == [] :
			print >> sys.stderr , "[WARNING] For chomosome " + chr_id + " sequences were removed from the required list as found incompatible with the marker tiling path. Sequences IDs: " + ", ".join([str(x) for x in unused_seq_ids])
	return new_list


def get_mapped_ids(unique_hits) :
	# unique_hits[Tid] =
	# 	[ 	... ,
	# 		[ Qid , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ] ,
	#		... ]
	Qid_ids = []
	for chr_id in unique_hits:
		Qid_ids.append( [ x[0] for x in unique_hits[chr_id] ] )
	return Qid_ids


def remove_missing_from_list(list, all_mapped , list_type = "") :
	new_list = {}
	for chr_id in list.keys() :
		for element in list[chr_id] :
			if element in all_mapped :
				if chr_id not in new_list :
					new_list[chr_id] = []
				new_list[chr_id].append(element)
			else :
				print >> sys.stderr , "[WARNING] Sequence " + str(element) + " from the list of " + list_type + " has found no mapping position. It will be excluded from list"
	return new_list


def rejected_QC(out_dir, query_name , query_fasta_db, chr_id, fasta_db_1, fasta_db_2, chr_to_fasta_1, chr_to_fasta_2, hap2_on_hap1_coords_file, associated_input_seqid_file , associated_legacy_ids_file, agp_db, legacy_agp, input_agp , groups_by_sequence , marker_bed, marker_usage_db, marker_map, cores, paths) :
	# marker_usage_db --> clean_marker_set_by_seq --> clean_marker_set_by_seq[seq_id]["markers"] : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ]
	# groups_by_sequence[seq_id] = [ ... , group_id , ... ]
	hap1_id = chr_to_fasta_1[chr_id]
	hap2_id = chr_to_fasta_2[chr_id]
	query_id = query_name[:-2]
	query_orientation = query_name[-1]
	query_seq_from_db = query_fasta_db[query_id]
	query_len = len(query_fasta_db[query_id])
	# Generate files for plots
	structure_file = query_id + ".structure.tsv"
	structure_file_fullpath = out_dir + "/" + structure_file
	agp_table = []
	for seq_id in agp_db :
		seq_agp = agp_db[seq_id]
		for start in sorted(seq_agp.keys()) :
			Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = seq_agp[start]
			if Compnt_Type == "W" :
				if CompntId in groups_by_sequence :
					groups = groups_by_sequence[CompntId]
					if len(groups) == 1 :
						group = groups[0]
					else :
						group = "multiple_groups"
				else :
					group = "none"
				agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation , group ])
	structure_file_fullpath = write_table(agp_table, structure_file_fullpath)

	if not legacy_agp == "" :
		legacy_structure_file = query_id + ".legacy_structure.tsv"
		legacy_structure_file_fullpath = out_dir + "/" + legacy_structure_file
		legacy_agp_table = []
		for seq_id in legacy_agp :
			seq_agp = legacy_agp[seq_id]
			for start in sorted(seq_agp.keys()) :
				Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = seq_agp[start]
				if Compnt_Type == "W" :
					if CompntId in groups_by_sequence :
						groups = groups_by_sequence[CompntId]
						if len(groups) == 1 :
							group = groups[0]
						else :
							group = "multiple_groups"
					else :
						group = "none"
					agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation , group ])
					legacy_agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation , group ])
		if not input_agp == "" :
			for seq_id in input_agp :
				seq_agp = input_agp[seq_id]
				for start in sorted(seq_agp.keys()) :
					Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = seq_agp[start]
					if Compnt_Type == "W" :
						if CompntId in groups_by_sequence :
							groups = groups_by_sequence[CompntId]
							if len(groups) == 1 :
								group = groups[0]
							else :
								group = "multiple_groups"
						else :
							group = "none"
						agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation , group ])
						legacy_agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation , group ])
		legacy_structure_file_fullpath = write_table(legacy_agp_table, legacy_structure_file_fullpath)
	else :
		legacy_structure_file = "0"
		# No legacy  structure to load

	if not associated_legacy_ids_file == "" :
		associated_legacy_ids = read_table(associated_legacy_ids_file)
		# associated_legacy_ids = [ ... , [ "hap1_legacy" , Tid , Tstart , Tstop , "Query_legacy" , Qid , Qstart , Qstop , group_id ] , ... ]
		# Coordinates are on pseudomolecules for hap1 and hap2 regions, on HS for "Query_legacy"
		associated_seqid_file = query_id + ".legacy_associations.tsv"
		associated_seqid_file_fullpath = out_dir + "/" + associated_seqid_file
		associated_seqid_file_connection = open( associated_seqid_file_fullpath , 'w')
		# make table file for plotting polygons
		category = polygons_from_ranges(associated_seqid_file_connection ,associated_legacy_ids , 1 , 0 , "hap1_to_hap2" , hap1_id , hap2_id)
		category = polygons_from_ranges(associated_seqid_file_connection ,associated_legacy_ids , 1 , 0 , "hap1_to_unplaced" , hap1_id , query_name)
		category = polygons_from_ranges(associated_seqid_file_connection ,associated_legacy_ids , 0 , 1 , "hap2_to_unplaced" , hap2_id , query_name)
		associated_seqid_file_connection.close()
	else :
		if not associated_input_seqid_file == "" :
			associated_seqid = read_table(associated_input_seqid_file)
			# associated_seqid = [ ... , [ "hap1_HS" , Tid , Tstart , Tstop , "Query_HS" , Qid , Qstart , Qstop , group_id ] , ... ]
			# Coordinates are on pseudomolecules for hap1 and hap2 regions, on HS for "Query_legacy
			associated_seqid_file = query_id + ".associations.tsv"
			associated_seqid_file_fullpath = out_dir + "/" + associated_seqid_file
			associated_seqid_file_connection = open( associated_seqid_file_fullpath , 'w')
			# make table file for plotting polygons
			category = polygons_from_ranges(associated_seqid_file_connection ,associated_seqid , 1 , 0 , "hap1_to_hap2" , hap1_id , hap2_id)
			category = polygons_from_ranges(associated_seqid_file_connection ,associated_seqid , 1 , 0 , "hap1_to_unplaced" , hap1_id , query_name)
			category = polygons_from_ranges(associated_seqid_file_connection ,associated_seqid , 0 , 1 , "hap2_to_unplaced" , hap2_id , query_name)
			associated_seqid_file_connection.close()
		else :
			# No association info based on sequences to use in the plot
			# Map unplaced and analyse all mappings (also hap1 to hap2) to define collinear regions
			hit_target_ranges = []
			associated_seqid_file = query_id + ".hits.tsv"
			# map
			seq_len_db = {}
			hap1_seq = { hap1_id : fasta_db_1[hap1_id] }
			seq_len_db[hap1_id] = len(hap1_seq[hap1_id])
			hap2_seq = {hap2_id : fasta_db_2[hap2_id] }
			seq_len_db[hap2_id] = len(hap2_seq[hap2_id])
			unplaced_seq = {}
			if query_name[-1] == "+" :
				unplaced_seq[query_name] = query_seq_from_db
			else :
				unplaced_seq[query_name] = str(Seq(query_seq_from_db).reverse_complement()).upper()
			seq_len_db[query_name] = len(query_seq_from_db)

			hap1_fasta = out_dir + "/" + hap1_id + ".fasta"
			write_fasta_from_db( hap1_seq , hap1_fasta )
			hap2_fasta = out_dir + "/" + hap2_id + ".fasta"
			write_fasta_from_db( hap2_seq , hap2_fasta )
			query_fasta = out_dir + "/" + query_id + ".fasta"
			write_fasta_from_db( unplaced_seq , query_fasta )

			unplaced_on_hap1_prefix = query_id + ".on.hap1"
			unplaced_on_hap1_prefix_fullpath =  out_dir + "/" + query_id + ".on.hap1"
			unplaced_on_hap1_coords = unplaced_on_hap1_prefix_fullpath + ".coords"
			unplaced_on_hap1_coords = map_nucmer( hap1_fasta , query_fasta , int(cores) , unplaced_on_hap1_coords , paths["nucmer"] , paths["show-coords"] , " --forward " , " -l -r -T -H ")
			unplaced_on_hap1_hits = read_nucmer_coords(unplaced_on_hap1_coords)
			unplaced_on_hap1_best_alignment = hits_best_tiling_path(unplaced_on_hap1_hits, seq_len_db)
			# unplaced_on_target_best_alignment = [ ... , [ Tid , int(Tstart) , int(Tstop) , Qid, int(Qstart) , int(Qstop) , int(align_length) ,  int(match_length) ] , ... ]
			# Convert to cords table format
			for hit in unplaced_on_hap1_best_alignment :
				Tid , Tstart , Tstop , Qid , Qstart , Qstop , align_length ,  match_length = hit
				# hit format: ["hap1_legacy" , Tid , Tstart , Tstop , "Query_legacy" , Qid , Qstart , Qstop]
				hit_target_ranges.append( [ "hap1" , Tid , str(min(int(Tstart) , int(Tstop))) , str(max(int(Tstart) , int(Tstop))) , "Unplace" , Qid , str(min(int(Qstart) , int(Qstop))) , str(max(int(Qstart) , int(Qstop))) , "map" ] )

			unplaced_on_hap2_prefix = query_id + ".on.hap2"
			unplaced_on_hap2_prefix_fullpath =  out_dir + "/" + query_id + ".on.hap2"
			unplaced_on_hap2_coords = unplaced_on_hap2_prefix_fullpath + ".coords"
			unplaced_on_hap2_coords = map_nucmer( hap2_fasta , query_fasta , int(cores) , unplaced_on_hap2_coords , paths["nucmer"] , paths["show-coords"] , " --forward " , " -l -r -T -H ")
			unplaced_on_hap2_hits = read_nucmer_coords(unplaced_on_hap2_coords)
			unplaced_on_hap2_best_alignment = hits_best_tiling_path(unplaced_on_hap2_hits, seq_len_db)
			# unplaced_on_target_best_alignment = [ ... , [ Tid , int(Tstart) , int(Tstop) , Qid, int(Qstart) , int(Qstop) , int(align_length) ,  int(match_length) ] , ... ]
			# Convert to cords table format
			for hit in unplaced_on_hap2_best_alignment :
				Tid , Tstart , Tstop , Qid , Qstart , Qstop , align_length ,  match_length = hit
				# hit format: ["hap2_legacy" , Tid , Tstart , Tstop , "Query_legacy" , Qid , Qstart , Qstop]
				hit_target_ranges.append( [ "hap2" , Tid , str(min(int(Tstart) , int(Tstop))) , str(max(int(Tstart) , int(Tstop))) , "Unplace" , Qid , str(min(int(Qstart) , int(Qstop))) , str(max(int(Qstart) , int(Qstop))) , "map"] )

			hit_target_ranges_file = out_dir + "/" + query_id + ".mappign_hits.txt"
			# Read Hap1 to Hap2 hits
			target_coords = read_table(hap2_on_hap1_coords_file)
			# target_coords_file =[ ... , [ tID , tLen , tStart , tStop , qID , qLen , qStart , qStop , identity , match ] , ... ]
			# filter hits and convert to ranges
			for hit in target_coords :
				tID , tLen , tStart , tStop , qID , qLen , qStart , qStop , identity , match = hit
				if tID == hap1_id and qID == hap2_id :
					hit_target_ranges.append( [ "hap1" , tID , str(min(int(tStart) , int(tStop))) , str(max(int(tStart) , int(tStop))) , "hap2" , qID , str(min(int(qStart) , int(qStop))) , str(max(int(qStart) , int(qStop))) , "map" ] )

			hit_target_ranges_file = write_table( hit_target_ranges , hit_target_ranges_file )

			# generate associations and make table file for plotting polygons
			associated_seqid_file_fullpath = out_dir + "/" + associated_seqid_file
			associated_seqid_file_connection = open( associated_seqid_file_fullpath , 'w')
			category = polygons_from_ranges(associated_seqid_file_connection ,hit_target_ranges , 1 , 0 , "hap1_to_hap2" , hap1_id , hap2_id)
			category = polygons_from_ranges(associated_seqid_file_connection ,hit_target_ranges , 1 , 0 , "hap1_to_unplaced" , hap1_id , query_name)
			category = polygons_from_ranges(associated_seqid_file_connection ,hit_target_ranges , 0 , 1 , "hap2_to_unplaced" , hap2_id , query_name)
			associated_seqid_file_connection.close()

	if not marker_bed == "" :
		# marker_bed == markers_db >> merged seq_id keys for pseudomolecules and input sequences
		# 	markers_db[seq_id][marker_id] = [ ... , [seq_id , start , stop , marker_id] , ... ]
		# marker_usage_db == clean_marker_set_by_seq
		# 	clean_marker_set_by_seq[seq_id] = {
		#		clean_marker_set_by_seq[seq_id]["id"] : seq_id
		# 		clean_marker_set_by_seq[seq_id]["chr"] : chr_id
		# 		clean_marker_set_by_seq[seq_id]["markers"] : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ]
		# 		clean_marker_set_by_seq[seq_id]["range"] : [marker_pos_min , marker_pos_max] ,
		# 		clean_marker_set_by_seq[seq_id]["orientation"] : ["+" or "-" or "."] }
		# marker_map == marker_map_by_seq
		#	marker_map[chr_id] = [ ... , [ int(pos) , marker_id ] , ... ]
		marker_scale = {}
		for marker in marker_map[chr_id] :
			pos, marker_id = marker
			marker_scale[marker_id] = pos

		marker_all_sequence_table_file = query_id + ".marker_all_sequence.tsv"
		marker_all_sequence_table_file_fullpath = out_dir + "/" + marker_all_sequence_table_file
		associated_markers_file = query_id + ".marker_associations.tsv"
		associated_markers_file_fullpath = out_dir + "/" + associated_markers_file
		# marker_all_sequence_table contains
		# 	a) all markers positions according to hap1 ,hap2 and unplaced query coordinates
		#	b) info on being used or not
		#  format associated_markers sets for segment plot >> [ x = t_start, y = t_height, xend = Q_start, yend = q_height, position(for color) , marker_id , category]
		marker_all_sequence_table = []
		associated_markers = []
		used_markers_by_seq = {}
		markers_ranges = {}

		#print >> sys.stderr , "marker_usage_db keys: " + str(marker_usage_db.keys())

		for component in agp_table :
			#print >> sys.stderr , component
			Obj_Name , Obj_start , Obj_End , CompntId , Orientation, group = component
			#print >> sys.stderr , Obj_Name
			if Obj_Name not in [ hap1_id , hap2_id , query_id ] :
				#print >> sys.stderr , ">>> " + Obj_Name + " not listed"
				continue
			else :
				#print >> sys.stderr , Obj_Name + " listed:"
				if Obj_Name == query_id :
					if Obj_Name in marker_usage_db :
						# Do not search for markers in legacy contig, info not available, only pseudomoelcules and input sequences
						if Obj_Name not in used_markers_by_seq :
							used_markers_by_seq[Obj_Name] = []
							markers_ranges[Obj_Name] = []

						if not marker_usage_db[Obj_Name]["orientation"] == "." :
							used_markers_by_seq[Obj_Name] += [str(x[2]) for x in marker_usage_db[Obj_Name]["markers"] ]
							#print >> sys.stderr , "---- Obj_Name: " + Obj_Name + " has direction"
						else :
							unordered_marker_list = [str(x[2]) for x in marker_usage_db[Obj_Name]["markers"]["+"] ]
							unordered_marker_list += [str(x[2]) for x in marker_usage_db[Obj_Name]["markers"]["-"] ]
							used_markers_by_seq[Obj_Name] += list(set(unordered_marker_list))
							#print >> sys.stderr , "---- Obj_Name: " + Obj_Name + " does not have a direction"
						markers_ranges[Obj_Name].append( [ Obj_Name , marker_usage_db[Obj_Name]["range"][0] , marker_usage_db[Obj_Name]["range"][1] , marker_usage_db[Obj_Name]["orientation"] , marker_usage_db[Obj_Name]["id"] , Obj_Name ] )
					#else :
					#	print >> sys.stderr , "---- CompntId: " + Obj_Name + " not in marker_usage_db"
				else :
					if CompntId in marker_usage_db :
						# Do not search for markers in legacy contig, info not available, only pseudomoelcules and input sequences
						if Obj_Name not in used_markers_by_seq :
							used_markers_by_seq[Obj_Name] = []
							markers_ranges[Obj_Name] = []

						if not marker_usage_db[CompntId]["orientation"] == "." :
							used_markers_by_seq[Obj_Name] += [str(x[2]) for x in marker_usage_db[CompntId]["markers"] ]
							#print >> sys.stderr , "---- CompntId: " + CompntId + " has direction"
						else :
							unordered_marker_list = [str(x[2]) for x in marker_usage_db[CompntId]["markers"]["+"] ]
							unordered_marker_list += [str(x[2]) for x in marker_usage_db[CompntId]["markers"]["-"] ]
							used_markers_by_seq[Obj_Name] += list(set(unordered_marker_list))
							#print >> sys.stderr , "---- CompntId: " + CompntId + " does not have a direction"
						markers_ranges[Obj_Name].append( [ Obj_Name , marker_usage_db[CompntId]["range"][0] , marker_usage_db[CompntId]["range"][1] , marker_usage_db[CompntId]["orientation"] , marker_usage_db[CompntId]["id"] , Obj_Name ] )
					#else :
					#	print >> sys.stderr , "---- CompntId: " + CompntId + " not in marker_usage_db"

		#print >> sys.stderr , "used_markers_by_seq keys: " + str(used_markers_by_seq.keys())
		# generate table of marker ranges -> coordinates on marker map -> pick the used markers regions in clean_marker_set_by_seq
		markers_ranges_file = query_id + ".used_markers_range.tsv"
		markers_ranges_file_fullpath = out_dir + "/" + markers_ranges_file
		markers_ranges_file_fullpath_connection = open( markers_ranges_file_fullpath , 'w')
		for chr_id in markers_ranges.keys() :
			for range in markers_ranges[chr_id] :
				print >> markers_ranges_file_fullpath_connection , "\t".join([ str(x) for x in range])
		markers_ranges_file_fullpath_connection.close()

		# Table format: seq_id , start ,stop , orientation , component_id , group
		#  format associated_markers sets for segment plot >> [ x = t_start, y = t_height, xend = Q_start, yend = q_height, position(for color) , marker_id , category]
		if hap1_id in marker_bed :
			for marker_id in sorted(marker_bed[hap1_id].keys()) :
				for hit in marker_bed[hap1_id][marker_id] :
					chr_id , start , stop , marker_name = hit
					# Add to marker_all_sequence_table
					if marker_id in marker_scale :
						marker_pos = marker_scale[marker_id]
						if hap1_id in used_markers_by_seq :
							if marker_id in used_markers_by_seq[hap1_id] :
								usage = "TRUE"
							else :
								usage = "FALSE"
						else :
							usage = "TRUE"
						marker_all_sequence_table.append([chr_id , start , marker_id , marker_pos, usage])
						# Find matches with unplaced
						if (query_id in marker_bed )  and (marker_id in marker_bed[query_id]) :
							for hit2 in marker_bed[query_id][marker_id] :
								un_chr_id , un_start , un_stop , un_marker_id = hit2
								associated_markers.append( [ start , 1 , un_start , 0 , marker_pos , marker_id , "hap1_to_unplaced" ] )
						# Find matches with hap2
						if marker_id in marker_bed[hap2_id] :
							for hit3 in marker_bed[hap2_id][marker_id] :
								hap2_chr_id , hap2_start , hap2_stop , hap2_marker_id = hit3
								associated_markers.append( [ start , 1 , hap2_start , 0 , marker_pos , marker_id , "hap1_to_hap2" ] )
					else :
						continue

		if query_id in marker_bed :
			# query_id has markers on it
			for marker_id in sorted(marker_bed[query_id].keys()) :
				for hit in marker_bed[query_id][marker_id] :
					chr_id , start , stop , marker_id = hit
					if marker_id in marker_scale :
						marker_pos = marker_scale[marker_id]
						if query_id in used_markers_by_seq :
							if marker_id in used_markers_by_seq[query_id] :
								usage = "TRUE"
							else :
								usage = "FALSE"
						else :
							usage = "FALSE"

						# Translate coordinates >> if "-" direction
						if query_orientation == "-" :
							start = str(int(query_len) - int(start))
						# Add to marker_all_sequence_table
						marker_all_sequence_table.append([query_id , start , marker_id , marker_pos , usage])
						# Find matches with hap2
						if marker_id in marker_bed[hap2_id] :
							for hit2 in marker_bed[hap2_id][marker_id] :
								hap2_chr_id , hap2_start , hap2_stop , hap2_marker_id = hit2
								associated_markers.append( [ start , 1 , hap2_start , 0 , marker_pos , marker_id , "hap2_to_unplaced" ] )
					else :
						continue

		if hap2_id in marker_bed :
			for marker_id in sorted(marker_bed[hap2_id].keys()) :
				for hit in marker_bed[hap2_id][marker_id] :
					chr_id , start , stop , marker_id = hit
					if marker_id in marker_scale :
						marker_pos = marker_scale[marker_id]
						if hap2_id in used_markers_by_seq :
							if marker_id in used_markers_by_seq[hap2_id] :
								usage = "TRUE"
							else :
								usage = "FALSE"
						else :
							usage = "TRUE"
						marker_all_sequence_table.append([chr_id , start , marker_id , marker_pos, usage])
					else :
						continue

		marker_all_sequence_table_file_fullpath = write_table(marker_all_sequence_table , marker_all_sequence_table_file_fullpath)
		associated_markers_file_fullpath = write_table(associated_markers , associated_markers_file_fullpath)

	else :
		marker_all_sequence_table_file = "0"
		associated_markers_file = "0"
		markers_ranges_file = "0"


	# Select the necessary plot Rmd file to render
	# 	all scripts share the same input structure, missing/unused elements are substituted with 0
	if not marker_bed == "" :
		if not legacy_structure_file == "0" :
			script=scriptDirectory + "/unplaced_qc.Rmd"
		else :
			script=scriptDirectory + "/unplaced_qc.no_legacy.Rmd"
	else :
		if not legacy_structure_file == "0" :
			script=scriptDirectory + "/unplaced_qc.no_markers.Rmd"
		else :
			script=scriptDirectory + "/unplaced_qc.no_markers.no_legacy.Rmd"

	## TOD0: Perform the rendering for each rejected sequence
	out_file_name_prefix = query_id + "_qc"
	output_file = out_dir + "/" + out_file_name_prefix
	log_connection = open( out_dir + "/." + query_id + "_qc.log" , 'w')
	err_connection = open( out_dir + "/." + query_id + "_qc.err" , 'w')
	command = "Rscript -e 'library(rmarkdown) ; rmarkdown::render(\"" + os.path.realpath(script) + "\" , knit_root_dir = \"" + os.path.realpath(out_dir) + "\" , output_file = \"" + out_file_name_prefix + "\" , output_dir = \"" + os.path.realpath(out_dir) + "\" , params=list( filename = \"" + os.path.realpath(output_file) + "\" , Hap1= \"" + hap1_id + "\" , Hap2= \"" + hap2_id + "\" , unplacedID= \"" + query_id + "\" , structure = \"" + structure_file + "\" , legacy = \"" + legacy_structure_file + "\" , markers = \"" + marker_all_sequence_table_file + "\" , seq_relationships = \"" + associated_seqid_file + "\" , marker_relationship= \"" + associated_markers_file + "\" , markers_ranges= \"" + markers_ranges_file +"\"))'"
	# Rscript -e 'library("rmarkdown") ; 		rmarkdown::render( "unplaced_qc.Rmd" ,                   knit_root_dir = ""                                    , output_file = "test"                           , output_dir = ""                                    , params=list( filename = "test"                                    , Hap1=   "NEW_Hap1_chr10"  , Hap2=   "NEW_Hap2_chr10"  , unplacedID= "seq99"            , structure = "seq99.structure.tsv" ,      legacy = "seq99.legacy_structure.tsv"      , markers = "seq99.marker_all_sequence.tsv"            , seq_relationships = "seq99.legacy_associations.tsv"   , marker_relationship = "seq99.marker_associations.tsv"    , markers_ranges = "seq99.used_markers_range.tsv"))'
	print >> sys.stderr, "#### Running command: " + command
	reportProcess = subprocess.Popen( command , shell=True , stdout=log_connection , stderr=err_connection )
	output, error = reportProcess.communicate()
	log_connection.close()
	err_connection.close()
	outfiles = {}
	outfiles["html"] = query_id + "_qc.html"
	outfiles["pdf"] = query_id + "_qc.pdf"
	outfiles["png"] = query_id + "_qc.png"
	outfiles["size"] = query_len
	# Index must be produced inside out_dir
	return outfiles


def compare_structures(out_dir, target_name , query_name , prefix , fasta_db_1, fasta_db_2, coords_file, associated_input_seqid_file , associated_legacy_ids_file, agp_db, legacy_agp, input_agp , groups_by_sequence , marker_bed, marker_usage_db, marker_map, cores, paths) :
	## TODO: Update to periwise comparion only
	## marker_usage_db --> clean_marker_set_by_seq --> clean_marker_set_by_seq[seq_id]["markers"] : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ]
	## groups_by_sequence[seq_id] = [ ... , group_id , ... ]

	## Generate files for plots
	#structure_file = prefix + ".structure.tsv"
	#structure_file_fullpath = out_dir + "/" + structure_file
	#agp_table = []
	#for seq_id in agp_db :
	#	seq_agp = agp_db[seq_id]
	#	for start in sorted(seq_agp.keys()) :
	#		Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = seq_agp[start]
	#		if Compnt_Type == "W" :
	#			if CompntId in groups_by_sequence :
	#				groups = groups_by_sequence[CompntId]
	#				if len(groups) == 1 :
	#					group = groups[0]
	#				else :
	#					group = "multiple_groups"
	#			else :
	#				group = "none"
	#			agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation , group ])
	#structure_file_fullpath = write_table(agp_table, structure_file_fullpath)

	#if not legacy_agp == "" :
	#	legacy_structure_file = prefix + ".legacy_structure.tsv"
	#	legacy_structure_file_fullpath = out_dir + "/" + legacy_structure_file
	#	legacy_agp_table = []
	#	for seq_id in legacy_agp :
	#		seq_agp = legacy_agp[seq_id]
	#		for start in sorted(seq_agp.keys()) :
	#			Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = seq_agp[start]
	#			if Compnt_Type == "W" :
	#				if CompntId in groups_by_sequence :
	#					groups = groups_by_sequence[CompntId]
	#					if len(groups) == 1 :
	#						group = groups[0]
	#					else :
	#						group = "multiple_groups"
	#				else :
	#					group = "none"
	#				agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation , group ])
	#				legacy_agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation , group ])
	#	if not input_agp == "" :
	#		for seq_id in input_agp :
	#			seq_agp = input_agp[seq_id]
	#			for start in sorted(seq_agp.keys()) :
	#				Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = seq_agp[start]
	#				if Compnt_Type == "W" :
	#					if CompntId in groups_by_sequence :
	#						groups = groups_by_sequence[CompntId]
	#						if len(groups) == 1 :
	#							group = groups[0]
	#						else :
	#							group = "multiple_groups"
	#					else :
	#						group = "none"
	#					agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation , group ])
	#					legacy_agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation , group ])
	#	legacy_structure_file_fullpath = write_table(legacy_agp_table, legacy_structure_file_fullpath)
	#else :
	#	legacy_structure_file = "0"
	#	# No legacy  structure to load

	#if not associated_legacy_ids_file == "" :
	#	associated_legacy_ids = read_table(associated_legacy_ids_file)
	#	# associated_legacy_ids = [ ... , [ "hap1_legacy" , Tid , Tstart , Tstop , "Query_legacy" , Qid , Qstart , Qstop , group_id ] , ... ]
	#	# Coordinates are on pseudomolecules for hap1 and hap2 regions, on HS for "Query_legacy"
	#	associated_seqid_file = prefix + ".legacy_associations.tsv"
	#	associated_seqid_file_fullpath = out_dir + "/" + associated_seqid_file
	#	associated_seqid_file_connection = open( associated_seqid_file_fullpath , 'w')
	#	# make table file for plotting polygons
	#	category = polygons_from_ranges(associated_seqid_file_connection ,associated_legacy_ids , 1 , 0 , "hap1_to_hap2" , target_name , query_name)
	#	category = polygons_from_ranges(associated_seqid_file_connection ,associated_legacy_ids , 1 , 0 , "hap1_to_unplaced" , target_name , query_name)
	#	category = polygons_from_ranges(associated_seqid_file_connection ,associated_legacy_ids , 0 , 1 , "hap2_to_unplaced" , query_name , query_name)
	#	associated_seqid_file_connection.close()
	#else :
	#	if not associated_input_seqid_file == "" :
	#		associated_seqid = read_table(associated_input_seqid_file)
	#		# associated_seqid = [ ... , [ "hap1_HS" , Tid , Tstart , Tstop , "Query_HS" , Qid , Qstart , Qstop , group_id ] , ... ]
	#		# Coordinates are on pseudomolecules for hap1 and hap2 regions, on HS for "Query_legacy
	#		associated_seqid_file = prefix + ".associations.tsv"
	#		associated_seqid_file_fullpath = out_dir + "/" + associated_seqid_file
	#		associated_seqid_file_connection = open( associated_seqid_file_fullpath , 'w')
	#		# make table file for plotting polygons
	#		category = polygons_from_ranges(associated_seqid_file_connection ,associated_seqid , 1 , 0 , "hap1_to_hap2" , target_name , query_name)
	#		category = polygons_from_ranges(associated_seqid_file_connection ,associated_seqid , 1 , 0 , "hap1_to_unplaced" , target_name , query_name)
	#		category = polygons_from_ranges(associated_seqid_file_connection ,associated_seqid , 0 , 1 , "hap2_to_unplaced" , query_name , query_name)
	#		associated_seqid_file_connection.close()
	#	else :
	#		# No association info based on sequences to use in the plot
	#		# Map unplaced and analyse all mappings (also hap1 to hap2) to define collinear regions
	#		hit_target_ranges = []
	#		associated_seqid_file = prefix + ".hits.tsv"
	#		# map
	#		seq_len_db = {}
	#		hap1_seq = { target_name : fasta_db_1[target_name] }
	#		seq_len_db[target_name] = len(hap1_seq[target_name])
	#		hap2_seq = {query_name : fasta_db_2[query_name] }
	#		seq_len_db[query_name] = len(hap2_seq[query_name])
	#		unplaced_seq = {}
	#		if query_name[-1] == "+" :
	#			unplaced_seq[query_name] = query_seq_from_db
	#		else :
	#			unplaced_seq[query_name] = str(Seq(query_seq_from_db).reverse_complement()).upper()
	#		seq_len_db[query_name] = len(query_seq_from_db)

	#		hap1_fasta = out_dir + "/" + target_name + ".fasta"
	#		write_fasta_from_db( hap1_seq , hap1_fasta )
	#		hap2_fasta = out_dir + "/" + query_name + ".fasta"
	#		write_fasta_from_db( hap2_seq , hap2_fasta )

	#		unplaced_on_hap1_prefix = prefix
	#		unplaced_on_hap1_prefix_fullpath =  out_dir + "/" + prefix
	#		unplaced_on_hap1_coords = unplaced_on_hap1_prefix_fullpath + ".coords"
	#		unplaced_on_hap1_coords = map_nucmer( hap1_fasta , query_fasta , int(cores) , unplaced_on_hap1_coords , paths["nucmer"] , paths["show-coords"] , " --forward " , " -l -r -T -H ")
	#		unplaced_on_hap1_hits = read_nucmer_coords(unplaced_on_hap1_coords)
	#		unplaced_on_hap1_best_alignment = hits_best_tiling_path(unplaced_on_hap1_hits, seq_len_db)
	#		# unplaced_on_target_best_alignment = [ ... , [ Tid , int(Tstart) , int(Tstop) , Qid, int(Qstart) , int(Qstop) , int(align_length) ,  int(match_length) ] , ... ]
	#		# Convert to cords table format
	#		for hit in unplaced_on_hap1_best_alignment :
	#			Tid , Tstart , Tstop , Qid , Qstart , Qstop , align_length ,  match_length = hit
	#			# hit format: ["hap1_legacy" , Tid , Tstart , Tstop , "Query_legacy" , Qid , Qstart , Qstop]
	#			hit_target_ranges.append( [ "hap1" , Tid , str(min(int(Tstart) , int(Tstop))) , str(max(int(Tstart) , int(Tstop))) , "Unplace" , Qid , str(min(int(Qstart) , int(Qstop))) , str(max(int(Qstart) , int(Qstop))) , "map" ] )

	#		hit_target_ranges_file = out_dir + "/" + prefix + ".mappign_hits.txt"
	#		# Read Hap1 to Hap2 hits
	#		target_coords = read_table(coords_file)
	#		# target_coords_file =[ ... , [ tID , tLen , tStart , tStop , qID , qLen , qStart , qStop , identity , match ] , ... ]
	#		# filter hits and convert to ranges
	#		for hit in target_coords :
	#			tID , tLen , tStart , tStop , qID , qLen , qStart , qStop , identity , match = hit
	#			if tID == target_name and qID == query_name :
	#				hit_target_ranges.append( [ "hap1" , tID , str(min(int(tStart) , int(tStop))) , str(max(int(tStart) , int(tStop))) , "hap2" , qID , str(min(int(qStart) , int(qStop))) , str(max(int(qStart) , int(qStop))) , "map" ] )

	#		hit_target_ranges_file = write_table( hit_target_ranges , hit_target_ranges_file )

	#		# generate associations and make table file for plotting polygons
	#		associated_seqid_file_fullpath = out_dir + "/" + associated_seqid_file
	#		associated_seqid_file_connection = open( associated_seqid_file_fullpath , 'w')
	#		category = polygons_from_ranges(associated_seqid_file_connection ,hit_target_ranges , 1 , 0 , "hap1_to_hap2" , target_name , query_name)
	#		category = polygons_from_ranges(associated_seqid_file_connection ,hit_target_ranges , 1 , 0 , "hap1_to_unplaced" , target_name , query_name)
	#		category = polygons_from_ranges(associated_seqid_file_connection ,hit_target_ranges , 0 , 1 , "hap2_to_unplaced" , query_name , query_name)
	#		associated_seqid_file_connection.close()

	#if not marker_bed == "" :
	#	# marker_bed == markers_db >> merged seq_id keys for pseudomolecules and input sequences
	#	# 	markers_db[seq_id][marker_id] = [ ... , [seq_id , start , stop , marker_id] , ... ]
	#	# marker_usage_db == clean_marker_set_by_seq
	#	# 	clean_marker_set_by_seq[seq_id] = {
	#	#		clean_marker_set_by_seq[seq_id]["id"] : seq_id
	#	# 		clean_marker_set_by_seq[seq_id]["chr"] : chr_id
	#	# 		clean_marker_set_by_seq[seq_id]["markers"] : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ]
	#	# 		clean_marker_set_by_seq[seq_id]["range"] : [marker_pos_min , marker_pos_max] ,
	#	# 		clean_marker_set_by_seq[seq_id]["orientation"] : ["+" or "-" or "."] }
	#	# marker_map == marker_map_by_seq
	#	#	marker_map[chr_id] = [ ... , [ int(pos) , marker_id ] , ... ]
	#	marker_scale = {}
	#	for chr_id in sorted(marker_map.keys()) :
	#		for marker in marker_map[chr_id] :
	#			pos, marker_id = marker
	#			marker_scale[marker_id] = pos

	#	marker_all_sequence_table_file = prefix + ".marker_all_sequence.tsv"
	#	marker_all_sequence_table_file_fullpath = out_dir + "/" + marker_all_sequence_table_file
	#	associated_markers_file = prefix + ".marker_associations.tsv"
	#	associated_markers_file_fullpath = out_dir + "/" + associated_markers_file
	#	# marker_all_sequence_table contains
	#	# 	a) all markers positions according to hap1 ,hap2 and unplaced query coordinates
	#	#	b) info on being used or not
	#	#  format associated_markers sets for segment plot >> [ x = t_start, y = t_height, xend = Q_start, yend = q_height, position(for color) , marker_id , category]
	#	marker_all_sequence_table = []
	#	associated_markers = []
	#	used_markers_by_seq = {}
	#	markers_ranges = {}

	#	for component in agp_table :
	#		#print >> sys.stderr , component
	#		Obj_Name , Obj_start , Obj_End , CompntId , Orientation, group = component
	#		#print >> sys.stderr , Obj_Name
	#		if Obj_Name not in [ target_name , query_name ] :
	#			#print >> sys.stderr , ">>> " + Obj_Name + " not listed"
	#			continue
	#		else :
	#			#print >> sys.stderr , Obj_Name + " listed:"
	#			if Obj_Name == query_id :
	#				if Obj_Name in marker_usage_db :
	#					# Do not search for markers in legacy contig, info not available, only pseudomoelcules and input sequences
	#					if Obj_Name not in used_markers_by_seq :
	#						used_markers_by_seq[Obj_Name] = []
	#						markers_ranges[Obj_Name] = []

	#					if not marker_usage_db[Obj_Name]["orientation"] == "." :
	#						used_markers_by_seq[Obj_Name] += [str(x[2]) for x in marker_usage_db[Obj_Name]["markers"] ]
	#						#print >> sys.stderr , "---- Obj_Name: " + Obj_Name + " has direction"
	#					else :
	#						unordered_marker_list = [str(x[2]) for x in marker_usage_db[Obj_Name]["markers"]["+"] ]
	#						unordered_marker_list += [str(x[2]) for x in marker_usage_db[Obj_Name]["markers"]["-"] ]
	#						used_markers_by_seq[Obj_Name] += list(set(unordered_marker_list))
	#						#print >> sys.stderr , "---- Obj_Name: " + Obj_Name + " does not have a direction"
	#					markers_ranges[Obj_Name].append( [ Obj_Name , marker_usage_db[Obj_Name]["range"][0] , marker_usage_db[Obj_Name]["range"][1] , marker_usage_db[Obj_Name]["orientation"] , marker_usage_db[Obj_Name]["id"] , Obj_Name ] )
	#				#else :
	#				#	print >> sys.stderr , "---- CompntId: " + Obj_Name + " not in marker_usage_db"
	#			else :
	#				if CompntId in marker_usage_db :
	#					# Do not search for markers in legacy contig, info not available, only pseudomoelcules and input sequences
	#					if Obj_Name not in used_markers_by_seq :
	#						used_markers_by_seq[Obj_Name] = []
	#						markers_ranges[Obj_Name] = []

	#					if not marker_usage_db[CompntId]["orientation"] == "." :
	#						used_markers_by_seq[Obj_Name] += [str(x[2]) for x in marker_usage_db[CompntId]["markers"] ]
	#						#print >> sys.stderr , "---- CompntId: " + CompntId + " has direction"
	#					else :
	#						unordered_marker_list = [str(x[2]) for x in marker_usage_db[CompntId]["markers"]["+"] ]
	#						unordered_marker_list += [str(x[2]) for x in marker_usage_db[CompntId]["markers"]["-"] ]
	#						used_markers_by_seq[Obj_Name] += list(set(unordered_marker_list))
	#						#print >> sys.stderr , "---- CompntId: " + CompntId + " does not have a direction"
	#					markers_ranges[Obj_Name].append( [ Obj_Name , marker_usage_db[CompntId]["range"][0] , marker_usage_db[CompntId]["range"][1] , marker_usage_db[CompntId]["orientation"] , marker_usage_db[CompntId]["id"] , Obj_Name ] )
	#				#else :
	#				#	print >> sys.stderr , "---- CompntId: " + CompntId + " not in marker_usage_db"
	#	# generate table of marker ranges -> coordinates on marker map -> pick the used markers regions in clean_marker_set_by_seq
	#	markers_ranges_file = query_id + ".used_markers_range.tsv"
	#	markers_ranges_file_fullpath = out_dir + "/" + markers_ranges_file
	#	markers_ranges_file_fullpath_connection = open( markers_ranges_file_fullpath , 'w')
	#	for chr_id in markers_ranges.keys() :
	#		for range in markers_ranges[chr_id] :
	#			print >> markers_ranges_file_fullpath_connection , "\t".join([ str(x) for x in range])
	#	markers_ranges_file_fullpath_connection.close()

	#	# Table format: seq_id , start ,stop , orientation , component_id , group
	#	#  format associated_markers sets for segment plot >> [ x = t_start, y = t_height, xend = Q_start, yend = q_height, position(for color) , marker_id , category]
	#	for marker_id in sorted(marker_bed[target_name].keys()) :
	#		for hit in marker_bed[target_name][marker_id] :
	#			chr_id , start , stop , marker_name = hit
	#			# Add to marker_all_sequence_table
	#			if marker_id in marker_scale :
	#				marker_pos = marker_scale[marker_id]
	#				if marker_id in used_markers_by_seq[target_name] :
	#					usage = "TRUE"
	#				else :
	#					usage = "FALSE"
	#				marker_all_sequence_table.append([chr_id , start , marker_id , marker_pos, usage])
	#				# Find matches with unplaced
	#				if marker_id in marker_bed[query_id] :
	#					for hit2 in marker_bed[query_id][marker_id] :
	#						un_chr_id , un_start , un_stop , un_marker_id = hit2
	#						associated_markers.append( [ start , 1 , un_start , 0 , marker_pos , marker_id , "hap1_to_unplaced" ] )
	#				# Find matches with hap2
	#				if marker_id in marker_bed[query_name] :
	#					for hit3 in marker_bed[query_name][marker_id] :
	#						hap2_chr_id , hap2_start , hap2_stop , hap2_marker_id = hit3
	#						associated_markers.append( [ start , 1 , hap2_start , 0 , marker_pos , marker_id , "hap1_to_hap2" ] )
	#			else :
	#				continue

	#	for marker_id in sorted(marker_bed[query_id].keys()) :
	#		for hit in marker_bed[query_id][marker_id] :
	#			chr_id , start , stop , marker_id = hit
	#			if marker_id in marker_scale :
	#				marker_pos = marker_scale[marker_id]
	#				try :
	#					a = used_markers_by_seq[query_id]
	#				except :
	#					print >> sys.stderr ,  "used_markers_by_seq.keys()"
	#					print >> sys.stderr , sorted(used_markers_by_seq.keys())
	#					sys.exit(33)
	#				if marker_id in used_markers_by_seq[query_id] :
	#					usage = "TRUE"
	#				else :
	#					usage = "FALSE"
	#				# Translate coordinates >> if "-" direction
	#				if query_orientation == "-" :
	#					start = str(int(query_len) - int(start))
	#				# Add to marker_all_sequence_table
	#				marker_all_sequence_table.append([query_id , start , marker_id , marker_pos , usage])
	#				# Find matches with hap2
	#				if marker_id in marker_bed[query_name] :
	#					for hit2 in marker_bed[query_name][marker_id] :
	#						hap2_chr_id , hap2_start , hap2_stop , hap2_marker_id = hit2
	#						associated_markers.append( [ start , 1 , hap2_start , 0 , marker_pos , marker_id , "hap2_to_unplaced" ] )
	#			else :
	#				continue

	#	for marker_id in sorted(marker_bed[query_name].keys()) :
	#		for hit in marker_bed[query_name][marker_id] :
	#			chr_id , start , stop , marker_id = hit
	#			if marker_id in marker_scale :
	#				marker_pos = marker_scale[marker_id]
	#				if marker_id in used_markers_by_seq[query_name] :
	#					usage = "TRUE"
	#				else :
	#					usage = "FALSE"
	#				marker_all_sequence_table.append([chr_id , start , marker_id , marker_pos, usage])
	#			else :
	#				continue

	#	marker_all_sequence_table_file_fullpath = write_table(marker_all_sequence_table , marker_all_sequence_table_file_fullpath)
	#	associated_markers_file_fullpath = write_table(associated_markers , associated_markers_file_fullpath)

	#else :
	#	marker_all_sequence_table_file = "0"
	#	associated_markers_file = "0"
	#	markers_ranges_file = "0"


	## Select the necessary plot Rmd file to render
	## 	all scripts share the same input structure, missing/unused elements are substituted with 0
	## TODO
	#if not marker_bed == "" :
	#	if not legacy_structure_file == "0" :
	#		script=scriptDirectory + "/structure_comparison.Rmd"
	#	else :
	#		script=scriptDirectory + "/structure_comparison.no_legacy.Rmd"
	#else :
	#	if not legacy_structure_file == "0" :
	#		script=scriptDirectory + "/structure_comparison.no_markers.Rmd"
	#	else :
	#		script=scriptDirectory + "/structure_comparison.no_markers.no_legacy.Rmd"

	### TOD0: Perform the rendering for each rejected sequence
	out_file_name_prefix = prefix + "_structure_comparison"
	output_file = out_dir + "/" + out_file_name_prefix
	#log_connection = open( out_dir + "/." + prefix + "_structure_comparison.log" , 'w')
	#err_connection = open( out_dir + "/." + prefix + "_structure_comparison.err" , 'w')
	#command = "Rscript -e 'library(rmarkdown) ; rmarkdown::render(\"" + os.path.realpath(script) + "\" , knit_root_dir = \"" + os.path.realpath(out_dir) + "\" , output_file = \"" + out_file_name_prefix + "\" , output_dir = \"" + os.path.realpath(out_dir) + "\" , params=list( filename = \"" + os.path.realpath(output_file) + "\" , Hap1= \"" + target_name + "\" , Hap2= \"" + query_name + "\" , unplacedID= \"" + query_id + "\" , structure = \"" + structure_file + "\" , legacy = \"" + legacy_structure_file + "\" , markers = \"" + marker_all_sequence_table_file + "\" , seq_relationships = \"" + associated_seqid_file + "\" , marker_relationship= \"" + associated_markers_file + "\" , markers_ranges= \"" + markers_ranges_file +"\"))'"
	## Rscript -e 'library("rmarkdown") ; 		rmarkdown::render( "unplaced_structure_comparison.Rmd" ,                   knit_root_dir = ""                                    , output_file = "test"                           , output_dir = ""                                    , params=list( filename = "test"                                    , Hap1=   "NEW_Hap1_chr10"  , Hap2=   "NEW_Hap2_chr10"  , unplacedID= "seq99"            , structure = "seq99.structure.tsv" ,      legacy = "seq99.legacy_structure.tsv"      , markers = "seq99.marker_all_sequence.tsv"            , seq_relationships = "seq99.legacy_associations.tsv"   , marker_relationship = "seq99.marker_associations.tsv"    , markers_ranges = "seq99.used_markers_range.tsv"))'
	#print >> sys.stderr, "#### Running command: " + command
	#reportProcess = subprocess.Popen( command , shell=True , stdout=log_connection , stderr=err_connection )
	#output, error = reportProcess.communicate()
	#log_connection.close()
	#err_connection.close()
	outfiles = {}
	outfiles["html"] = prefix + "_structure_comparison.html"
	outfiles["pdf"] = prefix + "_structure_comparison.pdf"
	outfiles["png"] = prefix + "_structure_comparison.png"
	# Index must be produced inside out_dir
	return outfiles


def polygons_from_ranges( file_connection , list , target_height , query_height , category , filter_target = "" , filter_query = "") :
	# Points format:
	# [ ... , [ x , y , group_id , category  ] , ...]
	for element in list :
		Tcat , Tid , Tstart , Tstop , Qcat , Qid , Qstart , Qstop , group_id = element
		if not filter_target == "" :
			if not filter_target == Tid :
				continue
		if not filter_query == "" :
			if not filter_query == Qid :
				continue
		print >> file_connection , "\t".join([str(x) for x in [ Tstart , target_height , group_id , category ] ])
		print >> file_connection , "\t".join([str(x) for x in [ Tstop  , target_height , group_id , category ] ])
		print >> file_connection , "\t".join([str(x) for x in [ Qstop  , query_height  , group_id , category ] ])
		print >> file_connection , "\t".join([str(x) for x in [ Qstart , query_height  , group_id , category ] ])
		print >> file_connection , "\t".join([str(x) for x in [ Tstart , target_height , group_id , category ] ])

	return category


def make_seq_pair_from_groups( matching_regions_file , group_file , agp_structure_db , hap1_list , hap2_list , input_list , component_len_db , agp_origin ) :
	matching_regions_file_connection = open(matching_regions_file , "w")
	components_positions = {}
	# agp_structure_db[seq_id][start] = Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation
	# with legacy_agp: seq_id in chr_id and input_seq_id and CompntId in legacy_seq_id >> component_len_db empty so no CompntId_to_CompntId relationship is produced
	# with pseudomolecules_agp_db : seq_id in chr_id and input_seq_id and CompntId in input_seq_id >> component_len_db has input_seq_id lengths so input_seq_id_to_input_seq_id relationship are produced for use on unplaced
	for chr_id in agp_structure_db.keys() :
		for start in agp_structure_db[chr_id].keys() :
			#print >> sys.stderr , agp_structure_db[chr_id][start]
			Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = agp_structure_db[chr_id][start]
			if CompntId not in components_positions :
				components_positions[CompntId] = []
				if CompntId in component_len_db :
					CompntLen = component_len_db[CompntId]
					components_positions[CompntId].append( [ CompntId , "1" , str(CompntLen) ] )
			components_positions[CompntId].append([Obj_Name , Obj_start , Obj_End , CompntId] )
			#print >> sys.stderr, "Adding: " + str([Obj_Name , Obj_start , Obj_End , CompntId])
			# components_positions[CompntId_1] = [ ... , [Obj_Name , Obj_start , Obj_End , CompntId] , ... ]

	# Identify all sequences belonging to the same group
	group_table = read_table(group_file)
	# group_table = [ ... , [ CompntId_1 , group_a ] , [ CompntId_1 , group_b ] , ... ]
	sequence_groups = {}
	groups_by_sequence = {}
	for relationship in group_table :
		seq_id , group_id = relationship
		if group_id not in sequence_groups :
			sequence_groups[group_id] = []
		sequence_groups[group_id].append(seq_id)
		if seq_id not in groups_by_sequence :
			groups_by_sequence[seq_id] = []
		groups_by_sequence[seq_id].append(group_id)

		# sequence_groups[group_a] = [ ... , CompntId_1 , ... ]

	# parse each group and create a match for each sequence pair of sequences for each region they hit + add orientation to sequences in input_list
	for group_id in sequence_groups.keys() :
		#print >> sys.stderr , "sequence_groups[group_id] :"
		#print >> sys.stderr , sequence_groups[group_id]
		grouped_pairs = list(itertools.product(sequence_groups[group_id], sequence_groups[group_id]))
		for pair in grouped_pairs :
			#print >> sys.stderr , "pair :"
			#print >> sys.stderr , pair
			component_1 , component_2 = pair
			if component_1 == component_2 :
				continue
			else :
				component_1_regions = components_positions[component_1]
				match_1_list = []
				# component_1_regions = [ ... , [ seq_1_id , component_1_start , component_1_stop , component_1] , ... ]
				component_2_regions = components_positions[component_2]
				match_2_list = []
				for region_1 in component_1_regions :
					seq_1_id , component_1_start , component_1_stop , component_1_id = region_1
					if seq_1_id in hap1_list :
						component_1_category = "hap1_" + agp_origin
						match_1_list.append( [ component_1_category , seq_1_id , component_1_start , component_1_stop ] )
					elif seq_1_id in hap2_list :
						component_1_category = "hap2_" + agp_origin
						match_1_list.append( [ component_1_category , seq_1_id , component_1_start , component_1_stop ] )
					elif seq_1_id in input_list :
						component_1_category = "input_" + agp_origin
						match_1_list.append( [ component_1_category , seq_1_id , component_1_start , component_1_stop ] )
						match_1_list.append( [ component_1_category , seq_1_id + "|+" , component_1_start , component_1_stop ] )
						match_1_list.append( [ component_1_category , seq_1_id + "|-" , component_1_start , component_1_stop ] )
						match_1_list.append( [ component_1_category , seq_1_id + "|." , component_1_start , component_1_stop ] )
					else :
						print >> sys.stdout , "[WARNING] Sequence " + seq_1_id + " assigned to group " + group_id + " is unknown, ignored"
				for region_2 in component_2_regions :
					seq_2_id , component_2_start , component_2_stop , component_2_id = region_2
					if seq_2_id in hap1_list :
						component_2_category = "hap1_" + agp_origin
						match_2_list.append( [ component_2_category , seq_2_id , component_2_start , component_2_stop ] )
					elif seq_2_id in hap2_list :
						component_2_category = "hap2_" + agp_origin
						match_2_list.append( [ component_2_category , seq_2_id , component_2_start , component_2_stop ] )
					elif seq_2_id in input_list :
						component_2_category = "input_" + agp_origin
						match_2_list.append( [ component_2_category , seq_2_id , component_2_start , component_2_stop ] )
						match_2_list.append( [ component_2_category , seq_2_id + "|+" , component_2_start , component_2_stop ] )
						match_2_list.append( [ component_2_category , seq_2_id + "|-" , component_2_start , component_2_stop ] )
						match_2_list.append( [ component_2_category , seq_2_id + "|." , component_2_start , component_2_stop ] )
					else :
						print >> sys.stdout , "[WARNING] Sequence " + seq_2_id + " assigned to group " + group_id + " is unknown, ignored"
				
				matching_pairs = list(itertools.product( match_1_list , match_2_list ))
				#print >> matching_regions_file_connection , "match_1_list:"
				#print >> matching_regions_file_connection , match_1_list
				#print >> matching_regions_file_connection , "match_2_list:"
				#print >> matching_regions_file_connection , match_2_list
				#print >> matching_regions_file_connection , "matching_pairs:"
				#print >> matching_regions_file_connection , matching_pairs

				for regions_pair in matching_pairs :
					#print >> matching_regions_file_connection , "regions_pair:"
					#print >> matching_regions_file_connection , regions_pair
					new_element = regions_pair[0][:]
					#print >> matching_regions_file_connection , new_element
					new_element += regions_pair[1]
					#print >> matching_regions_file_connection , new_element
					new_element.append(group_id)
					#print >> matching_regions_file_connection , new_element
					print >> matching_regions_file_connection , "\t".join([str(x) for x in new_element])

	matching_regions_file_connection.close()
	# matching_regions = [ ... , ["hap1_Input" , Tid , Tstart , Tstop , "Query_Input" , Qid , Qstart , Qstop , group_id] , ... ]
	# matching_regions = [ ... , ["hap1_Legacy" , Tid , Tstart , Tstop , "Query_Legacy" , Qid , Qstart , Qstop , group_id] , ... ]
	return matching_regions_file , groups_by_sequence


def make_seq_pair_from_constrains(matching_regions_file, known_input_groups, unwanted_input_pairs, alternative_pairs , agp_structure_db, hap1_list, hap2_list, input_list, component_len_db, agp_origin) :
	matching_regions_file_connection = open(matching_regions_file , 'w')
	# TODO: add alternative_pairs to the info for matching regions

	components_positions = {}
	# agp_structure_db[seq_id][start] = Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation
	# with pseudomolecules_agp_db : seq_id in chr_id and input_seq_id and CompntId in input_seq_id >> component_len_db has input_seq_id lengths so input_seq_id_to_input_seq_id relationship are produced for use on unplaced
	for chr_id in agp_structure_db.keys() :
		for start in agp_structure_db[chr_id].keys() :
			Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = agp_structure_db[chr_id][start]
			if CompntId not in components_positions :
				components_positions[CompntId] = []
				if CompntId in component_len_db :
					CompntLen = component_len_db[CompntId]
					components_positions[CompntId].append( [ CompntId , "1" , str(CompntLen) ] )
			components_positions[CompntId].append([Obj_Name , Obj_start , Obj_End , CompntId] )
			# components_positions[CompntId_1] = [ ... , [Obj_Name , Obj_start , Obj_End , CompntId] , ... ]

	# known_input_groups[group_a] = [ ... , CompntId_1 , ... ]
	# seq_id has no orientation
	# TODO: Generate relationships from group info
	for group_id in known_input_groups.keys() :
		grouped_pairs = list(itertools.product(known_input_groups[group_id], known_input_groups[group_id]))
		for pair in grouped_pairs :
			component_1 , component_2 = pair
			if component_1 == component_2 :
				continue
			else :
				component_1_regions = components_positions[component_1]
				match_1_list = []
				# component_1_regions = [ ... , [ seq_1_id , component_1_start , component_1_stop , component_1] , ... ]
				component_2_regions = components_positions[component_2]
				match_2_list = []
				for region_1 in component_1_regions :
					seq_1_id , component_1_start , component_1_stop , component_1_id = region_1
					if seq_1_id in hap1_list :
						component_1_category = "hap1_" + agp_origin
						match_1_list.append( [ component_1_category , seq_1_id , component_1_start , component_1_stop ] )
					elif seq_1_id in hap2_list :
						component_1_category = "hap2_" + agp_origin
						match_1_list.append( [ component_1_category , seq_1_id , component_1_start , component_1_stop ] )
					elif seq_1_id in input_list :
						component_1_category = "input_" + agp_origin
						match_1_list.append( [ component_1_category , seq_1_id , component_1_start , component_1_stop ] )
						match_1_list.append( [ component_1_category , seq_1_id + "|+" , component_1_start , component_1_stop ] )
						match_1_list.append( [ component_1_category , seq_1_id + "|-" , component_1_start , component_1_stop ] )
						match_1_list.append( [ component_1_category , seq_1_id + "|." , component_1_start , component_1_stop ] )
					else :
						print >> sys.stdout , "[WARNING] Sequence " + seq_1_id + " assigned to group " + group_id + " is unknown, ignored"
				for region_2 in component_2_regions :
					seq_2_id , component_2_start , component_2_stop , component_2_id = region_2
					if seq_2_id in hap1_list :
						component_2_category = "hap1_" + agp_origin
						match_2_list.append( [ component_2_category , seq_2_id , component_2_start , component_2_stop ] )
					elif seq_2_id in hap2_list :
						component_2_category = "hap2_" + agp_origin
						match_2_list.append( [ component_2_category , seq_2_id , component_2_start , component_2_stop ] )
					elif seq_2_id in input_list :
						component_2_category = "input_" + agp_origin
						match_2_list.append( [ component_2_category , seq_2_id , component_2_start , component_2_stop ] )
						match_2_list.append( [ component_2_category , seq_2_id + "|+" , component_2_start , component_2_stop ] )
						match_2_list.append( [ component_2_category , seq_2_id + "|-" , component_2_start , component_2_stop ] )
						match_2_list.append( [ component_2_category , seq_2_id + "|." , component_2_start , component_2_stop ] )
					else :
						print >> sys.stdout , "[WARNING] Sequence " + seq_2_id + " assigned to group " + group_id + " is unknown, ignored"

				matching_pairs = list(itertools.product( match_1_list , match_2_list ))
				for regions_pair in matching_pairs :
					new_element = regions_pair[0][:]
					new_element += regions_pair[1]
					new_element.append(group_id)
					print >> matching_regions_file_connection , "\t".join([str(x) for x in new_element])

	# For each pair of ids in unwanted_input_pairs generate a new matching region
	# unwanted_input_pairs[seq_1_id|+] = [ seq_2_id|+ , seq_2_id|- , seq_2_id|. , seq_3_id|+ , ... ]
	# sequences are reported in all orientations
	group_id = "exclusion"
	for component_1 in unwanted_input_pairs.keys() :
		matching_components = unwanted_input_pairs[component_1]
		component_1_regions = components_positions[component_1]

		for region_1 in component_1_regions :
			seq_1_id , component_1_start , component_1_stop , component_1_id = region_1
			if seq_1_id in hap1_list :
				component_1_category = "hap1_" + agp_origin
			elif seq_1_id in hap2_list :
				component_1_category = "hap2_" + agp_origin
			elif seq_1_id in input_list :
				component_1_category = "input_" + agp_origin
			else :
				print >> sys.stdout , "[WARNING] Sequence " + seq_1_id + " assigned to group " + group_id + " is unknown, ignored"
				continue
			line_part_1 = [ component_1_category , seq_1_id , component_1_start , component_1_stop ]

			for component_2 in matching_components :
				component_2_regions = components_positions[component_2]
				seq_2_id , component_2_start , component_2_stop , component_2_id = region_1
				if seq_2_id in hap1_list :
					component_2_category = "hap1_" + agp_origin
				elif seq_2_id in hap2_list :
					component_2_category = "hap2_" + agp_origin
				elif seq_2_id in input_list :
					component_2_category = "input_" + agp_origin
				else :
					print >> sys.stdout , "[WARNING] Sequence " + seq_2_id + " assigned to group " + group_id + " is unknown, ignored"
					continue
				line_part_2 = [ component_2_category , seq_2_id , component_2_start , component_2_stop ]

				new_line = line_part_1[:]
				new_line += line_part_2
				new_element.append(group_id)
				print >> matching_regions_file_connection , "\t".join([str(x) for x in new_element])

	matching_regions_file_connection.close()
	return matching_regions_file


def is_list_overlapping_list( list1 , list2) :
	overlapping = False
	for element in list1 :
		if element in list2 :
			overlapping = True
	return overlapping


def report_marker_usage( markers_bed_file , marker_map_by_seq , marker_map_by_id , agp_db , legacy_agp , hap1_to_chr , hap2_to_chr , unpl_to_chr , temp_folder) :
	# agp_db[chr][start]=agp_line.split("\t")
	hits_on_seq = {}

	bed_regions = read_bed_sorted_list(markers_bed_file)
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Translating marker coordinates"
	seq_to_chr = {}
	marker_hits_by_seq = {}
	#marker_hits_by_id = {}

	for seq_id in hap1_to_chr.keys() :
		seq_to_chr[seq_id] = hap1_to_chr[seq_id]
	for seq_id in hap2_to_chr.keys() :
		seq_to_chr[seq_id] = hap2_to_chr[seq_id]
	for seq_id in unpl_to_chr.keys() :
		seq_to_chr[seq_id] = unpl_to_chr[seq_id]

	for line in open( markers_bed_file ) :
		if line == "" :
			continue
		if line[0] == "#" :
			continue
		seq_id , pos_1 , pos_2 , marker_id = line.rstrip().split("\t")
		if marker_id in marker_map_by_id :
			marker_pos = marker_map_by_id[marker_id][1]
			marker_chr = marker_map_by_id[marker_id][0]
			# Keep markers only in the expected chromosome
			if seq_id in seq_to_chr :
				if marker_chr == seq_to_chr[seq_id] :
					if seq_id not in marker_hits_by_seq :
						marker_hits_by_seq[seq_id] = []
					start = min( int(pos_1) , int(pos_2) )
					stop = max( int(pos_1) , int(pos_2) )
					marker_hits_by_seq[seq_id].append([ marker_chr , marker_pos , marker_id , seq_id , int(start) , int(stop) ] )
			#else :
			#	if seq_id not in marker_hits_by_seq :
			#		marker_hits_by_seq[seq_id] = []
			#	start = min( int(pos_1) , int(pos_2) )
			#	stop = max( int(pos_1) , int(pos_2) )
			#	marker_hits_by_seq[seq_id].append([ marker_chr , marker_pos , marker_id , seq_id , int(start) , int(stop) ] )


	for chr in agp_db.keys() :
		for start in agp_db[chr].keys() :
			component_ID = agp_db[chr][start][5]
			component_orientation = agp_db[chr][start][8]
			seq_to_chr[component_ID] = [ chr , component_orientation ]

	for chr in legacy_agp.keys() :
		for start in legacy_agp[chr].keys() :
			component_ID = legacy_agp[chr][start][5]
			component_orientation = legacy_agp[chr][start][8]
			seq_to_chr[component_ID] = [ chr , component_orientation ]

	out_bed_file_name = temp_folder + "/input.markers.bed"
	out_bed_file = open(out_bed_file_name , 'w')
	for line in open(markers_bed_file , 'r') :
		print >> out_bed_file, line.rstrip()
	out_bed_file.close()

	# Read marker position
	for line in open( out_bed_file_name ) :
		if line == "" :
			continue
		if line[0] == "#" :
			continue
		seq_id , pos_1 , pos_2 , marker_id = line.rstrip().split("\t")
		if marker_id in marker_map_by_id :
			marker_pos = marker_map_by_id[marker_id][1]
			marker_chr = marker_map_by_id[marker_id][0]
			# Keep markers only in the expected chromosome
			if seq_id in seq_to_chr :
				if marker_chr == seq_to_chr[seq_id] :
					if seq_id not in marker_hits_by_seq :
						marker_hits_by_seq[seq_id] = []
					start = min( int(pos_1) , int(pos_2) )
					stop = max( int(pos_1) , int(pos_2) )
					marker_hits_by_seq[seq_id].append([ marker_chr , marker_pos , marker_id , seq_id , int(start) , int(stop) ] )


	legacy_seq_bed = translate_bed_sorted_list( bed_regions , invert_agp(legacy_agp) )
	legacy_seq_bed_file_name = temp_folder + "/legacy.markers.bed"
	out_bed_file = open(legacy_seq_bed_file_name , 'w')
	for line in sorted(legacy_seq_bed) :
		print >> out_bed_file, "\t".join([str(x) for x in line])
	out_bed_file.close()

	# Read marker position
	for line in open( legacy_seq_bed_file_name ) :
		if line == "" :
			continue
		if line[0] == "#" :
			continue
		seq_id , pos_1 , pos_2 , marker_id = line.rstrip().split("\t")
		if marker_id in marker_map_by_id :
			marker_pos = marker_map_by_id[marker_id][1]
			marker_chr = marker_map_by_id[marker_id][0]
			# Keep markers only in the expected chromosome
			if seq_id not in marker_hits_by_seq :
				marker_hits_by_seq[seq_id] = []
			start = min( int(pos_1) , int(pos_2) )
			stop = max( int(pos_1) , int(pos_2) )
			marker_hits_by_seq[seq_id].append([ marker_chr , marker_pos , marker_id , seq_id , int(start) , int(stop) ] )

	for seq_id in marker_hits_by_seq.keys() :
		hits_on_seq[seq_id] = {}
		hits_on_seq[seq_id]["id"] = seq_id
		hits_on_seq[seq_id]["chr"] = seq_to_chr[seq_id][0]
		hits_on_seq[seq_id]["orientation"] = seq_to_chr[seq_id][1]
		hits_on_seq[seq_id]["markers"] = marker_hits_by_seq[seq_id]
		marker_pos_min = min( [ int(x[1]) for x in marker_hits_by_seq[seq_id] ] )
		marker_pos_max = max( [ int(x[1]) for x in marker_hits_by_seq[seq_id] ] )
		hits_on_seq[seq_id]["range"] = [marker_pos_min , marker_pos_max]

	# Output:
		# 	hits_on_seq[seq_id]["id"] : seq_id
		# 	hits_on_seq[seq_id]["chr"] : chr_id
		# 	hits_on_seq[seq_id]["markers"] : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ]
		# 	hits_on_seq[seq_id]["range"] : [marker_pos_min , marker_pos_max] ,
		# 	hits_on_seq[seq_id]["orientation"] : ["+" or "-" or "."] }

	return hits_on_seq

