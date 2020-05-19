#!/usr/bin/env python

import sys
import subprocess
import argparse
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import gc
import os
import pipes


gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


###### Functions ######


def get_subgraph_from_path( original_graph , selected_path ):

	edges = []
	edges_names=[]

	for i in range( 0 , len(selected_path) -1 ):
		from_node = selected_path[i]
		to_node = selected_path[i+1]
		element = original_graph[from_node][to_node]
		edges_names.append(element["name"])
		edges.append([ element["name"] , "(" + str(from_node) + ":" + str(to_node) + ")" ,  element["region"] , element["align"] , element["match"] ])

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

		print >> sys.stderr , "#### " + Qid + " blacklisted: " + str(Qid in blacklist)

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
				hit_graph.add_edge( T_from , T_to , align=0 , match=0, name=T_gap, region="gap")
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
		hit_graph.add_edge( (Tstart,Qstart) , (Tstop,Qstop) , align=int(align) , match=int(matches))

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
			if T_gap >= 0  and  T_gap <= max_distance and Q_gap >= 0 and Q_gap <= max_distance :
				if not hit_graph.has_edge( (f_Tstop , f_Qstop) , (t_Tstart , t_Qstart) ) :
					if (f_Tstop , f_Qstop) != (t_Tstart , t_Qstart) :
						hit_graph.add_edge( (f_Tstop , f_Qstop) , (t_Tstart , t_Qstart) , align=0 , match=0)
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
				print >> sys.stderr , "###### Adding " + Qid + " to graph"
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

				edge = hit_graph[int(Tstart)][int(Tstop)]
				print >> sys.stderr , "####### edge: " + str(Tstart) + " " + str(Tstop) + " - " + ", ".join( str(key) + ": " + str(edge[key]) for key in sorted(edge.keys()) )

			else:
				print >> sys.stdout , "----- Error: " + Qid + " present bot as required and blacklisted."
				print >> sys.stderr , "##### Error: " + Qid + " present bot as required and blacklisted."
				sys.exit(3)

	# Forced Qids mappings sanity check
	print >> sys.stdout , "---- Testing required query sequences mapping positions order"
	prev_id , prev_start , prev_stop = [ "" , "" ,"" ]

	for id in forced_sorted_list :

		if id not in forced_nodes :
			print >> sys.stdout , "----- " + id + " not mapped, removing from required"
			print >> sys.stderr , "##### " + id + " not mapped, removing from required"
			continue

		start = int(forced_nodes[id][0])
		stop = int(forced_nodes[id][1])

		if not prev_id == "" :
			# Not the first element, test regions
			if (stop <= prev_stop) or (start <= prev_start) :
				# Mapping order not compatible
				print >> sys.stdout , "----- Incompatibility issue: (" + prev_id + " -> " + id + ") order not coherent with mapping results. Check " +  err_filename + " for further information"
				print >> sys.stderr , "##### ERROR: (" + prev_id + " -> " + id + ") order not coherent with mapping results. Check " +  err_filename + " for further information"
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
					print >> sys.stdout , "----- (" + prev_id + " -> " + id + ") overlap in mapping results."
					new_start = prev_stop
					hit_graph.add_edge( new_start , stop)
					hit_graph[new_start][stop].update(hit_graph[start][stop])
					hit_graph.remove_node(start)

		# Update
		prev_id = id
		prev_start = start
		prev_stop = stop

	print >> sys.stderr , "##### Forced query sequences mappings"
	print >> sys.stderr , "##### ID\tAlignment_Start\tAlignment_Stop"
	for id in forced_sorted_list :
		if id in forced_nodes :
			print >> sys.stderr ,"##### " + id + "\t" + str(forced_nodes[id][0]) + "\t" + str(forced_nodes[id][1])


	# Make nodes and add mapping edges of Qid compatible with forced list

	for hit in hit_list :
		Qid , Tstart , Tstop , Qstart , Qstop , matches , hitLen = hit

		print >> sys.stderr , "#### " + Qid + " Required: " +  str(Qid in forced_sorted_list) + " ; blacklisted: " + str(Qid in blacklist)
		start_nodes.sort()
		stop_nodes.sort()

		if  Qid not in forced_sorted_list:
			if Qid not in blacklist :
				# Compatibility test
				add = True
				for query_range in forced_nodes_ranges :
					if ( query_range[0] <= Tstart <= query_range[1] ) or (query_range[0] <= Tstop <= query_range[1]) : add = False

				if add :
					print >> sys.stderr , "###### Adding " + Qid + " to graph"
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
						print >> sys.stderr , "####### edge: " + str(Tstart) + " " + str(Tstop) + " - " + ", ".join( str(key) + ": " + str(edge[key]) for key in sorted(edge.keys()) )

					else :
						print >> sys.stderr , "###### Chromosome extremity"
						if int(Tstart) == 0 :
							print >> sys.stderr , "####### Chromosome start"
							print >> sys.stderr, sorted(stop_nodes)
							stop_nodes.append(int(Tstop))
							stop_nodes.sort()
							print >> sys.stderr, sorted(stop_nodes)
						else :
							print >> sys.stderr , "####### Chromosome end"
							print >> sys.stderr, sorted(start_nodes, reverse=True)
							start_nodes.append(int(Tstart))
							start_nodes.sort()
							print >> sys.stderr, sorted(start_nodes, reverse=True)
				else:
					print >> sys.stderr , "###### Refused " + Qid + " to graph, incompatible"

	### Add allowed gap edges between nodes

	stop_nodes_len = int(len(stop_nodes))
	start_nodes_len = int(len(start_nodes))

	print >> sys.stderr, "##### Graph mapping contigs edges"
	print >> sys.stderr, sorted(hit_graph.edges())
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
					hit_graph.add_edge( T_from , T_to , align=0 , match=0, name=T_gap, region="gap")
					print >> sys.stderr , "(" + str(T_from) + ") -> (" + str(T_to) + ") | Gap: "  + str(T_gap)

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


def make_agp_from_list( querylist , queryfastalen , gaplen , seqoutname , outfilename) :

	############################################################################################################################################################
	# AGP file format, tab delimited, 1-base coordinates #######################################################################################################
	############################################################################################################################################################
	#	##agp-version	2
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

	for CompntId in querylist :

		Id, T_range , Q_range , alignment , matches = CompntId

		if not Q_range == "gap" :
			# Skip alignment gaps in the list

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


def make_fasta_from_list( querylist , queryfasta , gaplen , seqoutname , outfilename) :

	### Query list element format
	# Gap:	[61252	,	(0:61252)		,	gap		,	61252	,	0]
	#		[length	,	(T_start:Tstop)	, 	"gap" 	, 	length 	, 	0]
	# Object:	[b40-14.HS_iter4_seq3733|+	,	(61252:6463804)	,	93612:7595148	,	-6402552			,	4526208]
	#			[ID|strand					,	(T_start:Tstop)	,	Q_start:Q_stop	,	-(alignment length)	,	matches]

	gaplen = int(gaplen)

	seq = Seq("")
	seq.id = seqoutname

	for CompntId in querylist :

		Id, T_range , Q_range , alignment , matches = CompntId

		if not Q_range == "gap" :

			CompntId_name = Id[:-2]
			Orientation = Id[-1]

			# Add gap between Components
			if str(seq) != "" :
				seq = seq + "N" * gaplen

			if Orientation == "-" :
				my_sub_seq = queryfasta[CompntId_name].reverse_complement()
			else :
				my_sub_seq = queryfasta[CompntId_name]

			seq = seq + my_sub_seq

	# Print the entire sequence
	seq.id = seqoutname
	seq.description = ""
	print >> outfilename, seq.format('fasta')


def hit_mu( paf_file , max_gap_size , rlen , qlen ) :
	print >> sys.stdout, "- Merge hits in blocks and select best alignments"
	original_hits = {}
	merged_hits = {}
	unique_hits = {}

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
		except :
			print >> sys.stderr, "#### Error merging hits of " + Qid + " on " + Tid + " (" + str(len(original_hits[entry])) + " hits)"
			print >> sys.stderr, sorted(map_graph.edges.data())
			print >> sys.stderr, nx.find_cycle(map_graph)
			exit(5)
		longest_graph = get_subgraph_from_path_tuples( map_graph , longest_merged_path )
		longest_graph_matches = longest_graph.size(weight='match')

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

	unique_temp_file = open( paf_file + ".uniq" , "w")

	print >> sys.stdout, "-- Selecting best alignment"

	for Qid in sorted(merged_hits.keys()) : # For each absolute id find the best alignment
		merged_list = sorted(merged_hits[Qid] , key=lambda item: int(item[10]) , reverse = True )
		#print >> sys.stderr, merged_list
		unique_hits[Qid] = merged_list[0]
		print >> unique_temp_file , "\t".join(str(x) for x in merged_list[0])

	return unique_hits


