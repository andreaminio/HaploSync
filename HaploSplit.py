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
from HaploFunct import *

gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

###### Options and help ######

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="hits",
				help="File of hits. Used as input if already mapped, output result of mapping procedure if -m/--map", metavar="file.paf [Required]")
parser.add_argument("-r", "--reference", dest="reference",
				help="FASTA file with reference sequence", metavar="reference.fasta [Required]")
parser.add_argument("-q", "--query", dest="query",
				help="FASTA file with query sequence", metavar="query.fasta [Required]")
parser.add_argument("-v", "--dry", dest="dry", action="store_true",
				help="Dry run: find best tiling paths and excluded nodes but no FASTA or AGP exported")
parser.add_argument("-o", "--out", dest="out", default="out",
				help="Output NAME prefix for files [default: out]", metavar="NAME")
parser.add_argument("-p", "--prefix", dest="prefix", default="NEW",
				help="Output NAME prefix for files [default: NEW]", metavar="PREFIX")
parser.add_argument("-a", "--agp", dest="agp", action="store_true",
				help="Export paths in AGP format to NAME.{1,2,Un}.agp")
parser.add_argument("-f", "--fasta", dest="sequence", action="store_true",
				help="Export paths sequences in NAME.{1,2,Un}.fasta")
parser.add_argument("-g", "--gapsize", dest="gap", default="1000",
				help="Minimum gap size placeholder for AGP and FASTA (in bp) [default: 1,000bp]" , metavar="N")
parser.add_argument("-c", "--concatenate", dest="conc", default="",
				help="Set the gap size used to concatenate unplaced sequences (in bp) [default: do not concatenated]", metavar="N")
parser.add_argument("--distance1", dest="distance1", default="2000000",
				help="Set the maximum distance (bp) between two hits to be considered adjacent [default: 2,000,000bp]", metavar="N")
parser.add_argument("--distance2", dest="distance2", default="4000000",
				help="Set the maximum distance (bp) between two hits to be considered adjacent at second round [default: 4,000,000bp]",metavar="N" )
parser.add_argument("-m", "--map" , dest="map", action="store_true",
				help="Do run mapping [overwrite input file.paf if existing]")
parser.add_argument("-u", "--uniq" , dest="unique", action="store_true",
				help="Map only" )
parser.add_argument("-t", "--tool", dest="mapping", default="minimap2 --cs -t 4 -x asm20 -r 1000",
				help="Mapping command for minimap2", metavar="\"minimap2 --cs -t 4 -x asm20 -r 1000\"")
parser.add_argument("-n" , "--hitgap" , dest="hitgap", default="500000",
				help="Allowed gap between hits to be merged [default: 500000]", metavar="N")
parser.add_argument("-1", "--path1", dest="use1",
				help="Use ONLY the following list of (stranded) contigs to create the first haplotype. Tab separated file with Target_ID followed by comma separated list of contigs (contig1,contig2,[...],contigN)", metavar="1st.txt")
parser.add_argument("-2", "--path2", dest="use2",
				help="Use ONLY the following list of (stranded) contigs to create the second haplotype. Tab separated file with Target_ID followed by comma separated list of contigs (contig1,contig2,[...],contigN)", metavar="2nd.txt")
parser.add_argument("--R1", dest="Require1",
				help="Require to use the following list of (stranded) contigs in the first haplotype. Tab separated file with Target_ID followed by comma separated list of contigs (contig1,contig2,[...],contigN)", metavar="1st.txt")
parser.add_argument("--min1", dest="minR1", default="0" ,
				help="Minimum length of sequence allowed to be a requirement for the first haplotype (default: 0)", metavar="N")
parser.add_argument("--R2", dest="Require2",
				help="Require to use the following list of (stranded) contigs in the second haplotype. Tab separated file with Target_ID followed by comma separated list of contigs (contig1,contig2,[...],contigN)", metavar="2nd.txt")
parser.add_argument("--min2", dest="minR2", default="0" ,
				help="Minimum length of sequence allowed to be a requirement for the second haplotype (default: 0)", metavar="N")
parser.add_argument("--B1", dest="Blacklist1",
				help="Blacklisted (stranded) contigs NOT to be used in the first haplotype. Tab separated file with Target_ID followed by comma separated list of contigs (contig1,contig2,[...],contigN)", metavar="2nd.txt")
parser.add_argument("--B2", dest="Blacklist2",
				help="Blacklisted (stranded) contigs NOT to be used in the second haplotype. Tab separated file with Target_ID followed by comma separated list of contigs (contig1,contig2,[...],contigN)", metavar="2nd.txt")
parser.add_argument("--N2", dest="No2", action="store_true",
				help="Don't run the search for the 2nd path")




# Sanity Check

if len(sys.argv) < 2:
	parser.print_help()
	sys.exit(1)

options = parser.parse_args()

if not options.hits :
	print >> sys.stderr , "[ERROR] Hit file FASTA file missing"
	parser.print_help()
	sys.exit(1)

if not options.reference :
	print >> sys.stderr , "[ERROR] Reference genome FASTA file missing"
	parser.print_help()
	sys.exit(1)

if not options.query :
	print >> sys.stderr , "[ERROR] Query genome FASTA file missing"
	parser.print_help()
	sys.exit(1)

#best_1_index = int( options.use1 ) - 1
#best_2_index = int( options.use2 ) - 1

###### Main ######
### Print intro
print >> sys.stdout, "Running HaploSplit - HaploSync v0.1beta 
print >> sys.stdout, "- Running tiling path search"
print >> sys.stdout, "-- To reproduce this run use the following command: " + " ".join( pipes.quote(x) for x in sys.argv)

### Read reference and compute lengths
print >> sys.stdout, "- Reading reference sequences"
reference = dict((seq.id, seq) for seq in SeqIO.parse(open(options.reference), "fasta"))
reference_len = {}
print >> sys.stdout, "-- Calculating reference sequences lengths"
for ref_name in reference :
	reference_len[ref_name] = int(len(reference[ref_name]))
print >> sys.stdout, "--- Number of reference sequences: " + str(len(reference_len.keys()))

### Read query and compute lengths
print >> sys.stdout, "- Reading query sequences"
query = dict((seq.id, seq) for seq in SeqIO.parse(open(options.query), "fasta"))
query_len = {}
print >> sys.stdout, "-- Calculating query sequences lengths"
for query_name in query :
	query_len[query_name] = int(len(query[query_name]))
print >> sys.stdout, "--- Number of query sequences: " + str(len(query_len.keys()))

### For each chromosome, Create graph
best_1_paths = {}
best_1_paths_edges = {}
best_2_paths = {}
best_2_paths_edges ={}
all_used = []


# Force path if necessary

fixed_path_1 = defaultdict(list)
fixed_path_2 = defaultdict(list)

if options.use1 or options.use2 :
	print >> sys.stdout , "- Reading given paths"

if options.use1 :
	print >> sys.stdout , "-- Sequences used in 1st path"
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
			print >> sys.stdout , "--- Seq " + str(Target_sequence) + " : " + ",".join(best_1_paths_edges[Target_sequence])
		except:
			pass

if options.use2 :
	print >> sys.stdout , "-- Sequences used in 2nd path"
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
			print >> sys.stdout , "--- Seq " + str(Target_sequence) + " : " + ",".join(best_2_paths_edges[Target_sequence])
		except :
			pass


# Create required lists
forced_list_1 = defaultdict(list)
forced_list_2 = defaultdict(list)
discarded=[]

if options.Require1 or options.Require2 :
	print >> sys.stdout , "- Reading required sequence lists"

if options.Require1 :
	print >> sys.stdout , "-- Sequence required in 1st path"
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
				else :
					discarded.append(id_name)
			print >> sys.stdout , "--- Seq " + str(Target_sequence) + ": " + ",".join(forced_list_1[Target_sequence])
		except :
			pass

if options.Require2 :
	print >> sys.stdout , "-- Sequence required in 2nd path"
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
				else :
					discarded.append(id_name)
			print >> sys.stdout , "--- Seq " + str(Target_sequence) + ": " + ",".join(forced_list_2[Target_sequence])
		except :
			pass

if not discarded == [] :
	print >> sys.stdout , "-- Discarded required sequences because of short length (hap1 <" + str(options.minR1) +" or hap2 <" + str(options.minR2) + ")"
	print >> sys.stdout , "--- " + ",".join(discarded)

# Create Blacklists
blacklist_1 = defaultdict(list)
blacklist_2 = defaultdict(list)

for ref_name in reference :
	blacklist_1[ref_name]=[]
	blacklist_2[ref_name]=[]

print >> sys.stdout , "- Blacklists"
if options.Blacklist1 :
	print >> sys.stdout , "-- Rejected for 1st path"
	for line in open(options.Blacklist1,"r") :
		try :
			Tid , Qids = line.rstrip().split("\t")
			for name in Qids.split(",") :
				blacklist_1[Tid].append(name[:-1]+"+")
				blacklist_1[Tid].append(name[:-1]+"-")
			print >> sys.stdout , "--- Seq " + str(Tid) + " : " + ",".join(blacklist_1[Tid])
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
		blacklist_1[key_in].append(name[:-1]+"+")
		blacklist_1[key_in].append(name[:-1]+"-")

	print >> sys.stderr , "### To reject for Path1@" + str(key_in) + " : " + ",".join(blacklist_1[key_in])


if options.Blacklist2 :
	print >> sys.stdout , "-- Rejected for 2nd path"
	for line in open(options.Blacklist2,"r") :
		try :
			Tid , Qids = line.rstrip().split("\t")
			for name in Qids.split(",") :
				blacklist_2[Tid].append(name[:-1]+"+")
				blacklist_2[Tid].append(name[:-1]+"-")
			print >> sys.stdout , "--- Seq " + str(Tid) + " : " + ",".join(blacklist_2[Tid])
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
		blacklist_2[key_in].append(name[:-1]+"+")
		blacklist_2[key_in].append(name[:-1]+"-")

	print >> sys.stderr , "### To reject for Path2@" + str(key_in) + " : " + ",".join(blacklist_2[key_in])



if options.map :
	doubleQuery = {}
	doubleQuery_len = {}
	print >> sys.stdout, "- Mapping query sequences on reference"

	query_for = "tmp." + options.query + ".for"

	qf_file = open(query_for , 'w' )

	for query_seq in query :
		print >> qf_file , ">" + query_seq + "|+"
		print >> qf_file , str(query[query_seq].seq).upper()
		doubleQuery[query_seq+ "|+"]=str(query[query_seq].seq).upper()
		doubleQuery_len[query_seq+ "|+"]=int(len(str(query[query_seq].seq)))
	qf_file.close()

	print >> sys.stdout, "-- Mapping forward query sequences"
	map_for_file = open("tmp." + options.hits + ".for" , "w")
	mappingCommand = options.mapping + " --for-only " + options.reference + " " + query_for + " | awk \'$5==\"+\"\'"
	print >> sys.stdout, "--- Command line: " + mappingCommand + " > tmp." + options.hits + ".for"
	mapProcess = subprocess.Popen(mappingCommand, shell=True, stdout=map_for_file)
	output, error = mapProcess.communicate()

	print >> sys.stdout, "-- Reversing query sequences"
	query_rev = "tmp." + options.query + ".rev"

	qr_file = open(query_rev , 'w' )

	for query_seq in query :
		print >> qr_file , ">" + query_seq + "|-"
		print >> qr_file , str(query[query_seq].reverse_complement().seq).upper()
		doubleQuery[query_seq+ "|-"]=str(query[query_seq].reverse_complement().seq).upper()
		doubleQuery_len[query_seq+ "|-"]=int(len(str(query[query_seq].seq)))
	qr_file.close()

	print >> sys.stdout, "-- Mapping reverse query sequences"
	map_rev_file = open("tmp." + options.hits + ".rev" , "w")
	mappingCommand = options.mapping + " --for-only " + options.reference + " " + query_rev + " | awk \'$5==\"+\"\'"
	print >> sys.stdout, "--- Command line: " + mappingCommand + " > tmp." + options.hits + ".rev"
	mapProcess = subprocess.Popen(mappingCommand, shell=True, stdout=map_rev_file)
	output, error = mapProcess.communicate()

	map_multi_file = open("tmp." + options.hits + ".multimapping" , "w")
	print >> sys.stdout, "-- Merge forward and reverse alignments"
	mapProcess = subprocess.Popen("cat tmp." + options.hits + ".for tmp." + options.hits + ".rev", shell=True, stdout=map_multi_file)
	output, error = mapProcess.communicate()

	map_multi = "tmp." + options.hits + ".multimapping"
	mapProcess = subprocess.Popen("cp " + map_multi + " " + options.hits , shell=True)
	output, error = mapProcess.communicate()

else :
	map_multi = "tmp." + options.hits + ".multimapping"
	mapProcess = subprocess.Popen("cp " + options.hits + " " + map_multi , shell=True)
	output, error = mapProcess.communicate()

	doubleQuery = {}
	doubleQuery_len = {}
	for query_seq in query :
		doubleQuery[query_seq+ "|+"]=str(query[query_seq].seq).upper()
		doubleQuery_len[query_seq+ "|+"]=int(len(str(query[query_seq].seq)))
		doubleQuery[query_seq+ "|-"]=str(query[query_seq].reverse_complement().seq).upper()
		doubleQuery_len[query_seq+ "|-"]=int(len(str(query[query_seq].seq)))

#### Map done, check if
if options.unique :
	print >> sys.stdout , "Mapping done, quitting"
	sys.exit(0)

#### Merge and unique ####

unique_hits = hit_mu("tmp." + options.hits+ ".multimapping" , int(options.hitgap) , reference_len , doubleQuery_len)

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
print >> sys.stdout, "- Generating graphs and best paths"

for Target_sequence in sorted(hits.keys()) :
	hits_1 = hits[Target_sequence]
	hits_1.append(["ChrStart" , 0 , 0 , 0 , 0 , 0 , 0 ])
	hits_1.append(["ChrStop" , reference_len[Target_sequence] , reference_len[Target_sequence] , 0 , 0 , 0 , 0 ])

	if not options.use1 :
		print >> sys.stdout, "-- Graphing: " + Target_sequence
		print >> sys.stderr, Target_sequence
		map_total= len(hits_1)

		print >> sys.stdout, "--- Mapping query sequences: " + str(map_total)
		#print >> sys.stderr, len(hits_1)
		if options.Require1 :
			hit_graph = make_forced_graph(hits_1, int(options.distance1) , forced_list_1[Target_sequence] , blacklist_1[Target_sequence]  , options.Require1 + ".err.txt" )
		else :
			hit_graph = make_graph(hits_1, int(options.distance1) , blacklist_1[Target_sequence] )
		print >> sys.stderr, hit_graph.edges.data()

		distance1 = int(options.distance1)

		while not nx.has_path(hit_graph , source=0 , target=reference_len[Target_sequence] ) :
				distance1 += 100000
				print >> sys.stdout, "---- Increased allowed gap to " + str(distance1)
				print >> sys.stderr, "#### Rebuilding increased allowed gap to " + str(distance1)
				if options.Require2 :
					hit_graph = make_forced_graph(hits_1, distance1 , forced_list_1[Target_sequence] , blacklist_1[Target_sequence]  , options.Require1 + ".err.txt" )
				else :
					hit_graph = make_graph(hits_1, distance1 , blacklist_1[Target_sequence] )


		# Select first
		print >> sys.stdout, "--- Searching 1st path maximizing coverage: " + Target_sequence + " (0:" + str(reference_len[Target_sequence]) + ")"
		print >> sys.stdout, "--- Graph order: " + str(hit_graph.order()) + "; Graph size: " + str(hit_graph.size())
		#best_1_nodes = nx.bellman_ford_path( hit_graph , source=0 , target=reference_len[Target_sequence] , weight="align" )
		best_1_nodes = nx.dag_longest_path(hit_graph, weight='align')
		print >> sys.stderr , best_1_nodes
		best_1_edges , best_1_edges_names = get_subgraph_from_path( hit_graph , best_1_nodes )
		print >> sys.stderr , best_1_edges
		#print >> sys.stderr , best_1_edges_names

		#print >> sys.stderr , best_1_edges
		#print >> sys.stderr ,"### BEST @ 1st iter: " + " -> ".join( str(x) for x in best_1_edges_names)
		print >> sys.stderr ,"### BEST mapping @ 1st iter: [" + "] -> [".join( ",".join(str(r) for r in x) for x in best_1_edges) + "]"

		best_1_paths[Target_sequence] = best_1_edges

	used , best_1_paths_edges[Target_sequence] = make_list_from_path( best_1_paths[Target_sequence] )
	all_used += used

	print >> sys.stdout , "--- Used for " + Target_sequence + ": " + ",".join(best_1_paths_edges[Target_sequence])
	print >> sys.stdout , "--- Used so far: " + ",".join(all_used)


### Remove nodes form best path and search for second best tiling
for Target_sequence in sorted(hits.keys()) :
	hits_2 = [x for x in hits[Target_sequence] if x[0] not in all_used ]
	unmap_2 = len(hits_2)

	if not options.use2 :
		print >> sys.stdout, "--- Updating graph removing best path for " + Target_sequence

		#print >> sys.stdout, "---- Mapping query sequences used in 1st path: " + str(map_total + 2 - unmap_1)
		#print >> sys.stderr, len(hits_2)
		if not options.No2 :
			### Make new graph
			print >> sys.stdout, "--- Recreating the graph"
			if options.Require2 :
				hit_graph_2 = make_forced_graph(hits_2, int(options.distance2) , forced_list_2[Target_sequence] , blacklist_2[Target_sequence] , options.Require2 + ".err.txt" )
			else :
				hit_graph_2 = make_graph(hits_2, int(options.distance2) , blacklist_2[Target_sequence] )

			# Check if graph is connected, otherwise increase allowed distance between nodes
			distance2 = int(options.distance2)

			while not nx.has_path(hit_graph_2 , source=0 , target=reference_len[Target_sequence] ) :
				distance2 += 1000000
				print >> sys.stdout, "---- Increased allowed gap to " + str(distance2)
				if options.Require2 :
					hit_graph_2 = make_forced_graph(hits_2, distance2 , forced_list_2[Target_sequence], blacklist_2[Target_sequence] , options.Require2 + ".err.txt" )
				else :
					hit_graph_2 = make_graph(hits_2, distance2 , blacklist_2[Target_sequence] )

			### Extract best 2nd
			print >> sys.stdout, "--- Graph order: " + str(hit_graph_2.order()) + "; Graph size: " + str(hit_graph_2.size())
			print >> sys.stdout, "--- Searching 2nd path maximizing coverage: " + Target_sequence + " (0:" + str(reference_len[Target_sequence]) + ")"
			#best_2_nodes = nx.bellman_ford_path( hit_graph_2 , source=0 , target=reference_len[Target_sequence] , weight="align" )
			best_2_nodes = nx.dag_longest_path(hit_graph_2, weight='align')
			#print >> sys.stderr , best_2_nodes
			best_2_edges , best_2_edges_names = get_subgraph_from_path( hit_graph_2 , best_2_nodes )
			#print >> sys.stderr , best_2_edges
			#print >> sys.stderr , best_2_edges_names

			#print >> sys.stderr , best_2_edges
			#print >> sys.stderr ,"### BEST @ 2nd iter: " + " -> ".join( str(x) for x in best_2_edges_names)
			print >> sys.stderr ,"### BEST mapping @ 2nd iter: [" + "] -> [".join( ",".join(str(r) for r in x) for x in best_2_edges) + "]"
			best_2_paths[Target_sequence] = best_2_edges

		else :
			best_2_edges_names = ["ChrStart" , "ChrStop"]


	used , best_2_paths_edges[Target_sequence] = make_list_from_path( best_2_paths[Target_sequence] )
	all_used += used

	print >> sys.stdout , "----- Used for " + Target_sequence + ": " + ",".join(best_2_paths_edges[Target_sequence])
	print >> sys.stdout , "----- Used so far: " + ",".join(all_used)


print >> sys.stdout, "---- Find leftover query sequences"
#### Extract list of unused sequences

all_unused = []

for id in sorted(query_len.keys()) :
	if str(id) + "|+" in all_used :
		# Query sequence, forward or reversed, used for assembling
		continue
	else :
		all_unused.append(id)

print >> sys.stdout, "--- " + str(len(all_unused)) + " mapping sequences left unused"
print >> sys.stderr ,"### Leftover sequences: "+ " , ".join(all_unused)
print >> sys.stderr ,"------------------------------"

#### Print output files

gap_length = int(options.gap)


#### Print lists

if not options.dry :

	list_1_file = open(options.out + ".1" + ".list" , "w")
	if not options.No2 :
		list_2_file = open(options.out + ".2" + ".list" , "w")
	list_Un_file = open(options.out + ".Un" + ".list" , "w")
	print >> sys.stdout , "- Printing sorted list of sequences in tiling paths"

	for Target_sequence in sorted(best_1_paths.keys()):
		used , Ids = make_list_from_path( best_1_paths[Target_sequence] )
		print >> list_1_file , Target_sequence + "\t" + ",".join(Ids)

	if not options.No2 :
		for Target_sequence in sorted( best_2_paths.keys() ):
			used , Ids = make_list_from_path( best_2_paths[Target_sequence] )
			print >> list_2_file , Target_sequence + "\t" + ",".join( Ids )

	print >> list_Un_file , ",".join(all_unused)

	list_1_file.close()
	if not options.No2 :
		list_2_file.close()
	list_Un_file.close()


##### Print path in AGP format
if options.agp and not options.dry :

	### Query list element format
	# Gap:	[61252	,	(0:61252)		,	gap		,	61252	,	0]
	#		[length	,	(T_start:Tstop)	, 	"gap" 	, 	length 	, 	0]
	# Object:	[b40-14.HS_iter4_seq3733|+	,	(61252:6463804)	,	93612:7595148	,	-6402552			,	4526208]
	#			[ID|strand					,	(T_start:Tstop)	,	Q_start:Q_stop	,	-(alignment length)	,	matches]

	agp_1_file = open(options.out + ".1" + ".agp" , "w")
	if not options.No2 :
		agp_2_file = open(options.out + ".2" + ".agp" , "w")
	agp_Un_file = open(options.out + ".Un" + ".agp" , "w")
	print >> sys.stdout , "- Printing the tiling path in AGP file (gap size " + str(gap_length) + "bp)"
	print >> sys.stdout , "-- Hap1"

	for Target_sequence in sorted(best_1_paths_edges.keys()):
		Obj_name = options.prefix + "_Hap1_" + Target_sequence
		print >> sys.stdout , "--- " + Obj_name
		make_agp_from_list( best_1_paths[Target_sequence] , query_len , gap_length , Obj_name , agp_1_file)
	agp_1_file.close()

	print >> sys.stdout , "-- Hap2"

	if not options.No2 :
		for Target_sequence in sorted(best_2_paths_edges.keys()):
			Obj_name = options.prefix + "_Hap2_" + Target_sequence
			print >> sys.stdout , "--- " + Obj_name
			make_agp_from_list( best_2_paths[Target_sequence] , query_len , gap_length , Obj_name , agp_2_file)
		agp_2_file.close()

	if not options.conc == "" :
		conc_gap = int(options.conc)
		Obj_name = options.prefix + "_Un"
		print >> sys.stdout , "-- Unplaced sequences will be concatenated (gap size" + options.conc + "bp)"
		all_unused_path = [ [x+"|+" , 0 ,0 ,0 ,0] for x in all_unused ]
		make_agp_from_list( all_unused_path , query_len , conc_gap , Obj_name , agp_Un_file)
	else :
		id = 0
		print >> sys.stdout , "-- Unplaced sequences"
		for comp in all_unused :
			id += 1
			Obj_name = options.prefix + "_Un_" + str(id)
			make_agp_from_list( [ [comp+"|+", 0 ,0 ,0 ,0] ] , query_len , 0 , Obj_name , agp_Un_file)
	agp_Un_file.close()


#### Print sequence FASTA
if options.sequence and not options.dry :

	fasta_1_file = open(options.out + ".1" + ".fasta" , "w")
	if not options.No2 :
		fasta_2_file = open(options.out + ".2" + ".fasta" , "w")
	fasta_Un_file = open(options.out + ".Un" + ".fasta" , "w")
	print >> sys.stdout , "- Printing sequences of tiling path in fasta files (gap size " + str(gap_length) + "bp)"
	print >> sys.stdout , "-- Hap1"

	for Target_sequence in sorted(best_1_paths.keys()):
		Obj_name = options.prefix + "_Hap1_" + Target_sequence
		print >> sys.stdout , "--- " + Obj_name
		make_fasta_from_list( best_1_paths[Target_sequence] , query ,  gap_length , Obj_name ,  fasta_1_file)
	fasta_1_file.close()

	print >> sys.stdout , "-- Hap2"
	if not options.No2 :
		for Target_sequence in sorted(best_2_paths.keys()):
			Obj_name = options.prefix + "_Hap2_" + Target_sequence
			print >> sys.stdout , "--- " + Obj_name
			make_fasta_from_list( best_2_paths[Target_sequence] , query ,  gap_length , Obj_name ,  fasta_2_file)
		fasta_2_file.close()

	if not options.conc == "" :
		conc_gap = int(options.conc)
		Obj_name = options.prefix + "_Un"
		print >> sys.stdout , "-- Unplaced sequences will be concatenated (gap size" + options.conc + "bp)"
		all_unused_path = [ [x+"|+" , 0 ,0 ,0 ,0] for x in all_unused ]
		make_fasta_from_list( all_unused_path , query , conc_gap , Obj_name , fasta_Un_file)
	else :
		id = 0
		print >> sys.stdout , "-- Unplaced sequences"
		for comp in all_unused :
			id += 1
			Obj_name = options.prefix + "_Un_" + str(id)
			make_fasta_from_list( [ [comp+"|+", 0 ,0 ,0 ,0] ] , query , 0 , Obj_name , fasta_Un_file)
	fasta_Un_file.close()

### Finished

print >> sys.stdout , "------------------------------"
print >> sys.stdout , "- Done"
print >> sys.stdout , "------------------------------"

