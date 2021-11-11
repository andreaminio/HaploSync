#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
import gzip
import sys
import re
import numpy as np



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


def get_fasta_lengths_from_file( fasta_file ) :
	length_db = {}

	for record in SeqIO.parse(open(fasta_file), 'fasta'):
		length_db[record.id] = len(record)

	return length_db


def get_fasta_ids( fasta_file ) :
	ids_list = []

	for line in open( fasta_file , 'r' ) :
		if line[0] == ">" :
			seq_id = line.rstrip()[1:].split()[0]
			ids_list.append(seq_id)
		else :
			continue

	return sorted(ids_list)


def get_gap_from_fasta_db( fasta_dict , mask_ranges = {} ) :
	gap_db = {}
	for id in fasta_dict :
		gap_db[id] = []
		for match in re.finditer('N+', fasta_dict[id]):
			start = int(match.start())
			stop = int(match.end())
			to_mask = False
			if id in mask_ranges :
				for range in mask_ranges[id] :
					try :
						int(range[1])
					except :
						print >> sys.stderr, range
					if ( int(range[1]) <= start <= int(range[2]) ) or ( int(range[1]) <= stop <= int(range[2]) ):
						to_mask = True
						break
			if not to_mask :
				gap_db[id].append([ id , start , stop ])

	return gap_db


def get_length_from_fasta_db(fasta_dict) :
	sequences_len = {}

	for ref_name in fasta_dict :
		sequences_len[ref_name] = int(len(fasta_dict[ref_name]))

	return sequences_len


def split_fasta(fasta_db , folder) :
	fasta_files_db = {}

	for seq_id in fasta_db :
		filename = folder + "/" + seq_id + ".fasta"
		out_fasta_file = open(filename , 'w')
		print >> out_fasta_file , ">" + str(seq_id)
		seq_len = len(fasta_db[seq_id])
		for i in range( 0 , seq_len , 60 ):
			print >> out_fasta_file , fasta_db[seq_id][i : min( i + 60 , seq_len) ]
		out_fasta_file.close()
		fasta_files_db[seq_id] = filename

	return fasta_files_db


def remove_islets_fasta( sequence_string ) :
	new_string = ""
	prev_start = 0
	end = len(sequence_string)

	for match in re.finditer('N+', sequence_string) :
		gap_start = int(match.start())
		gap_stop = int(match.end())
		gap_seq = sequence_string[gap_start:gap_stop]
		prev_stop = gap_start
		if prev_start == 0 or (prev_stop - prev_start) >= 10 :
			chunk_string = sequence_string[prev_start:prev_stop]
		else :
			chunk_string = "N"*(prev_stop-prev_start)
		new_string += chunk_string
		new_string += gap_seq
		prev_start = gap_stop

	new_string += sequence_string[prev_start:end]

	return new_string


def get_block_extremities(fasta_db, fasta_length_db, chr, start, stop, orientation, gene_position) :
	region_start = int(start)
	region_stop = int(stop)
	#if not start < 5 :
	#	# Region do not start from the beginning of the sequence
	#	# Evaluate possibility to add overhang
	#	# 	try to extend upstream up to overhang
	#	#	stop the extension as a stretch of 5 Ns is found
	#	delta = 4
	#	max_delta = min( start - 5 , overhang - 5 )
	#	while delta < max_delta :
	#		delta += 1
	#		upstream_base = fasta_db[chr][ start - delta : start - delta + 5]
	#		if str(upstream_base).upper() == "NNNNN":
	#			delta = delta - 5
	#			break
	#	region_start = start - delta
	#else :
	#	region_start = 0
	#if not stop > int(fasta_length_db[chr]) - 5 :
	#	# Region do not stop at the end of the sequence
	#	# Evaluate possibility to add overhang
	#	# 	try to extend downstream up to overhang
	#	#	stop the extension as a stretch of 5 Ns is found
	#	delta = 4
	#	max_delta = min( fasta_length_db[chr] - stop - 5 , overhang - 5 )
	#	while delta < max_delta :
	#		delta += 1
	#		upstream_base = fasta_db[chr][ stop + delta - 5 : stop + delta]
	#		if str(upstream_base).upper() == "NNNNN":
	#			delta = delta - 5
	#			break
	#	region_stop = stop + delta
	#else :
	#	region_stop = fasta_length_db[chr]

	# Extract genes in the range
	matching_genes = []
	gene_list = []
	if chr in gene_position :
		do_search_genes = True
		while do_search_genes :
			range_extended = False
			for gene in gene_position[chr] :

				# gene_position format:
				# gene_position[chr] = [
				#	...
				# 	[chr, int(start) , int(end) , attributes_dict["ID"]]
				#	...
				#	]
				if region_start <= gene[1] <= region_stop or region_start <= gene[2] <= region_stop :
					if gene not in matching_genes :
						matching_genes.append(gene)
					# if the gene is partially outside of the region, expand to cover it entirely
					if gene[1] < region_start :
						region_start = gene[1]
						range_extended = True
					if gene[2] > region_stop :
						region_stop = gene[2]
						range_extended = True

			if not range_extended :
				do_search_genes=False

	# Extract sequence
	print >> sys.stderr , "### Generating region. Region required: " + chr + ":" + str(start) + ":" + str(stop) + "(" + orientation + ") | Region in use: " + chr + ":" + str(region_start) + ":" + str( region_stop) + "(" + orientation + ")"
	coordinates = [ chr , region_start , region_stop , orientation ]
	region_sequence = fasta_db[chr][region_start : region_stop]
	if orientation == "+" :
		sequence = region_sequence
	else :
		sequence = str( Seq( region_sequence ).reverse_complement() )

	# Convert annotation to region
	for gene in sorted(matching_genes) :
		if orientation == "+" :
			gene_start = gene[1] - region_start
			gene_stop = gene[2] - region_start
			gene_id = gene[3]
			gene_list.append([ gene_start , gene_stop , gene_id ])
			#print >> sys.stderr , "#### " +  gene_id + " | " + str(gene[1]) + ":" + str(gene[2])
		else :
			# Reverse complement
			gene_start = region_stop - gene[2]
			gene_stop = region_stop - gene[1]
			gene_id = gene[3]
			gene_list.append([ gene_start , gene_stop , gene_id ])
			#print >> sys.stderr , "#### " + gene_id + " | " + str(gene[1]) + ":" + str(gene[2]) + " >>> " + str(gene_start) + ":" + str(gene_stop)

	return coordinates , sequence , gene_list


def make_fasta_from_list( querylist , queryfasta , gaplen , seqoutname ) :
	### Query list element format
	# Gap:	[61252	,	(0:61252)		,	gap		,	61252	,	0]
	#		[length	,	(T_start:Tstop)	, 	"gap" 	, 	length 	, 	0]
	# Object:	[b40-14.HS_iter4_seq3733|+	,	(61252:6463804)	,	93612:7595148	,	-6402552			,	4526208]
	#			[ID|strand					,	(T_start:Tstop)	,	Q_start:Q_stop	,	-(alignment length)	,	matches]
	gaplen = int(gaplen)
	seq = ""

	for Id in querylist :
		CompntId_name = Id[:-2]
		Orientation = Id[-1]
		# Add gap between Components
		if str(seq) != "" :
			seq = seq + "N" * gaplen
		if Orientation == "-" :
			my_sub_seq = str(Seq(queryfasta[CompntId_name]).reverse_complement())
		else :
			my_sub_seq = queryfasta[CompntId_name]
		seq = seq + my_sub_seq

	return seq


def Nvalue( n , len_list):
	curr_sum = 0
	count = 0
	tot_len = sum(len_list)
	for length in sorted(len_list, reverse=True):
		curr_sum += length
		count += 1
		if curr_sum >= float(tot_len*n)/float(100):
			return length , count
	return 0 , 0


def NGvalue( n , len_list , genome_size):
	curr_sum = 0
	count = 0
	for length in sorted(len_list, reverse=True):
		curr_sum += length
		count += 1
		if curr_sum >= float(genome_size*n)/float(100) :
			return length , count
	return 0 , 0



def len_stats( len_list ) :
	lengths = {}
	# generate a DB with all the length descriptive statistics
	# save the values as formatted string ready to print
	total_sum = sum(len_list)
	total_count = len(len_list)
	lengths["Cumulative length"] = str('%.0f' % total_sum)
	lengths["Number of sequences"] = str('%.0f' % total_count)
	lengths["Average sequence length"] = str('%.2f' % np.mean(len_list))
	lengths["Median sequence length"] = str('%.0f' % np.median(len_list))
	lengths["Minimum sequence length"] = str('%.0f' % np.min(len_list))
	lengths["Maximum sequence length"] = str('%.0f' % np.max(len_list))

	lengths["Sequences > 100b - Count" ] = str('%.0f' % len([x for x in len_list if x > 100 ]))
	lengths["Sequences > 500b - Count" ] = str('%.0f' % len([x for x in len_list if x > 500 ])) 
	lengths["Sequences > 1Kb - Count" ] = str('%.0f' % len([x for x in len_list if x > 1000 ])) 
	lengths["Sequences > 5Kb - Count" ] = str('%.0f' % len([x for x in len_list if x > 5000 ])) 
	lengths["Sequences > 10Kb - Count" ] = str('%.0f' % len([x for x in len_list if x > 10000 ])) 
	lengths["Sequences > 50Kb - Count" ] = str('%.0f' % len([x for x in len_list if x > 50000 ])) 
	lengths["Sequences > 100Kb - Count" ] = str('%.0f' % len([x for x in len_list if x > 100000 ])) 
	lengths["Sequences > 500Kb - Count" ] = str('%.0f' % len([x for x in len_list if x > 500000 ])) 
	lengths["Sequences > 1Mb - Count" ] = str('%.0f' % len([x for x in len_list if x > 1000000 ])) 
	lengths["Sequences > 5Mb - Count" ] = str('%.0f' % len([x for x in len_list if x > 5000000 ])) 
	lengths["Sequences > 100b - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 100 ]))
	lengths["Sequences > 500b - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 500 ]))
	lengths["Sequences > 1Kb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 1000 ]))
	lengths["Sequences > 5Kb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 5000 ]))
	lengths["Sequences > 10Kb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 10000 ]))
	lengths["Sequences > 50Kb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 50000 ]))
	lengths["Sequences > 100Kb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 100000 ]))
	lengths["Sequences > 500Kb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 500000 ]))
	lengths["Sequences > 1Mb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 1000000 ]))
	lengths["Sequences > 5Mb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 5000000 ]))

	lengths["Sequences > 100b ratio - Count" ] = str('%.2f' % ( float(100 * len([x for x in len_list if x > 100 ]) ) / float(total_count) ) )
	lengths["Sequences > 500b ratio - Count" ] = str('%.2f' % ( float(100 * len([x for x in len_list if x > 500 ]) ) / float(total_count) ) )
	lengths["Sequences > 1Kb ratio - Count" ] = str('%.2f' % ( float(100 * len([x for x in len_list if x > 1000 ]) ) / float(total_count) ) )
	lengths["Sequences > 5Kb ratio - Count" ] = str('%.2f' % ( float(100 * len([x for x in len_list if x > 5000 ]) ) / float(total_count) ) )
	lengths["Sequences > 10Kb ratio - Count" ] = str('%.2f' % ( float(100 * len([x for x in len_list if x > 10000 ]) ) / float(total_count) ) )
	lengths["Sequences > 50Kb ratio - Count" ] = str('%.2f' % ( float(100 * len([x for x in len_list if x > 50000 ]) ) / float(total_count) ) )
	lengths["Sequences > 100Kb ratio - Count" ] = str('%.2f' % ( float(100 * len([x for x in len_list if x > 100000 ]) ) / float(total_count) ) )
	lengths["Sequences > 500Kb ratio - Count" ] = str('%.2f' % ( float(100 * len([x for x in len_list if x > 500000 ]) ) / float(total_count) ) )
	lengths["Sequences > 1Mb ratio - Count" ] = str('%.2f' % ( float(100 * len([x for x in len_list if x > 1000000 ]) ) / float(total_count) ) )
	lengths["Sequences > 5Mb ratio - Count" ] = str('%.2f' % ( float(100 * len([x for x in len_list if x > 5000000 ]) ) / float(total_count) ) )
	lengths["Sequences > 100b ratio - Cumulative length"] = str('%.2f' % ( float(100 * sum([x for x in len_list if x > 100 ]) ) / float(total_sum) ) )
	lengths["Sequences > 500b ratio - Cumulative length"] = str('%.2f' % ( float(100 * sum([x for x in len_list if x > 500 ]) ) / float(total_sum) ) )
	lengths["Sequences > 1Kb ratio - Cumulative length"] = str('%.2f' % ( float(100 * sum([x for x in len_list if x > 1000 ]) ) / float(total_sum) ) )
	lengths["Sequences > 5Kb ratio - Cumulative length"] = str('%.2f' % ( float(100 * sum([x for x in len_list if x > 5000 ]) ) / float(total_sum) ) )
	lengths["Sequences > 10Kb ratio - Cumulative length"] = str('%.2f' % ( float(100 * sum([x for x in len_list if x > 10000 ]) ) / float(total_sum) ) )
	lengths["Sequences > 50Kb ratio - Cumulative length"] = str('%.2f' % ( float(100 * sum([x for x in len_list if x > 50000 ]) ) / float(total_sum) ) )
	lengths["Sequences > 100Kb ratio - Cumulative length"] = str('%.2f' % ( float(100 * sum([x for x in len_list if x > 100000 ]) ) / float(total_sum) ) )
	lengths["Sequences > 500Kb ratio - Cumulative length"] = str('%.2f' % ( float(100 * sum([x for x in len_list if x > 500000 ]) ) / float(total_sum) ) )
	lengths["Sequences > 1Mb ratio - Cumulative length"] = str('%.2f' % ( float(100 * sum([x for x in len_list if x > 1000000 ]) ) / float(total_sum) ) )
	lengths["Sequences > 5Mb ratio - Cumulative length"] = str('%.2f' % ( float(100 * sum([x for x in len_list if x > 5000000 ]) ) / float(total_sum) ) )

	return lengths


def gc_count( fasta_db ) :
	# Report values as formatted text
	global_counts = { "A_count":0 , "C_count":0 , "T_count":0 , "G_count":0 , "U_count":0 , "N_count":0 , "Other_count":0 , "GC_perc":0 , "GC_perc_corr":0}
	global_total = 0
	for seq_id in fasta_db :
		actual_seq = str(fasta_db[seq_id]).upper()
		total = len(actual_seq)
		A_count = actual_seq.count("A")
		C_count = actual_seq.count("C")
		T_count = actual_seq.count("T")
		G_count = actual_seq.count("G")
		U_count = actual_seq.count("U")
		N_count = actual_seq.count("N")
		Other_count = total - (A_count + C_count + T_count + G_count + U_count + N_count)
		global_total += total
		global_counts["A_count"] += A_count
		global_counts["C_count"] += C_count
		global_counts["T_count"] += T_count
		global_counts["G_count"] += G_count
		global_counts["U_count"] += U_count
		global_counts["N_count"] += N_count
		global_counts["Other_count"] += Other_count

	counts = {}
	counts["A_count"]	   = str('%.0f' % global_counts["A_count"] )
	counts["C_count"]	   = str('%.0f' % global_counts["C_count"] )
	counts["T_count"]	   = str('%.0f' % global_counts["T_count"] )
	counts["G_count"]	   = str('%.0f' % global_counts["G_count"] )
	counts["U_count"]	   = str('%.0f' % global_counts["U_count"] )
	counts["N_count"]	   = str('%.0f' % global_counts["N_count"] )
	counts["Other_count"]  = str('%.0f' % global_counts["Other_count"] )
	counts["GC_perc"]      = str('%.2f' % ( float( 100*( global_counts["C_count"] + global_counts["G_count"] ) ) / float(global_total) ) )
	counts["GC_perc_corr"] = str('%.2f' % ( float( 100*( global_counts["C_count"] + global_counts["G_count"] ) ) / float( global_total - ( global_counts["N_count"] + global_counts["Other_count"] ) ) ) )

	return counts


def gap_stats( gap_db , sequences_length , min_gap_size ) :
	# gap_db[seq_id] = [ [ id , start , stop ] , .. , [ ] ]
	len_list = []
	for seq_id in gap_db.keys() :
		for element in gap_db[seq_id] :
			gap_len = int(element[2]) - int(element[1])
			if gap_len >= min_gap_size :
				len_list.append(gap_len)

	lengths = {}
	total_sum = sum(len_list)
	total_count = len(len_list)
	lengths["Cumulative length"] = str('%.0f' % total_sum)
	lengths["Genome percentage"] = str('%.0f' % ( float(100 * total_sum) / float(sequences_length) ) )
	lengths["Number of sequences"] = str('%.0f' % total_count)
	if len_list == [] :
		lengths["Average sequence length"] = str('%.2f' % 0)
		lengths["Median sequence length"] = str('%.0f' % 0)
		lengths["Minimum sequence length"] = str('%.0f' % 0)
		lengths["Maximum sequence length"] = str('%.0f' % 0)
	else :
		lengths["Average sequence length"] = str('%.2f' % np.mean(len_list))
		lengths["Median sequence length"] = str('%.0f' % np.median(len_list))
		lengths["Minimum sequence length"] = str('%.0f' % np.min(len_list))
		lengths["Maximum sequence length"] = str('%.0f' % np.max(len_list))

	lengths["Sequences > 100b - Count" ] = str('%.0f' % len([x for x in len_list if x > 100 ]))
	lengths["Sequences > 500b - Count" ] = str('%.0f' % len([x for x in len_list if x > 500 ]))
	lengths["Sequences > 1Kb - Count" ] = str('%.0f' % len([x for x in len_list if x > 1000 ]))
	lengths["Sequences > 5Kb - Count" ] = str('%.0f' % len([x for x in len_list if x > 5000 ]))
	lengths["Sequences > 10Kb - Count" ] = str('%.0f' % len([x for x in len_list if x > 10000 ]))
	lengths["Sequences > 50Kb - Count" ] = str('%.0f' % len([x for x in len_list if x > 50000 ]))
	lengths["Sequences > 100Kb - Count" ] = str('%.0f' % len([x for x in len_list if x > 100000 ]))
	lengths["Sequences > 500Kb - Count" ] = str('%.0f' % len([x for x in len_list if x > 500000 ]))
	lengths["Sequences > 1Mb - Count" ] = str('%.0f' % len([x for x in len_list if x > 1000000 ]))
	lengths["Sequences > 5Mb - Count" ] = str('%.0f' % len([x for x in len_list if x > 5000000 ]))
	lengths["Sequences > 100b - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 100 ]))
	lengths["Sequences > 500b - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 500 ]))
	lengths["Sequences > 1Kb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 1000 ]))
	lengths["Sequences > 5Kb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 5000 ]))
	lengths["Sequences > 10Kb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 10000 ]))
	lengths["Sequences > 50Kb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 50000 ]))
	lengths["Sequences > 100Kb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 100000 ]))
	lengths["Sequences > 500Kb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 500000 ]))
	lengths["Sequences > 1Mb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 1000000 ]))
	lengths["Sequences > 5Mb - Cumulative length"] = str('%.0f' % sum([x for x in len_list if x > 5000000 ]))

	return lengths


def split_fasta_on_gap( fasta_db , gap_db , min_gap_size )  :
	# gap_db[seq_id] = [ [ id , start , stop ] , .. , [ ] ]
	contigs_db = {}
	for seq_id in sorted(fasta_db.keys()) :
		#print >> sys.stderr, seq_id
		sequence = fasta_db[seq_id]
		sequence_len = len(sequence)
		gaps_to_break = []

		for element in gap_db[seq_id] :
			gap_len = int(element[2]) - int(element[1])
			if gap_len >= int(min_gap_size) :
				gaps_to_break.append(element)

		#print >> sys.stderr, gaps_to_break
		contigs = reverse_regions( seq_id , gaps_to_break , sequence_len)
		#print >> sys.stderr, contigs
		# Extract new fasta
		component_id = 0
		for contig in sorted(contigs):
			component_id += 1
			contig_id = seq_id + "_Contig" + str(component_id)
			#print >> sys.stderr, seq_id + "_Contig" + str(component_id) + " = " + seq_id + ":" + str(contig[1]) + "-" + str((contig[2]))
			contigs_db[contig_id] = sequence[ int(contig[1]) : int(contig[2]) ]

	return contigs_db


def reverse_regions( seq_id , positive_regions , seq_len ) :
	negative_regions = []
	negative_start = 0
	negative_id = seq_id
	for element in sorted(positive_regions):
		positive_start = int(element[1])
		positive_stop = int(element[2])
		if positive_start == negative_start :
			continue
		else :
			negative_stop = positive_start
			negative_regions.append([ negative_id , negative_start , negative_stop ] )
			negative_start = positive_stop
	if not negative_start == int(seq_len) :
		negative_regions.append([ negative_id , negative_start , int(seq_len) ] )

	return sorted(negative_regions)


def list_from_stats_db( stats_db , id_order , feature , subfeature) :
	#results_db[sequence_file]["len_stats"]["Cumulative length"]
	value_list = []
	for file_id in id_order :
		if (file_id in stats_db) and (feature in stats_db[file_id]) and (subfeature in stats_db[file_id][feature]) :
			value_list.append(stats_db[file_id][feature][subfeature])
		else :
			value_list.append("-")

	return value_list


def make_chunks_from_fasta( seq_name , fasta_sequence , chunk_length=1000) :
	chunk_db = {}
	fasta_len = len(fasta_sequence)
	for start in range(0, fasta_len, chunk_length) :
		stop = min(start + chunk_length , fasta_len)
		chunk_id = seq_name + "#" + str(start) + "#" + str(stop)
		chunk_db[chunk_id] = fasta_sequence[start : stop]
	return chunk_db

