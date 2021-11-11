#!/usr/bin/env python

import argparse
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *

gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


def main():
	#### Main

	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fasta", dest="fasta", nargs='+',
					help="FASTA file(s) of sequences to analyse. [Required]", metavar="seq.fasta")

	parser.add_argument("-g", "--genomesize" , dest="genome_size",
					help="Expected genome size. If given, it is used to calculate NG values", metavar="N")
	parser.add_argument("-s", "--nstep" , dest="nstep", default="10",
					help="Step between N values reported [Default: 10]", metavar="N")
	parser.add_argument("--ngstep" , dest="ngstep", default="1",
					help="Step between NG values reported [Default: 1]", metavar="N")
	parser.add_argument("-c", "--contigs", dest="contigs_stats", default=True, action="store_false",
					help="Avoid performing contigs identification and statistics [Default: perform]")

	parser.add_argument("--nolengths", dest="print_len", default=True, action="store_false",
					help="Avoid printing sequence length file [Default: print]")
	parser.add_argument("--printcontigs", dest="print_contigs", default=False, action="store_true",
					help="Print contigs sequences in [input_fasta].contigs.fasta [Default: do not print]")

	parser.add_argument("--gap", dest="gap_size", default="10",
					help="Minimum size of a Ns stretch to be considered a gap between contigs [Default: 10bp]", metavar="N")

	if len(sys.argv) < 3:
		parser.print_help()
		sys.exit(1)

	options = parser.parse_args()
	gap_size = int(options.gap_size)

	if not options.fasta :
		print >> sys.stderr , "[ERROR] Genome FASTA file missing"
		parser.print_help()
		sys.exit(1)

	fasta_files_list = options.fasta

	Nsize_list = range( 0 , 100 , int(options.nstep) )
	NGsize_list = range( 0 , 100 , int(options.ngstep) )

	results_db = {}
	for sequence_file in fasta_files_list :
		print >> sys.stderr, "# Processing " + sequence_file
		print >> sys.stderr, "## Importing sequences"
		fasta_in = read_fasta( sequence_file )
		reference = {}
		for id in fasta_in :
			reference[id] = remove_islets_fasta(fasta_in[id].upper())
		results_db[sequence_file] = {}

		print >> sys.stderr, "## Analysing base content"
		results_db[sequence_file]["base_content"] = gc_count( reference )

		print >> sys.stderr, "## Analysing input sequences"
		genome_seq_lengths = sorted([ len(reference[x]) for x in reference ] ,reverse=True)
		if options.print_len :
			print >> sys.stderr, "## Printing sequences lengths in: " + sequence_file + ".len"
			try :
				out_file = open( sequence_file + ".len" ,'w')
				for seq_id in sorted(reference.keys()) :
					print >> out_file, seq_id + "\t" + str(len(reference[seq_id]))
				out_file.close()
			except IOError , e :
				print >> sys.stderr, "[ERROR] Impossible to write on " + sequence_file + ".len" + "(Error " + str(e[0]) + ": " + e[1] + ")"

		results_db[sequence_file]["len_stats"] = len_stats(genome_seq_lengths)

		results_db[sequence_file]["Nvalue_Length"] = {}
		results_db[sequence_file]["Nvalue_Index"] = {}
		for threshold in Nsize_list :
			N , L = Nvalue( int(threshold) , genome_seq_lengths)
			results_db[sequence_file]["Nvalue_Length"][threshold] = str('%.0f' % N)
			results_db[sequence_file]["Nvalue_Index"][threshold] = str('%.0f' % L)

		# If there the expected genome size is given, calculate also NG values
		if options.genome_size :
			results_db[sequence_file]["NGvalue_Length"] = {}
			results_db[sequence_file]["NGvalue_Index"] = {}
			for threshold in NGsize_list :
				N , L = NGvalue( int(threshold) , genome_seq_lengths , int(options.genome_size) )
				results_db[sequence_file]["NGvalue_Length"][threshold] = str('%.0f' % N)
				results_db[sequence_file]["NGvalue_Index"][threshold] = str('%.0f' % L)

		# Get gaps
		print >> sys.stderr, "## Analysing gaps in input sequences"
		gaps_in = get_gap_from_fasta_db(reference)
		# format:
		# gaps_in == gap_db[seqid] = [ [ id , start , stop ] , .. , [ ] ]
		results_db[sequence_file]["Gaps_stats"] = gap_stats( gaps_in , sum(genome_seq_lengths) , int(options.gap_size) )

		if options.contigs_stats :
			print >> sys.stderr, "## Extracting contigs"
			contigs_in = split_fasta_on_gap( reference , gaps_in , int(options.gap_size) )
			print >> sys.stderr, "## Analysing contigs sequences"
			contigs_seq_lengths = sorted([ len(contigs_in[x]) for x in contigs_in ] ,reverse=True)
			results_db[sequence_file]["contig_len_stats"] = len_stats(contigs_seq_lengths)

			results_db[sequence_file]["contig_Nvalue_Length"] = {}
			results_db[sequence_file]["contig_Nvalue_Index"] = {}
			for threshold in Nsize_list :
				N , L = Nvalue( int(threshold) , contigs_seq_lengths)
				results_db[sequence_file]["contig_Nvalue_Length"][threshold] = str('%.0f' % N)
				results_db[sequence_file]["contig_Nvalue_Index"][threshold] = str('%.0f' % L)

			# If there the expected genome size is given, calculate also NG values
			if options.genome_size :
				results_db[sequence_file]["contig_NGvalue_Length"] = {}
				results_db[sequence_file]["contig_NGvalue_Index"] = {}
				for threshold in NGsize_list :
					N , L = NGvalue( int(threshold) , contigs_seq_lengths , int(options.genome_size) )
					results_db[sequence_file]["contig_NGvalue_Length"][threshold] = str('%.0f' % N)
					results_db[sequence_file]["contig_NGvalue_Index"][threshold] = str('%.0f' % L)

			# Print contigs
			if options.print_contigs :
				print >> sys.stderr, "## Printing contigs sequences in: " + sequence_file + ".contigs.fasta"
				try :
					write_fasta_from_db(contigs_in , sequence_file + ".contigs.fasta" , False)
				except IOError , e :
					print >> sys.stderr, "[ERROR] Impossible to write on " + sequence_file + ".contigs.fasta" + "(Error " + str(e[0]) + ": " + e[1] + ")"

	print >> sys.stderr, "# Printing statistics"

	# print ids
	id_list = fasta_files_list
	print >> sys.stdout, "# Input sequences statistics"
	print >> sys.stdout , "File name" + "\t" + "\t".join([ str(x) for x in id_list ] )

	print >> sys.stdout , "Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Cumulative length" ) )
	print >> sys.stdout , "Number of sequences" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Number of sequences" ) )
	print >> sys.stdout , "Average sequence length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Average sequence length" ) )
	print >> sys.stdout , "Median sequence length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Median sequence length" ) )
	print >> sys.stdout , "Minimum sequence length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Minimum sequence length" ) )
	print >> sys.stdout , "Maximum sequence length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Maximum sequence length" ) )

	print >> sys.stdout, "# Sequence count by sequence size"
	print >> sys.stdout , "Sequences > 100b - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100b - Count" ) )
	print >> sys.stdout , "Sequences > 500b - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500b - Count" ) )
	print >> sys.stdout , "Sequences > 1Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Kb - Count" ) )
	print >> sys.stdout , "Sequences > 5Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Kb - Count" ) )
	print >> sys.stdout , "Sequences > 10Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 10Kb - Count" ) )
	print >> sys.stdout , "Sequences > 50Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 50Kb - Count" ) )
	print >> sys.stdout , "Sequences > 100Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100Kb - Count" ) )
	print >> sys.stdout , "Sequences > 500Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500Kb - Count" ) )
	print >> sys.stdout , "Sequences > 1Mb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Mb - Count" ) )
	print >> sys.stdout , "Sequences > 5Mb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Mb - Count" ) )

	print >> sys.stdout, "# Sequence cumulative length by sequence size"
	print >> sys.stdout , "Sequences > 100b - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100b - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 500b - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500b - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 1Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Kb - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 5Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Kb - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 10Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 10Kb - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 50Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 50Kb - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 100Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100Kb - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 500Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500Kb - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 1Mb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Mb - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 5Mb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Mb - Cumulative length" ) )

	print >> sys.stdout, "# Sequence count fraction by sequence size"
	print >> sys.stdout , "Sequences > 100b ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100b ratio - Count" ) )
	print >> sys.stdout , "Sequences > 500b ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500b ratio - Count" ) )
	print >> sys.stdout , "Sequences > 1Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Kb ratio - Count" ) )
	print >> sys.stdout , "Sequences > 5Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Kb ratio - Count" ) )
	print >> sys.stdout , "Sequences > 10Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 10Kb ratio - Count" ) )
	print >> sys.stdout , "Sequences > 50Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 50Kb ratio - Count" ) )
	print >> sys.stdout , "Sequences > 100Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100Kb ratio - Count" ) )
	print >> sys.stdout , "Sequences > 500Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500Kb ratio - Count" ) )
	print >> sys.stdout , "Sequences > 1Mb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Mb ratio - Count" ) )
	print >> sys.stdout , "Sequences > 5Mb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Mb ratio - Count" ) )

	print >> sys.stdout, "# Sequence cumulative length fraction by sequence size"
	print >> sys.stdout , "Sequences > 100b ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100b ratio - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 500b ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500b ratio - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 1Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Kb ratio - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 5Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Kb ratio - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 10Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 10Kb ratio - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 50Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 50Kb ratio - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 100Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100Kb ratio - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 500Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500Kb ratio - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 1Mb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Mb ratio - Cumulative length" ) )
	print >> sys.stdout , "Sequences > 5Mb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Mb ratio - Cumulative length" ) )

	# N values, length and index
	print >> sys.stdout, "# Sequence N-values length"
	for value in Nsize_list :
		# Print N length
		print >> sys.stdout , "N" + str(value) + " Length\t" + "\t".join( list_from_stats_db( results_db , id_list , "Nvalue_Length" , value ) )
	print >> sys.stdout, "# Sequence N-values index"
	for value in Nsize_list :
		# Print N index
		print >> sys.stdout , "N" + str(value) + " Index\t" + "\t".join( list_from_stats_db( results_db , id_list , "Nvalue_Index" , value ) )
		
	# NG values, length and index
	if options.genome_size :
		print >> sys.stdout, "# Sequence NG-values length (based on and expected genome size of " + options.genome_size + "bp)"
		for value in NGsize_list :
			# Print NG length
			print >> sys.stdout , "NG" + str(value) + " Length\t" + "\t".join( list_from_stats_db( results_db , id_list , "NGvalue_Length" , value ) )
		print >> sys.stdout, "# Sequence NG-values index (based on and expected genome size of " + options.genome_size + "bp)"
		for value in NGsize_list :
			# Print NG index
			print >> sys.stdout , "NG" + str(value) + " Index\t" + "\t".join( list_from_stats_db( results_db , id_list , "NGvalue_Index" , value ) )

	# Composition
	print >> sys.stdout, "# Sequences nucleotide content"
	print >> sys.stdout ,"A count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "A_count" ) )
	print >> sys.stdout ,"C count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "C_count" ) )
	print >> sys.stdout ,"T count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "T_count" ) )
	print >> sys.stdout ,"G count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "G_count" ) )
	print >> sys.stdout ,"U count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "U_count" ) )
	print >> sys.stdout ,"N count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "N_count" ) )
	print >> sys.stdout ,"Other bases count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "Other_count" ) )
	print >> sys.stdout ,"GC percentage" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "GC_perc" ) )
	print >> sys.stdout ,"GC percentage without unknown bases" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "GC_perc_corr" ) )

	print >> sys.stdout, "##############################"
	print >> sys.stdout, "# Gaps statistics"
	print >> sys.stdout , "File name" + "\t" + "\t".join([ str(x) for x in id_list ] )
	# Gaps
	print >> sys.stdout ,"Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Cumulative length" ) )
	print >> sys.stdout ,"Percentage of assembly length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Genome percentage" ) )
	print >> sys.stdout ,"Number of gaps" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Number of sequences" ) )
	print >> sys.stdout ,"Average gap length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Average sequence length" ) )
	print >> sys.stdout ,"Median gap length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Median sequence length" ) )
	print >> sys.stdout ,"Minimum gap length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Minimum sequence length" ) )
	print >> sys.stdout ,"Maximum gap length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Maximum sequence length" ) )

	print >> sys.stdout, "# Gap count by sequence size"
	print >> sys.stdout ,"Gaps > 100b - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 100b - Count" ) )
	print >> sys.stdout ,"Gaps > 500b - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 500b - Count" ) )
	print >> sys.stdout ,"Gaps > 1Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 1Kb - Count" ) )
	print >> sys.stdout ,"Gaps > 5Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 5Kb - Count" ) )
	print >> sys.stdout ,"Gaps > 10Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 10Kb - Count" ) )
	print >> sys.stdout ,"Gaps > 50Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 50Kb - Count" ) )
	print >> sys.stdout ,"Gaps > 100Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 100Kb - Count" ) )
	print >> sys.stdout ,"Gaps > 500Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 500Kb - Count" ) )
	print >> sys.stdout ,"Gaps > 1Mb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 1Mb - Count" ) )
	print >> sys.stdout ,"Gaps > 5Mb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 5Mb - Count" ) )

	print >> sys.stdout, "# Gap cumulative length by sequence size"
	print >> sys.stdout ,"Gaps > 100b - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 100b - Cumulative length" ) )
	print >> sys.stdout ,"Gaps > 500b - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 500b - Cumulative length" ) )
	print >> sys.stdout ,"Gaps > 1Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 1Kb - Cumulative length" ) )
	print >> sys.stdout ,"Gaps > 5Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 5Kb - Cumulative length" ) )
	print >> sys.stdout ,"Gaps > 10Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 10Kb - Cumulative length" ) )
	print >> sys.stdout ,"Gaps > 50Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 50Kb - Cumulative length" ) )
	print >> sys.stdout ,"Gaps > 100Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 100Kb - Cumulative length" ) )
	print >> sys.stdout ,"Gaps > 500Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 500Kb - Cumulative length" ) )
	print >> sys.stdout ,"Gaps > 1Mb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 1Mb - Cumulative length" ) )
	print >> sys.stdout ,"Gaps > 5Mb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 5Mb - Cumulative length" ) )

	# Contigs
	if options.contigs_stats :
		print >> sys.stdout, "##############################"
		print >> sys.stdout, "# Contigs statistics"
		print >> sys.stdout , "File name" + "\t" + "\t".join([ str(x) for x in id_list ] )

		print >> sys.stdout , "Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Cumulative length" ) )
		print >> sys.stdout , "Number of contigs" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Number of sequences" ) )
		print >> sys.stdout , "Average contig length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Average sequence length" ) )
		print >> sys.stdout , "Median contig length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Median sequence length" ) )
		print >> sys.stdout , "Minimum contig length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Minimum sequence length" ) )
		print >> sys.stdout , "Maximum contig length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Maximum sequence length" ) )

		print >> sys.stdout, "# Contig count by sequence size"
		print >> sys.stdout , "Contigs > 100b - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100b - Count" ) )
		print >> sys.stdout , "Contigs > 500b - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500b - Count" ) )
		print >> sys.stdout , "Contigs > 1Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Kb - Count" ) )
		print >> sys.stdout , "Contigs > 5Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Kb - Count" ) )
		print >> sys.stdout , "Contigs > 10Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 10Kb - Count" ) )
		print >> sys.stdout , "Contigs > 50Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 50Kb - Count" ) )
		print >> sys.stdout , "Contigs > 100Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100Kb - Count" ) )
		print >> sys.stdout , "Contigs > 500Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500Kb - Count" ) )
		print >> sys.stdout , "Contigs > 1Mb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Mb - Count" ) )
		print >> sys.stdout , "Contigs > 5Mb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Mb - Count" ) )

		print >> sys.stdout, "# Contig cumulative length by sequence size"
		print >> sys.stdout , "Contigs > 100b - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100b - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 500b - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500b - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 1Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Kb - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 5Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Kb - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 10Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 10Kb - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 50Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 50Kb - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 100Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100Kb - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 500Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500Kb - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 1Mb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Mb - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 5Mb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Mb - Cumulative length" ) )

		print >> sys.stdout, "# Contig count ratio by sequence size"
		print >> sys.stdout , "Contigs > 100b ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100b ratio - Count" ) )
		print >> sys.stdout , "Contigs > 500b ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500b ratio - Count" ) )
		print >> sys.stdout , "Contigs > 1Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Kb ratio - Count" ) )
		print >> sys.stdout , "Contigs > 5Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Kb ratio - Count" ) )
		print >> sys.stdout , "Contigs > 10Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 10Kb ratio - Count" ) )
		print >> sys.stdout , "Contigs > 50Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 50Kb ratio - Count" ) )
		print >> sys.stdout , "Contigs > 100Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100Kb ratio - Count" ) )
		print >> sys.stdout , "Contigs > 500Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500Kb ratio - Count" ) )
		print >> sys.stdout , "Contigs > 1Mb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Mb ratio - Count" ) )
		print >> sys.stdout , "Contigs > 5Mb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Mb ratio - Count" ) )

		print >> sys.stdout, "# Contig cumulative length ratio by sequence size"
		print >> sys.stdout , "Contigs > 100b ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100b ratio - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 500b ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500b ratio - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 1Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Kb ratio - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 5Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Kb ratio - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 10Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 10Kb ratio - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 50Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 50Kb ratio - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 100Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100Kb ratio - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 500Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500Kb ratio - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 1Mb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Mb ratio - Cumulative length" ) )
		print >> sys.stdout , "Contigs > 5Mb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Mb ratio - Cumulative length" ) )

		print >> sys.stdout, "# Contig N-values length"
		# N values, length and index
		for value in Nsize_list :
			# Print N length
			print >> sys.stdout , "N" + str(value) + " Length\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_Nvalue_Length" , value ) )
		print >> sys.stdout, "# Contig N-values index"
		for value in Nsize_list :
			# Print N index
			print >> sys.stdout , "N" + str(value) + " Index\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_Nvalue_Index" , value ) )
			
		# NG values, length and index
		if options.genome_size :
			print >> sys.stdout, "# Contig NG-values length (based on and expected genome size of " + options.genome_size + "bp)"
			for value in NGsize_list :
				# Print NG length
				print >> sys.stdout , "NG" + str(value) + " Length\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_NGvalue_Length" , value ) )
			print >> sys.stdout, "# Contig NG-values index (based on and expected genome size of " + options.genome_size + "bp)"
			for value in NGsize_list :
				# Print NG index
				print >> sys.stdout , "NG" + str(value) + " Index\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_NGvalue_Index" , value ) )


if __name__ == '__main__':
	main()
