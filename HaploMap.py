#!/usr/bin/env python

import argparse
from lib_files.HaploFunct import *
from lib_files.map_lib import *


def main() :

	###### Options and help ######
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fasta", default=False, dest="reference",
					help="FASTA file(s) with the genomic sequences to pair", metavar="genome.fasta [Required]")
	parser.add_argument("-c", "--sequence_correspondence", dest="corr",
					help="Tab text file(s) of corresponding sequence names in the two haplotypes. [REQUIRED]", metavar="Hap1_to_Hap2.txt")
	parser.add_argument("--tmpdir", dest="tmp", default="tmp_HaploMap" ,
					help="Path to temporary directory. [Default: ./tmp_HaploMap]", metavar="/path/to/tmp")
	parser.add_argument("-m" , "--mapper", dest="mapper", default="nucmer",
					help="Mapping tool to use [Default: nucmer]" , metavar="[minimap|nucmer|blat]")
	parser.add_argument("-t", "--threads", dest="cores", default=4,
					help="Cores used in mapping process [default: 4]", metavar="N")
	parser.add_argument("-o", "--out", dest="out", default="out",
					help="Output files prefix [default: out]", metavar="NAME")

	print >> sys.stdout, "Running HaploMaker tool from HaploSync version " + get_version()
	print >> sys.stdout, "To reproduce this run use the following command: " + " ".join( pipes.quote(x) for x in sys.argv)
	print >> sys.stdout, "----"
	scriptDirectory = os.path.dirname(os.path.realpath(__file__)) + "/support_scripts"
	# Sanity Check

	if len(sys.argv) < 3:
		parser.print_help()
		sys.exit(1)

	options = parser.parse_args()

	if not ( options.hap1 and options.hap2 ) :
		print >> sys.stderr , "[ERROR] FASTA file missing"
		parser.print_help()
		sys.exit(1)

	if not options.corr :
		print >> sys.stderr , "[ERROR] Sequence name correspondence file missing"
		parser.print_help()
		sys.exit(1)

	paths = set_paths(os.path.join(sys.path[0], 'HaploSync.conf.toml'))
	minimap_path = paths["minimap2"]
	samtools_path = paths["samtools"]
	bedtools_path = paths["bedtools"]
	nucmer_path = paths["nucmer"]
	showcoords_path = paths["show-coords"]
	blat_path = paths["blat"]

	temp_folder = options.tmp
	if not os.path.exists( options.tmp ) :
		mkdir(temp_folder)

	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] = Read inputs'
	print >> sys.stderr , '# Read inputs'
	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Loading FASTA sequences'
	print >> sys.stderr , '## Loading FASTA sequences'
	fasta_files = options.fasta
	fasta_files_list = fasta_files.split(",")
	fasta_dict = {}
	for file_name in fasta_files_list :
		if fasta_dict == {} :
			fasta_dict = read_fasta(file_name)
		else :
			fasta_dict.update(read_fasta(file_name))
	fasta_len_dict = get_length_from_fasta_db(fasta_dict)

	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Loading sequence correspondence information'
	print >> sys.stderr , '## Loading sequence correspondence information'
	corr_files = options.corr
	corr_files_list = corr_files.split(",")
	pairs = []
	for file in corr_files_list :
		for line in open(file) :
			seq1 , seq2 = line.rstrip().split("\t")
			pairs.append([seq1 , seq2])

	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Paring sequences'
	print >> sys.stderr , '## Paring sequences'
	uniq_alignment_file = {}
	for seq_pair in sorted(pairs) :
		print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] === Paring ' + seq_pair[0] + " and " + seq_pair[1]
		print >> sys.stderr , '### Paring ' + seq_pair[0] + " and " + seq_pair[1]
		print >> sys.stderr , '#### Mapping '
		alignment_file_prefix = seq_pair[1] + ".on." + seq_pair[0]
		seq1_id = seq_pair[0]
		seq1_fasta = fasta_dict(seq_pair[0])
		seq2_id = seq_pair[1]
		seq2_fasta = fasta_dict(seq_pair[1])
		all_hits = map_sequences( seq1_id , seq1_fasta , seq2_id , seq2_fasta , options.mapper , int(options.cores) , temp_folder , alignment_file_prefix , paths  )
		print >> sys.stderr , '#### Finding best tiling path'
		uniq_alignment_file[( seq1_id , seq2_id )] = hits_best_tiling_path(all_hits, fasta_len_dict)[( seq1_id , seq2_id )]
	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Exporting maps'
	print >> sys.stderr , '## Exporting maps'
	out_file_name = options.out + ".paired_regions.txt"
	out_file = open(out_file_name , 'w')
	for pair in sorted(uniq_alignment_file.keys()) :
		for hit in sorted(uniq_alignment_file[pair]) :
			print >> out_file , "\t".join([ str(x) for x in hit ])

	###### END ######

	print >> sys.stdout , "------------------------------"
	print >> sys.stdout , "- Done"
	print >> sys.stdout , "------------------------------"
	print >> sys.stderr , "##############################"
	print >> sys.stderr , "# Done"
	print >> sys.stderr , "##############################"



if __name__ == '__main__':
	main()