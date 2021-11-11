#!/usr/bin/env python

import argparse
from lib_files.HaploFunct import *
from lib_files.GFF_lib import *
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *

gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


def main() :

	###### Options and help ######
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fasta", dest="fasta",
					help="FASTA file(s) with all input sequences", metavar="genome.fasta [Required]")
	parser.add_argument("-s", "--structure", dest="structure",
					help="Text file describing the novel sequence structure to build. [REQUIRED]", metavar="structure.txt")
	parser.add_argument("--format", dest="format", default="BLOCK" ,
					help="Structure file format [Default: BLOCK]", metavar="[BLOCK|AGP|BED]")

	parser.add_argument("-o", "--out", dest="out", default="out" ,
					help="Output file names [Default: out]", metavar="outname")
	parser.add_argument("-p", "--prefix", dest="prefix", default="" ,
					help="New sequence names prefix to be added before the old name. Leave unset to keep the same names [Default: NONE]", metavar="new")

	parser.add_argument("-g", "--gff3", dest="gff3",
					help="Genome annotation in GFF3 format.", metavar="annotation.gff3")
	parser.add_argument("-b", "--bed", dest="bed",
					help="Region annotation in BED format. Accepted BED3 and BED6 formats, additional columns will be considered as annotation and reported unedited", metavar="annotation.bed")

	parser.add_argument("-a", "--agp", dest="agp",
					help="AGP file defining the coordinates of original sequences composing actual assembly", metavar="previous_to_actual.genome.agp")
	parser.add_argument("-m" , "--mapper", dest="mapper", default="blat",
					#help="Mapping tool to use [Default: blat]" , metavar="[blat|nucmer]")
					help=argparse.SUPPRESS )
	parser.add_argument("-c" , "--cores", dest="cores", default="1",
					#help="[For nucmer] Number of cores to use for mapping [Default: 1]" , metavar="N")
					help=argparse.SUPPRESS )
	parser.add_argument("--skipoverlap", dest="skipoverlap", default=False, action="store_true",
					help="Skip the search of overlap between flanking blocks")
	parser.add_argument("--overhang" , dest="overhang", default="10000" ,
					#help="Maximum overhang sequence to use for overlap search [default 10000]", metavar="N")
					help=argparse.SUPPRESS )
	parser.add_argument("--gap", dest="gap_size", default="1000",
					help="Size of the residual gap around inserted sequences if no overlap is detected between flanking regions of consecutive blocks [Default: 1,000bp]", metavar="N")
	parser.add_argument("--spacer", dest="spacer", default="10" ,
					help="Size of the spacer gap inserted between trimmed sequences when there is an overlap between flanking regions [Default: 10bp]", metavar="N")
	parser.add_argument("--ignoreids", dest="ignoreids", default=False, action="store_true",
					help="Ignore output sequence ids reported in structure file, use [-p|--prefix] with progressive numbers instead [Default for BED input]")

	parser.add_argument("--reverse", dest="mode", default=False, action="store_true",
					help="[For AGP structure files] Reverse direction of feature extraction, from new to old. FASTA input must contain the genomic sequences of the new version")

	parser.add_argument("-u" , "--unplaced" , dest="add_unplaced" , default=False, action="store_true",
					help="Add unplaced sequences into output files" )

	parser.add_argument("--noagp", dest="noagp", default=False, action="store_true",
					help="Avoid printing the AGP file output")
	parser.add_argument("--noprint", dest="noprint", default=False, action="store_true",
					help="Avoid printing the FASTA file output")

	print >> sys.stdout, "Running HaploMake tool from HaploSync version " + get_version()
	print >> sys.stdout, "To reproduce this run use the following command: " + " ".join( pipes.quote(x) for x in sys.argv)
	print >> sys.stdout, "----"

	# Sanity Check

	paths = set_paths(os.path.join(sys.path[0], 'HaploSync.conf.toml'))
	minimap_path = paths["minimap2"]
	samtools_path = paths["samtools"]
	bedtools_path = paths["bedtools"]
	nucmer_path = paths["nucmer"]
	showcoords_path = paths["show-coords"]
	blat_path = paths["blat"]
	scriptDirectory = os.path.dirname(os.path.realpath(__file__)) + "/support_scripts"

	if len(sys.argv) < 3:
		parser.print_help()
		sys.exit(1)

	options = parser.parse_args()
	gap_size = int(options.gap_size)

	if not options.fasta :
		print >> sys.stderr , "[ERROR] Genome FASTA file missing"
		parser.print_help()
		sys.exit(1)

	if not options.structure :
		print >> sys.stderr , "[ERROR] Structure file missing"
		parser.print_help()
		sys.exit(1)

	if (not options.skipoverlap) and options.mode :
		print >> sys.stderr , "[ERROR] Reverse use of AGP structure [--reverse] incompatible with overlap search and correction."
		sys.exit(1)

	if options.noprint :
		no_fasta = True
	else :
		no_fasta = False

	if not options.gff3 :
		print >> sys.stdout, "[MEMO] Annotation file missing: Coordinates conversion will not take place. You'll need to run it downstream of this process"

	#if not options.agp :
	#	print >> sys.stdout, "[MEMO] No AGP file was given. Structure of original sequences will not be available"

	if options.skipoverlap :
		print >> sys.stderr, "[MEMO] Adjacent sequences overlap analysis will be skipped. Sequence will be inserted entirely separated by gaps of " + str(gap_size) + "bp in length.\nThis procedure may allow the use of duplicated genomic content"

	# Read inputs
	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] = Read inputs'
	print >> sys.stderr , '# Read inputs'
	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Loading FASTA sequences'
	print >> sys.stderr , '## Loading FASTA sequences'
	fasta_files = options.fasta
	fasta_files_list = fasta_files.split(",")
	fasta_dict = {}
	for file_name in fasta_files_list :
		print >> sys.stderr , '### Loading ' + file_name
		if fasta_dict == {} :
			fasta_dict = read_fasta(file_name)
		else :
			fasta_dict.update(read_fasta(file_name))
	fasta_len_dict = get_length_from_fasta_db(fasta_dict)

	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Listing structure files'
	print >> sys.stderr , '## Listing structure files'
	structure_files = options.structure
	structure_files_list = structure_files.split(",")

	if options.gff3 :
		print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Loading GFF3 annotation'
		print >> sys.stderr , '## Loading GFF3 annotation'
		annotation_gff3 , mRNA_db = read_gff3(options.gff3)
	else :
		annotation_gff3 = ""

	if options.bed :
		print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Loading BED file'
		print >> sys.stderr , '## Loading BED file'
		bed_regions = read_bed_sorted_list( options.bed )
		# bed_regions[line] = [chrom , chromStart(0->) , chromEnd(1->) , name , score , strand , ... ]
	else :
		bed_regions = ""

	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Loading structure(s)'
	print >> sys.stderr , '## Loading structure(s)'
	structure_file_format = str(options.format).lower()
	print >> sys.stderr , '### Structure file(s) format: ' + structure_file_format

	if options.skipoverlap :
		# Convert input structures to AGP and export results
		if structure_file_format == "bed" :
			print >> sys.stdout , "[MEMO] Structure in BED format. A new sequence will be generated for each input file. Progressive numeric id will follow input order."
			id = 0
			agp_db = {}
			for file_name in structure_files_list :
				id += 1
				bed_db = read_bed_sorted_list(file_name)
				print >> sys.stderr , '### Input file id: ' + str(id) + " | Input BED File name: " + file_name + " | Corresponfing sequence IDs: " + options.prefix + "_" + str(id)
				if agp_db == {} :
					agp_db = bed_to_agp_onefile(bed_db , gap_size , options.prefix + "_" + str(id))
				else :
					agp_db.update(bed_to_agp_onefile(bed_db , gap_size , options.prefix + "_" + str(id)))

		elif structure_file_format == "agp" :
			#### Read AGP
			agp_db = {}
			for file_name in structure_files_list :
				if agp_db == {} :
					agp_db = read_agp(file_name)
				else :
					agp_db.update(read_agp(file_name))

			mode = "old_to_new"
			if options.mode :
				print >> sys.stdout , "[WARNING] Reverse use of structure AGP file requested [--reverse]"
				print >> sys.stdout , "[WARNING] FASTA file will not be produced"
				no_fasta = True
				mode = "new_to_old"
				# Invert AGP
				new_agp_db = invert_agp(agp_db)
				agp_db = clean_self_agp(new_agp_db)
				#agp_file_name = options.out + ".inverted.agp"
				#agp_file_name = write_agp( new_agp_db , agp_file_name )

		elif structure_file_format == "block" :
			agp_db = {}
			for file_name in structure_files_list :
				block_db = read_block(file_name)
				if agp_db == {} :
					agp_db = block_to_agp(block_db , options.gap_size )
				else :
					agp_db.update( block_to_agp(block_db) )

		else :
			print >> sys.stderr , "[ERROR] Structure file format " + str(options.format) + " unknown."
			print >> sys.stdout , "[ERROR] Structure file format " + str(options.format) + " unknown."
			parser.print_help()
			sys.exit(1)

		if options.add_unplaced :
			print >> sys.stderr , '### Adding unplaced sequences to output'
			used_sequences = []
			# agp_db[seq_id][int(start)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation ]
			for seq_id in agp_db.keys() :
				for start in agp_db[seq_id].keys() :
					if agp_db[seq_id][int(start)][4] == "W" :
						used_sequences.append(agp_db[seq_id][int(start)][5])

			for seq_id in sorted(fasta_len_dict.keys()) :
				if seq_id not in used_sequences :
					seq_length = fasta_len_dict(seq_id)
					agp_db[seq_id][1] = [ seq_id , 1 , seq_length , 1 , "W" , seq_id , 1 , seq_length ,  "+" ]

	else :
		# Convert structures to BLOCK,
		# identify overlaps between blocks
		# correct coordinates
		# convert to AGP
		# export results
		if structure_file_format == "block" :
			block_db = {}
			for file_name in structure_files_list :
				print >> sys.stderr , '### Loading ' + file_name
				if block_db == {} :
					block_db = read_block(file_name)
				else :
					block_db.update( read_block(file_name) )

		elif structure_file_format == "agp" :
			if options.mode :
				print >> sys.stdout , "[ERROR] Reverse use of AGP structure file incompatible with overlap search and correction."
				print >> sys.stdout , "[ERROR] Reverse use of AGP structure file incompatible with overlap search and correction."
				sys.exit(1)
			agp_db = {}
			for file_name in structure_files_list :
				print >> sys.stderr , '### Loading ' + file_name
				if agp_db == {} :
					agp_db = read_agp(file_name)
				else :
					agp_db.update(read_agp(file_name))
			block_db = agp_to_block( agp_db )

		elif structure_file_format == "bed" :
			id = 0
			agp_db = {}
			for file_name in structure_files_list :
				print >> sys.stderr , '### Loading ' + file_name
				id += 1
				bed_db = read_bed_sorted_list(file_name)
				print >> sys.stderr , '### Input file id: ' + str(id) + " | Input BED File name: " + file_name + " | Corresponfing sequence IDs: " + options.prefix + "_" + str(id)
				if agp_db == {} :
					agp_db = bed_to_agp_onefile(bed_db , gap_size , options.prefix + "_" + str(id))
				else :
					agp_db.update(bed_to_agp_onefile(bed_db , gap_size , options.prefix + "_" + str(id)))
			block_db = agp_to_block( agp_db )

		else :
			print >> sys.stderr , "[ERROR] Structure file format " + str(options.format) + " unknown."
			print >> sys.stdout , "[ERROR] Structure file format " + str(options.format) + " unknown."
			parser.print_help()
			sys.exit(1)

		# Correct agp regions with smart overlap dodging
		temp_dir = options.out + ".temp_dir"
		mkdir(temp_dir)
		print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] = Refining structure(s) based on overlap'
		print >> sys.stderr , '# Refining structure(s) based on overlap'
		agp_db , harmed_loci = dodge_overlaps( block_db , fasta_dict , fasta_len_dict , int(options.gap_size) , int(options.spacer) , int(options.cores) , options.mapper , annotation_gff3 , temp_dir , paths, options.add_unplaced)
		if not harmed_loci == [] :
			harmed_loci_file = open( options.out + ".loci_to_check.txt" , "w+")
			for name in harmed_loci :
				print >> harmed_loci_file , name
			harmed_loci_file.close()

		if options.ignoreids :
			agp_db = rename_agp_sequences(agp_db , options.prefix )


	# Convert legacy agp if given
	if options.agp :
		print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] = Translating legacy components position from legacy agp file'
		print >> sys.stderr , '# Translating legacy components position from legacy agp file'
		old_agp = read_agp(options.agp)
		new_to_legacy_agp_db = agp_translate_agp(agp_db , old_agp)

	# Convert BED file if given
	if options.bed :
		print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] = Translating coordinates of features in bed the file'
		print >> sys.stderr , '# Translating coordinates of features in bed the file'
		new_bed = translate_bed_sorted_list( bed_regions ,  agp_db )

	# Writing output files
	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] = Writing output files'
	print >> sys.stderr , '# Writing output files'

	if options.skipoverlap :
		if not options.noagp :
			agp_file_name = options.out + ".structure.agp"
			agp_file_name = write_agp( agp_db , agp_file_name )
		export_from_agp(options.out, no_fasta, agp_db, fasta_dict, "old_to_new" , fasta_len_dict, annotation_gff3)

	else :
		export_from_agp(options.out, no_fasta, agp_db, fasta_dict, "old_to_new", fasta_len_dict, annotation_gff3)
		agp_file_name = options.out + ".structure.agp"
		agp_file_name = write_agp( agp_db , agp_file_name )

	if options.agp :
		legacy_agp_file_name = options.out + ".legacy_structure.agp"
		legacy_agp_file_name = write_agp( new_to_legacy_agp_db , legacy_agp_file_name )

	if options.bed :
		out_bed_file_name = options.out + ".bed"
		out_bed_file = open(out_bed_file_name , 'w')
		for line in sorted(new_bed) :
			print >> out_bed_file, "\t".join([str(x) for x in line])
		out_bed_file.close()

	###### END ######

	print >> sys.stdout , "------------------------------"
	print >> sys.stdout , "- Done"
	print >> sys.stdout , "------------------------------"
	print >> sys.stderr , "##############################"
	print >> sys.stderr , "# Done"
	print >> sys.stderr , "##############################"


if __name__ == '__main__':
	main()
