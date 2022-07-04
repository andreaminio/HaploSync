#!/usr/bin/env python

import argparse
import gzip
from lib_files.HaploFunct import *
from lib_files.AGP_lib import *
from lib_files.GFF_lib import *
from lib_files.FASTA_lib import *
from lib_files.map_lib import *



def main() :

	###### Options and help ######
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fasta", default=False, dest="fasta",
						help="FASTA file(s) with genomic sequences", metavar="genome.fasta [Required]")
	parser.add_argument("-g", "--gff", default=False, dest="gff",
						help="Annotation file(s) in GFF3 format", metavar="annotation.gff3 [Required]")
	parser.add_argument("-a", "--annotation", default=False, dest="annotation",
						help="Text file with the functional annotation of each transcript", metavar="functional_annotation.txt")
	parser.add_argument("-c", "--correspondence", dest="corr", default=False ,
						help="Tab text file(s) of corresponding sequence names between haplotypes for each chromosome and, if present, the reference genome. Tab separated file with 3/4 columns: 1) Chromosome ID, 2) Haplotype 1 ID, 3) Haplotype 2 ID, 4) [optional] Reference ID [Required]", metavar="Hap1_to_Hap2.txt")

	# Add marker BED input
	parser.add_argument( "-b" , "--markers_bed" , dest="markers_hits" ,
						help="Marker position on input sequences", metavar="markers.bed")
	parser.add_argument("--markers_map", dest="marker_map",
					help="Table of markers order when sorted genomic position. Tab separated file with 3 columns: 1) Chromosome ID, 2) Marker sorting position, 3) Marker ID.", metavar="markers_list")

	# Add pseudomolecule structure AGP
	parser.add_argument("--agp" , dest="agp" ,
						help="AGP structure of input sequences", metavar="structure.agp")
	parser.add_argument("--input_groups" , dest="input_groups",
					help="Tab separated values file of group association of input sequences, it will be used for structure QC if '--legacy_groups' is not given. Two columns reporting the sequences ID and the association groups id. Sequence can be associated to multiple groups with multiple rows."  , metavar="input_groups.tsv")

	# Add pseudomolecule legacy structure AGP
	parser.add_argument("--legacy_agp" , dest="legacy_agp" ,
						help="AGP structure of input sequences based on older component sequences", metavar="legacy_structure.agp")
	parser.add_argument("--legacy_groups" , dest="legacy_groups",
					help="Tab separated values file of group association of components present the input sequences. Two columns reporting the component ID and the association group it belongs to. Sequence can be associated to multiple groups with multiple rows.", metavar="legacy_groups.tsv")


	# QC for rejected
	parser.add_argument("--rejected_list" , dest="rejected" ,
						help="List of rejected sequences for each chromosome, run QC of rejection", metavar="rejected.list")

	parser.add_argument("-o", "--out", dest="out", default="out",
						help="Output files prefix [default: out]", metavar="NAME")

	parser.add_argument("--hit_identity", dest="hit_identity", default="90",
						help="Mimimum genome mapping hit identity [default: 90]", metavar="N")
	parser.add_argument("--hit_coverage", dest="hit_length", default="3000",
						help="Mimimum genome mapping hit length [default: 3000]", metavar="N")

	parser.add_argument("--gene_identity", dest="gene_identity", default="95",
						help="Mimimum gene mapping identity [default: 95]", metavar="N")
	parser.add_argument("--gene_coverage", dest="gene_coverage", default="95",
						help="Mimimum gene mapping coverage [default: 95]", metavar="N")

	parser.add_argument("--unbalanced_ratio", dest="unbalanced_ratio" , default="0.33",
						help="Gene count ratio between haplotype to call the locus underrepresented [values range: 0-1]", metavar="N")

	parser.add_argument("-r", "--reference", default=False, dest="reference",
						help="FASTA file(s) of reference genomes sequences (haploid)", metavar="reference.fasta")

	parser.add_argument("--reuse_mappings" , dest="reuse_mappings", default=False, action="store_true",
						help="If set, alignments present in the output folder are reused and not overwritten by performing alignments again [Default: overwrite]")
	parser.add_argument("--reuse_dotplots" , dest="reuse_dotplots", default=False, action="store_true",
						help="If set, dotplots present in the output folder are reused and not overwritten [Default: overwrite]")
	parser.add_argument("--reuse_gmap" , dest="reuse_gmap", default=False, action="store_true",
						help="If set, CDS mapping with GMAP are reused and not overwritten by performing again the analysis [Default: overwrite]")
	parser.add_argument("--skip_dotplots_by_chr" , dest="skip_dotplots", default=False, action="store_true",
						help="If set, prevents the production of dotplots comparing each chromosome sequence to any other chromosome sequence. Whole genome dotplot is produced anyway. [Default: overwrite]")

	# TODO: allow to use a custom set of CDS sequences instead of annotations
	parser.add_argument("-c", "--cds", default=False, dest="cds",
						#help="CDS sequences to use to generate a temporary annotation", metavar="cds.fasta [Required]")
						help=argparse.SUPPRESS )

	# Mapper selection is hidden as only gmap alignment is supported
	parser.add_argument("-m" , "--mapper", dest="mapper", default="gmap",
						#help="Mapping tool to use [Default: gmap]" , metavar="[gmap|blat]" )
						help=argparse.SUPPRESS )
	parser.add_argument("-t", "--threads", dest="cores", default=4,
						help="Cores used in mapping process [default: 4]", metavar="N")
	parser.add_argument("--feature", dest="feature", default="CDS",
						help="If GFF is used, feature type to use for mapping. Choice of CDS or mRNA [default: CDS]", metavar="[CDS|mRNA]")
	parser.add_argument("-w", "--window", dest="window", default=10,
						help="Window size (number of genes) for search of blocks of genes with unbalanced count between alleles. 0 disables search [default: 10]", metavar="N")
	parser.add_argument("--allowed", dest="allowed", default=5,
						help="Allowed number of unbalance gene per window. [default: 5]", metavar="N")



	scriptDirectory = os.path.dirname(os.path.realpath(__file__)) + "/support_scripts"
	print >> sys.stdout, "Running HaploDup tool from HaploSync version " + get_version()
	print >> sys.stdout, "To reproduce this run use the following command: " + " ".join( pipes.quote(x) for x in sys.argv)
	print >> sys.stdout, "----"
	# Sanity Check

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)

	options = parser.parse_args()

	if not options.fasta :
		print >> sys.stderr , "[ERROR] Genome FASTA file(s) missing"
		parser.print_help()
		sys.exit(1)

	if not ( options.gff or options.cds ) :
		print >> sys.stderr , "[ERROR] Annotation file missing"
		parser.print_help()
		sys.exit(1)

	if not options.corr :
		print >> sys.stderr , "[ERROR] Annotation file or and CDS sequences missing"
		parser.print_help()
		sys.exit(1)

	haplodup_dir = options.out + ".HaploDup_dir"

	if (options.reuse_mappings or options.reuse_dotplots or options.reuse_gmap) and (not os.path.exists(haplodup_dir)):
		print >> sys.stderr , "[ERROR] Required to reuse output folder content but output folder do not exist (prefix used: " + options.out + ")"
		parser.print_help()
		sys.exit(1)

	paths = set_paths(os.path.join(sys.path[0], 'HaploSync.conf.toml'))

	iden_threshold = float(options.gene_identity)
	cov_threshold = float(options.gene_coverage)

	hit_iden = options.hit_identity
	hit_len = options.hit_length
	
	unbalanced_ratio = options.unbalanced_ratio 

	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] = Read inputs'
	print >> sys.stderr , '# Read inputs'
	# Make directory

	coord_tables = {}
	plot_files = {}
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Performing HaploDup"
	mkdir( "./" + haplodup_dir )

	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Loading FASTA sequences'
	print >> sys.stderr , '## Loading FASTA sequences'
	fasta_files = options.fasta
	fasta_files_list = fasta_files.split(",")
	fasta_dict = {}
	for file_name in fasta_files_list :
		print >> sys.stderr , '### Loading FASTA file: ' + file_name
		if fasta_dict == {} :
			fasta_dict = read_fasta(file_name)
		else :
			fasta_dict.update(read_fasta(file_name))
	fasta_len_dict = get_length_from_fasta_db(fasta_dict)
	complete_fasta_file_name = write_fasta_from_db( fasta_dict , haplodup_dir + "/" + options.out + ".input.fasta" , False )

	if options.reference :
		reference = read_fasta(options.reference)
		reference_len = get_length_from_fasta_db(reference)
		reference_file_name = write_fasta_from_db( reference , haplodup_dir + "/" + options.out + ".ref.fasta" , False)

	print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Loading sequence correspondence information'
	print >> sys.stderr , '## Loading sequence correspondence information'
	corr_files = options.corr
	corr_files_list = corr_files.split(",")
	pairs = []

	hap1_fasta = {}
	hap2_fasta = {}

	hap1_ids = []
	hap2_ids = []
	ref_ids = []
	ref_to_hap1 = {}
	hap1_to_ref = {}
	ref_to_hap2 = {}
	hap2_to_ref = {}
	hap2_to_hap1 = {}
	hap1_to_hap2 = {}

	chr_ids = []
	chr_to_hap1 = {}
	chr_to_hap2 = {}
	chr_to_ref = {}
	hap1_to_chr = {}
	hap2_to_chr = {}
	ref_to_chr = {}

	for file_name in corr_files_list :
		print >> sys.stderr , '### Loading correspondence file: ' + file_name
		if not options.reference :
			for line in open(file_name) :
				try :
					chr, seq1 , seq2 = line.rstrip().split("\t")
				except :
					print >> sys.stdout , '[ERROR] Quitting'
					print >> sys.stderr , '[ERROR] Correspondence file  ' + file_name + ' has an unexpected number of columns (expected 3, in order: Chr_id, Hap1_id and Hap2_id)'
					sys.exit(1)
				else :
					if seq1 not in fasta_dict or seq2 not in fasta_dict :
						print >> sys.stdout , '[ERROR] Quitting'
						print >> sys.stderr , '[ERROR] Correspondence file  ' + file_name + " contains sequences missing from the FASTA file:" + seq1
						sys.exit(1)
					if seq2 not in fasta_dict :
						print >> sys.stdout , '[ERROR] Quitting'
						print >> sys.stderr , '[ERROR] Correspondence file  ' + file_name + " contains sequences missing from the FASTA file:" + seq2
						sys.exit(1)
					hap1_ids.append(seq1)
					hap1_fasta[seq1] = fasta_dict[seq1]
					hap2_ids.append(seq2)
					hap2_fasta[seq2] = fasta_dict[seq2]
					pairs.append([seq1 , seq2])
					hap1_to_hap2[seq1] = seq2
					hap2_to_hap1[seq2] = seq1
					chr_ids.append(chr)
					chr_to_hap1[chr] = seq1
					chr_to_hap2[chr] = seq2
					hap1_to_chr[seq1] = chr
					hap2_to_chr[seq2] = chr
		else :
			for line in open(file_name) :
				try :
					chr , seq1 , seq2 , ref = line.rstrip().split("\t")
				except :
					print >> sys.stdout , '[ERROR] Quitting'
					print >> sys.stderr , '[ERROR] Correspondence file  ' + file_name + ' has an unexpected number of columns (expected 4, in order: Chr_id, Hap1_id, Hap2_id, Ref_id)'
					sys.exit(1)
				else :
					if seq1 not in fasta_dict or seq2 not in fasta_dict :
						print >> sys.stdout , '[ERROR] Quitting'
						print >> sys.stderr , '[ERROR] Correspondence file  ' + file_name + " contains sequences missing from the FASTA file:" + seq1
						sys.exit(1)
					if seq2 not in fasta_dict :
						print >> sys.stdout , '[ERROR] Quitting'
						print >> sys.stderr , '[ERROR] Correspondence file  ' + file_name + " contains sequences missing from the FASTA file:" + seq2
						sys.exit(1)
					hap1_ids.append(seq1)
					hap1_fasta[seq1] = fasta_dict[seq1]
					hap2_ids.append(seq2)
					hap2_fasta[seq2] = fasta_dict[seq2]
					ref_ids.append(ref)
					pairs.append([seq1 , seq2])
					hap1_to_hap2[seq1] = seq2
					hap2_to_hap1[seq2] = seq1
					hap1_to_ref[seq1] = ref
					hap2_to_ref[seq2] = ref
					ref_to_hap2[ref] = seq2
					ref_to_hap1[ref] = seq1
					chr_ids.append(chr)
					chr_to_hap1[chr] = seq1
					chr_to_hap2[chr] = seq2
					chr_to_ref[chr] = ref
					hap1_to_chr[seq1] = chr
					hap2_to_chr[seq2] = chr
					ref_to_chr[ref] = chr

	# Generate separate databases for the 2 haplotypes
	fasta_db_1 = {}
	fasta_db_2 = {}
	fasta_1_len = {}
	fasta_2_len = {}
	for hap1_seq_id in sorted(hap1_to_hap2.keys()) :
		hap2_seq_id = hap1_to_hap2[hap1_seq_id]
		fasta_db_1[hap1_seq_id] = fasta_dict[hap1_seq_id]
		fasta_1_len[hap1_seq_id] = len(fasta_dict[hap1_seq_id])
		fasta_db_2[hap2_seq_id] = fasta_dict[hap2_seq_id]
		fasta_2_len[hap2_seq_id] = len(fasta_dict[hap2_seq_id])

	# Read GFF files
	# For now the the annotation in GFF is mandatory
	# TODO: Use CDS sequences to build a temporary reference annotation

	annotation_db = {}

	if options.gff :
		if options.feature.upper() == "CDS" :
			feat_to_extract = "CDS"
		elif options.feature.upper() == "MRNA" :
			feat_to_extract = "mRNA"
		else :
			print >> sys.stdout , '[ERROR] Unknown mapping feature ('+options.feature+') requested'
			print >> sys.stderr , '[ERROR] Unknown mapping feature ('+options.feature+') requested'
			sys.exit(1)

		print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Loading GFF3 annotation'
		print >> sys.stderr , '## Loading GFF3 annotation'
		gff_files = options.gff
		gff_files_list = gff_files.split(",")
		gff_db = {}
		mRNA_db = {}
		for file_name in gff_files_list :
			print >> sys.stderr , '### Loading GFF3 file: ' + file_name
			if gff_db == {} :
				gff_db , mRNA_db = read_gff3(file_name)
			else :
				gff_tmp , mrna_tmp = read_gff3(file_name)
				gff_db.update(gff_tmp)
				mRNA_db.update(mrna_tmp)
		#else :
		#	# options.cds has been set instead
		#	# TODO: make one loci annotation if CDS sequences are given
		mRNA_to_gene_db = get_gene2mRNA_from_db(gff_db)
		# mRNA_to_gene_db[mRNA_id] = gene_id

		if options.annotation :
			# Read the functional annotation associated to the GFF file
			print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Loading functional annotation'
			print >> sys.stderr , '## Loading functional annotation'
			annotation_files = options.annotation
			annotation_files_list = annotation_files.split(",")
			for file_name in annotation_files_list :
				print >> sys.stderr , '### Loading annotation file: ' + file_name
				for line in open(file_name) :
					if line == "" or line[0] == "#" :
						continue
					else :
						feat_name , description = line.rstrip().split("\t")[0:2]
						description = description.upper()
						if feat_name in mRNA_to_gene_db :
							# Annotation associated to mRNAs -> add to gene annotation (uniquely)
							gene_id =  mRNA_to_gene_db[feat_name]
							if not gene_id in annotation_db :
								annotation_db[gene_id] = {}
							if description not in annotation_db[gene_id] :
								annotation_db[gene_id][description] = []
						else :
							if not feat_name in annotation_db :
								annotation_db[feat_name] = {}
							if description not in annotation_db[feat_name] :
								annotation_db[feat_name][description] = []

			annotation_db = add_counts_to_dict(annotation_db)

	# Read marker BED input
	if options.markers_hits :
		print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] == Loading markers position'
		print >> sys.stderr , '## Loading markers position'
		# Copy the file in haplodup_dir
		all_markers_file_name = "all_markers.bed"
		dup_markers_file_name = "duplicated_markers.bed"
		all_markers_file_name_fullpath = haplodup_dir + "/all_markers.bed"
		dup_markers_file_name_fullpath = haplodup_dir + "/duplicated_markers.bed"
		seq_all_markers_file = open(all_markers_file_name_fullpath , 'w')
		seq_duplicated_markers_file = open(dup_markers_file_name_fullpath , 'w')
		## Write markers files
		markers_db = {}
		for line in open(options.markers_hits) :
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
		dup_markers_file_name = ""


	# Pseudomolecule structure AGP
	if options.agp :
		agp_db = read_agp(options.agp)
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
	else :
		structure_file = ""
		structure_file_fullpath = ""

	# Pseudomolecule legacy structure AGP
	if options.legacy_agp :
		legacy_agp = read_agp(options.legacy_agp)
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
		legacy_structure_file_fullpath = ""

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
	hap2_ids = ",".join(sorted(fasta_db_2.keys()))

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Mapping sequences"
	query_1_file = haplodup_dir + "/" + options.out + ".1.fasta"
	write_fasta_from_db( fasta_db_1 , query_1_file )
	query_2_file = haplodup_dir + "/" + options.out + ".2.fasta"
	write_fasta_from_db( fasta_db_2 , query_2_file )

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap1 vs Hap1"
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Mapping"
	if options.reuse_mappings :
		outfile_prefix = "Hap1.on.Hap1"
	else :
		outfile_prefix = map_nucmer_dotplot("Hap1" , query_1_file , "Hap1" , query_1_file , haplodup_dir , options.cores , paths , False )
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Converting files"
	coord_tables["Hap1_vs_Hap1"] = make_coords_table( outfile_prefix , haplodup_dir , command_line)

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap2 vs Hap2"
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Mapping"
	if options.reuse_mappings :
		outfile_prefix = "Hap2.on.Hap2"
	else:
		outfile_prefix = map_nucmer_dotplot("Hap2" , query_2_file , "Hap2" , query_2_file , haplodup_dir , options.cores , paths , False  )
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Converting files"
	coord_tables["Hap2_vs_Hap2"] = make_coords_table( outfile_prefix , haplodup_dir , command_line)

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap2 vs Hap1"
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Mapping"
	if options.reuse_mappings :
		outfile_prefix = "Hap2.on.Hap1"
	else :
		outfile_prefix = map_nucmer_dotplot("Hap1" , query_1_file , "Hap2" , query_2_file , haplodup_dir , options.cores , paths , False  )
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Converting files"
	coord_tables["Hap2_vs_Hap1"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , coord_tables["Hap1_vs_Hap1"] ]

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap1 vs Hap2"
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Mapping"
	if options.reuse_mappings :
		outfile_prefix = "Hap1.on.Hap2"
	else :
		outfile_prefix = map_nucmer_dotplot("Hap2" , query_2_file , "Hap1" , query_1_file , haplodup_dir , options.cores , paths , False  )
	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Converting files"
	coord_tables["Hap1_vs_Hap2"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , coord_tables["Hap2_vs_Hap2"] ]


	if options.reference :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap1 vs Reference"
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Mapping"
		if options.reuse_mappings :
			outfile_prefix = "Hap1.on.Ref"
		else :
			outfile_prefix = map_nucmer_dotplot("Ref" , options.reference , "Hap1" , query_1_file , haplodup_dir , options.cores , paths , False )
			# outfile_prefix.delta and outfile_prefix.coords (show-cords -c) are generated
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Converting files"
		coord_tables["Hap1_vs_Reference"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , "" ]

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Reference vs Hap1 "
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Mapping"
		if options.reuse_mappings :
			outfile_prefix = "Ref.on.Hap1"
		else:
			outfile_prefix = map_nucmer_dotplot( "Hap1" , query_1_file , "Ref" , options.reference , haplodup_dir , options.cores , paths , False )
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Converting files"
		coord_tables["Reference_vs_Hap1"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , coord_tables["Hap1_vs_Hap1"] ]

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap2 vs Reference"
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Mapping"
		if options.reuse_mappings :
			outfile_prefix = "Hap2.on.Ref"
		else :
			outfile_prefix = map_nucmer_dotplot("Ref" , options.reference , "Hap2" , query_2_file , haplodup_dir , options.cores , paths , False )
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Converting files"
		coord_tables["Hap2_vs_Reference"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , "" ]

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Reference vs Hap2"
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Mapping"
		if options.reuse_mappings :
			outfile_prefix = "Ref.on.Hap2"
		else :
			outfile_prefix = map_nucmer_dotplot("Hap2" , query_2_file , "Ref" , options.reference , haplodup_dir , options.cores , paths , False )
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Converting files"
		coord_tables["Reference_vs_Hap2"] = [make_coords_table( outfile_prefix , haplodup_dir , command_line) , coord_tables["Hap2_vs_Hap2"] ]


	if options.gff :
		# Map genes, generate copy number counts, render by chromosome reports
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Generating files about gene map count for fusion dedup"

		## Generate CDS sequences from first haplotype
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Generating CDS sequences"
		if options.reuse_gmap :
			gmap_results = haplodup_dir + "/CDS.on.ref.gmap.gff3"

		else :
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Mapping CDSs on genome assembly"
			CDS_file = get_sequence( gff_db , fasta_dict , haplodup_dir + "/new" , "CDS")
			index_dir = haplodup_dir + "/gmap_index"

			# The new gmap

			mkdir(index_dir)

			# Hap1
			## index
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Hap1 "
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Indexing genome "
			indexing_out_file = open( haplodup_dir + "/gmap_index.1.log" ,"w" )
			indexing_err_file = open( haplodup_dir + "/gmap_index.1.err" ,"w" )
			hap1_name = haplodup_dir + "/new.hap1.fasta"
			hap1_file = write_fasta_from_db( hap1_fasta , hap1_name)
			indexing_command = "gmap_build -D " + index_dir + " -d hap1.fasta " + hap1_name
			indexProcess = subprocess.Popen(indexing_command, shell=True, stdout=indexing_out_file , stderr=indexing_err_file)
			output, error = indexProcess.communicate()
			indexing_out_file.close()
			indexing_err_file.close()

			# Gmap CDS on results
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Mapping with Gmap "
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
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Hap2 "
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Indexing genome "
			indexing_out_file = open( haplodup_dir + "/gmap_index.2.log" ,"w" )
			indexing_err_file = open( haplodup_dir + "/gmap_index.2.err" ,"w" )
			hap2_name = haplodup_dir + "/new.hap2.fasta"
			hap2_file = write_fasta_from_db( hap2_fasta , hap2_name)
			indexing_command = "gmap_build -D " + index_dir + " -d hap2.fasta " + hap2_name
			indexProcess = subprocess.Popen(indexing_command, shell=True, stdout=indexing_out_file , stderr=indexing_err_file)
			output, error = indexProcess.communicate()
			indexing_out_file.close()
			indexing_err_file.close()
			# Gmap CDS on results
			print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] ==== Mapping with Gmap "
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
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Extracting valid alignments (identity > " + str(iden_threshold) + "% , coverage> " + str(cov_threshold)  + "%)"
		gmap_hits_hap1 , gmap_hits_hap2 = read_gmap_results_Hap(gmap_results, cov_threshold , iden_threshold , "by_locus", mRNA_to_gene_db)
		test_1 = open( haplodup_dir + "/gmap_hits_hap1.txt" ,'w')
		for chr in sorted(gmap_hits_hap1.keys()) :
			for locus in sorted(gmap_hits_hap1[chr].keys()) :
				for hit in sorted(gmap_hits_hap1[chr][locus]) :
					print >> test_1 , chr + "\t" + "\t".join([ str(x) for x in hit ]) + "\t" + locus
		test_1.close()
		## Make Hap1 gene table
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Extracting loci positions on Hap1"
		hap1_genes = gff3_filter2table_Hap(gff_db, "gene", "_Hap1_")
		#print >> sys.stderr, "### hap1 gene"
		#print >> sys.stderr, hap1_genes

		test_2 = open( haplodup_dir + "/gmap_hits_hap2.txt" ,'w')
		for chr in sorted(gmap_hits_hap2.keys()) :
			for locus in sorted(gmap_hits_hap2[chr].keys()) :
				for hit in sorted(gmap_hits_hap2[chr][locus]) :
					print >> test_2 , chr + "\t" + "\t".join([ str(x) for x in hit ]) + "\t" + locus
		test_2.close()
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] === Extracting loci positions on Hap2"
		hap2_genes = gff3_filter2table_Hap(gff_db, "gene", "_Hap2_")
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



	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Generating dotplots"

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap1 vs Hap1"
	outfile_prefix = "Hap1.on.Hap1"
	plot_files["Hap1_vs_Hap1"] = {}
	plot_files["Hap1_vs_Hap1"]["Whole"] , plot_files["Hap1_vs_Hap1"]["All_Dotplots"] = whole_genome_dotplot( hap1_ids , hap1_ids , outfile_prefix, haplodup_dir, "Hap1_vs_Hap1", coord_tables["Hap1_vs_Hap1"] , options.reuse_dotplots , options.skip_dotplots)

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap2 vs Hap2"
	outfile_prefix = "Hap2.on.Hap2"
	plot_files["Hap2_vs_Hap2"] = {}
	plot_files["Hap2_vs_Hap2"]["Whole"] , plot_files["Hap2_vs_Hap2"]["All_Dotplots"] = whole_genome_dotplot( hap2_ids, hap2_ids , outfile_prefix, haplodup_dir, "Hap2_vs_Hap2", coord_tables["Hap2_vs_Hap2"] , options.reuse_dotplots , options.skip_dotplots)

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap2 vs Hap1"
	outfile_prefix = "Hap2.on.Hap1"
	plot_files["Hap2_vs_Hap1"] = {}
	plot_files["Hap2_vs_Hap1"]["Whole"] , plot_files["Hap2_vs_Hap1"]["All_Dotplots"] = whole_genome_dotplot( hap1_ids , hap2_ids , outfile_prefix, haplodup_dir, "Hap2_vs_Hap1", coord_tables["Hap2_vs_Hap1"][0] , options.reuse_dotplots , options.skip_dotplots)

	print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap1 vs Hap2"
	outfile_prefix = "Hap1.on.Hap2"
	plot_files["Hap1_vs_Hap2"] = {}
	plot_files["Hap1_vs_Hap2"]["Whole"], plot_files["Hap1_vs_Hap2"]["All_Dotplots"] = whole_genome_dotplot( hap2_ids , hap1_ids, outfile_prefix, haplodup_dir, "Hap1_vs_Hap2", coord_tables["Hap1_vs_Hap2"][0] , options.reuse_dotplots , options.skip_dotplots)

	if options.reference :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap1 vs Reference"
		outfile_prefix = "Hap1.on.Ref"
		plot_files["Hap1_vs_Reference"] = {}
		plot_files["Hap1_vs_Reference"]["Whole"] , plot_files["Hap1_vs_Reference"]["All_Dotplots"] = whole_genome_dotplot( reference_ids , hap1_ids , outfile_prefix, haplodup_dir, "Hap1_vs_Reference", coord_tables["Hap1_vs_Reference"][0] , options.reuse_dotplots , options.skip_dotplots)

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Reference vs Hap1 "
		outfile_prefix = "Ref.on.Hap1"
		plot_files["Reference_vs_Hap1"] = {}
		plot_files["Reference_vs_Hap1"]["Whole"] , plot_files["Reference_vs_Hap1"]["All_Dotplots"] = whole_genome_dotplot( hap1_ids , reference_ids , outfile_prefix, haplodup_dir, "Reference_vs_Hap1", coord_tables["Reference_vs_Hap1"][0], options.reuse_dotplots , options.skip_dotplots)

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Hap2 vs Reference"
		outfile_prefix = "Hap2.on.Ref"
		plot_files["Hap2_vs_Reference"] = {}
		plot_files["Hap2_vs_Reference"]["Whole"] , plot_files["Hap2_vs_Reference"]["All_Dotplots"] = whole_genome_dotplot( reference_ids , hap2_ids , outfile_prefix, haplodup_dir, "Hap2_vs_Reference", coord_tables["Hap2_vs_Reference"][0] , options.reuse_dotplots , options.skip_dotplots)

		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] == Reference vs Hap2"
		outfile_prefix = "Ref.on.Hap2"
		plot_files["Reference_vs_Hap2"] = {}
		plot_files["Reference_vs_Hap2"]["Whole"] , plot_files["Reference_vs_Hap2"]["All_Dotplots"] = whole_genome_dotplot( hap2_ids , reference_ids , outfile_prefix, haplodup_dir, "Reference_vs_Hap2", coord_tables["Reference_vs_Hap2"][0] , options.reuse_dotplots , options.skip_dotplots)

	if options.gff :
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
					plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, "" , "" , "" ,  "" , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio )
					#																			coords, 		coords_self, 			workdir, 				output_dir, 				                                     queryID, refID, structure , legacy , markers, dup_markers , hap1ID , hap2ID , hap1Len , hap2Len , counts_hap1 ="diploid_gene_count_trace.hap1.txt", counts_hap2 ="diploid_gene_count_trace.hap2.txt", min_align ="3000", similarity ="90", ratio="0.33") :
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
			elif comparison == "Reference_vs_Hap1" :
				for queryID in sorted(ref_ids) :
					refID = ref_to_hap1[queryID]
					hap1ID = refID
					hap2ID = hap1_to_hap2[hap1ID]
					hap1Len = fasta_1_len[hap1ID]
					hap2Len = fasta_2_len[hap2ID]
					plot_files[comparison]["Reports"][queryID] = {}
					plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio )
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_pair_pdf_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file   , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , "diploid_gene_count_trace.hap1.txt", hit_len , hit_iden, unbalanced_ratio )
					#													make_pair_pdf_report(coords						  , coords_self                       , workdir						  , output_dir				     , queryID, refID, structure        , legacy                , markers               , dup_markers           , hap1ID , hap2ID , counts_hap1                        , min_align , similarity , ratio ) :																																												                                    refID, structure = "" , legacy = "" , markers = "" , dup_markers = "" , hap1ID
			elif comparison == "Hap1_vs_Hap1" :
				for queryID in sorted(hap1_ids.split(",")) :
					refID = queryID
					hap1ID = queryID
					hap2ID = hap1_to_hap2[hap1ID]
					hap1Len = fasta_1_len[hap1ID]
					hap2Len = fasta_2_len[hap2ID]
					plot_files[comparison]["Reports"][queryID] = {}
					plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name ,  dup_markers_file_name , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio )
					#																			coords, 					coords_self, 						workdir, 							output_dir, 				queryID, refID, structure 	   , legacy                , markers,                dup_markers           ,  hap1ID , hap2ID , hap1Len , hap2Len  , counts_hap1                         , counts_hap2                      , min_align , similarity , ratio) :
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_pair_pdf_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , "diploid_gene_count_trace.hap1.txt", hit_len , hit_iden, unbalanced_ratio )
					#													make_pair_pdf_report(coords						  , coords_self                       , workdir						  , output_dir				     , queryID, refID, structure      , legacy                , markers               , dup_markers           , hap1ID , hap2ID , counts_hap1                        , min_align , similarity , ratio ) :
			elif comparison == "Hap2_vs_Reference" :
				for queryID in sorted(hap2_ids.split(",")) :
					refID = hap2_to_ref[queryID]
					hap1ID = hap2_to_hap1[queryID]
					hap2ID = queryID
					hap1Len = fasta_1_len[hap1ID]
					hap2Len = fasta_2_len[hap2ID]
					plot_files[comparison]["Reports"][queryID] = {}
					plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, "" , "" , "" ,  "" , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio )
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
			elif comparison == "Reference_vs_Hap2" :
				for queryID in sorted(ref_ids) :
					refID = ref_to_hap2[queryID]
					hap1ID = hap2_to_hap1[refID]
					hap2ID = refID
					hap1Len = fasta_1_len[hap1ID]
					hap2Len = fasta_2_len[hap2ID]
					plot_files[comparison]["Reports"][queryID] = {}
					plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name ,  dup_markers_file_name , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio )
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_pair_pdf_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio )
					#													make_pair_pdf_report(coords						  , coords_self                       , workdir						  , output_dir				     , queryID, refID, structure      , legacy                , markers               , dup_markers           , hap1ID , hap2ID , counts_hap1                        , min_align , similarity , ratio ) :
			elif comparison == "Hap2_vs_Hap1" :
				for hap2ID in sorted(hap2_ids.split(",")) :
					hap1ID = hap2_to_hap1[hap2ID]
					refID = hap1ID
					queryID = hap2ID
					hap1Len = fasta_1_len[hap1ID]
					hap2Len = fasta_2_len[hap2ID]
					plot_files[comparison]["Reports"][queryID] = {}
					plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name ,  dup_markers_file_name , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio )
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_pair_pdf_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , "diploid_gene_count_trace.hap1.txt", hit_len , hit_iden, unbalanced_ratio )
					#													make_pair_pdf_report(coords						  , coords_self                       , workdir						  , output_dir				     , queryID, refID, structure      , legacy                , markers               , dup_markers           , hap1ID , hap2ID , counts_hap1                        , min_align , similarity , ratio ) :
			elif comparison == "Hap1_vs_Hap2" :
				for hap1ID in sorted(hap1_ids.split(",")) :
					hap2ID = hap1_to_hap2[hap1ID]
					queryID = hap1ID
					refID = hap2ID
					hap1Len = fasta_1_len[hap1ID]
					hap2Len = fasta_2_len[hap2ID]
					plot_files[comparison]["Reports"][queryID] = {}
					plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name ,  dup_markers_file_name ,  hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio )
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_pair_pdf_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio )
					#													make_pair_pdf_report(coords						  , coords_self                       , workdir						  , output_dir				     , queryID, refID, structure      , legacy                , markers               , dup_markers           , hap1ID , hap2ID , counts_hap1                        , min_align , similarity , ratio ) :
			elif comparison == "Hap2_vs_Hap2" :
				for hap2ID in sorted(hap2_ids.split(",")) :
					hap1ID = hap2_to_hap1[hap2ID]
					refID = hap2ID
					queryID = hap2ID
					hap1Len = fasta_1_len[hap1ID]
					hap2Len = fasta_2_len[hap2ID]
					plot_files[comparison]["Reports"][queryID] = {}
					plot_files[comparison]["Reports"][queryID]["html"] = make_pair_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name ,  dup_markers_file_name , hap1ID , hap2ID , hap1Len , hap2Len , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio )
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_pair_pdf_report(os.path.basename(coords_file), os.path.basename(coords_file_self), os.path.realpath(haplodup_dir), os.path.realpath(outdir_name), queryID, refID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hap1ID , hap2ID , "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio )
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
					plot_files[comparison]["Reports"][queryID]["html"] = make_no_genes_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
			elif comparison == "Reference_vs_Hap1" :
				for queryID in sorted(ref_ids) :
					refID = ref_to_hap1[queryID]
					plot_files[comparison]["Reports"][queryID] = {}
					plot_files[comparison]["Reports"][queryID]["html"] = make_no_genes_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
			elif comparison == "Hap1_vs_Hap1" :
				continue
			elif comparison == "Hap2_vs_Reference" :
				for queryID in sorted(hap2_ids) :
					refID = hap2_to_ref[queryID]
					plot_files[comparison]["Reports"][queryID] = {}
					plot_files[comparison]["Reports"][queryID]["html"] = make_no_genes_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
			elif comparison == "Reference_vs_Hap2" :
				for queryID in sorted(ref_ids) :
					refID = ref_to_hap2[queryID]
					plot_files[comparison]["Reports"][queryID] = {}
					plot_files[comparison]["Reports"][queryID]["html"] = make_no_genes_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
			elif comparison == "Hap2_vs_Hap1" :
				for queryID in sorted(hap2_ids) :
					refID = hap2_to_hap1[queryID]
					plot_files[comparison]["Reports"][queryID] = {}
					plot_files[comparison]["Reports"][queryID]["html"] = make_no_genes_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
			elif comparison == "Hap1_vs_Hap2" :
				for queryID in sorted(hap1_ids) :
					refID = hap1_to_hap2[queryID]
					plot_files[comparison]["Reports"][queryID] = {}
					plot_files[comparison]["Reports"][queryID]["html"] = make_no_genes_html_report(os.path.basename(coords_file), os.path.basename(coords_file_self), haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
					plot_files[comparison]["Reports"][queryID]["pdf"] = make_no_genes_pdf_report( coords_file , coords_file_self, haplodup_dir, outdir_name, queryID, refID, hit_len , hit_iden)
			elif comparison == "Hap2_vs_Hap2" :
				continue
			else :
				print >> sys.stderr, "[ERROR] Report required for unknown comparison: " + comparison
				sys.exit(1)

	# Make Index
	html_index = make_index_from_report_db("index.html" , "." , haplodup_dir ,  plot_files  )


	# Windows with hotspots
	window_size = int(options.window)
	# On hit_counts_1 and hit_counts_2
	# hit_counts_1:
	# hit_counts_1[chr_hap1]= [
	#	...
	# 	[chr_hap1, start , end , id , h1_len , h2_len , ratio , descriptions , counts]
	#	...
	#	]
	if window_size > 0 :
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + '] === Searching for windows of genes with unbalanced count'
		print >> sys.stdout, '### Searching for windows of genes with unbalanced count'
		max_unbalanced = int(options.allowed)
		ratio_threshold = float(unbalanced_ratio)
		inverted_ratio_threshold = 1/float(unbalanced_ratio)

		hotspot_windows = {}
		for chr in sorted(hit_counts_1.keys()) :
			hotspot_windows[chr] = []
			gene_counts = sorted(hit_counts_1[chr])
			# Sorted by chr id (always the same) and then by gene start then by gene stop
			gene_counts_len = len(gene_counts)
			# Make windows
			for i in range(0,gene_counts_len + 1 - window_size ) :
				start = gene_counts[i][1]
				stop = gene_counts[i+window_size-1][2]
				count_off = 0
				for j in range(i , i + window_size ) :
					if int(gene_counts[j][4]) == 0 or int(gene_counts[j][5]) == 0 :
						if (int(gene_counts[j][4]) == 0 and int(gene_counts[j][5]) > 1 ) or (int(gene_counts[j][4])>1 and int(gene_counts[j][5]) == 0) :
							# One of the two counts is null and the other is 2 or more
							count_off += 1
					else :
						# None of the count is null, check the ratio
						ratio = float(gene_counts[j][5]) / float(gene_counts[j][4])
						if ratio <= ratio_threshold or ratio >= inverted_ratio_threshold :
							count_off += 1
				if count_off > max_unbalanced :
					hotspot_windows[chr].append( [ chr , int(start) , int(stop )] )
		hotspot_hap1_file = haplodup_dir + "/" + options.out + ".hotspots.windows.hap1.txt"
		hotspot_hap1 = open(hotspot_hap1_file , "w+")
		for chr in sorted(hotspot_windows.keys()) :
			for element in sorted(hotspot_windows[chr]) :
				print >> hotspot_hap1 , "\t".join([ str(x) for x in element ])
		hotspot_hap1.close()

		hotspot_windows = {}
		for chr in sorted(hit_counts_2.keys()) :
			hotspot_windows[chr] = []
			gene_counts = sorted(hit_counts_2[chr])
			# Sorted by chr id (always the same) and then by gene start then by gene stop
			gene_counts_len = len(gene_counts)
			# Make windows
			for i in range(0,gene_counts_len + 1 - window_size ) :
				start = gene_counts[i][1]
				stop = gene_counts[i+window_size-1][2]
				count_off = 0
				for j in range(i , i + window_size ) :
					if int(gene_counts[j][4]) == 0 or int(gene_counts[j][5]) == 0 :
						if (int(gene_counts[j][4]) == 0 and int(gene_counts[j][5]) > 1 ) or (int(gene_counts[j][4])>1 and int(gene_counts[j][5]) == 0) :
							# One of the two counts is null and the other is 2 or more
							count_off += 1
					else :
						# None of the count is null, check the ratio
						ratio = float(gene_counts[j][5]) / float(gene_counts[j][4])
						if ratio <= ratio_threshold or ratio >= inverted_ratio_threshold :
							count_off += 1
				if count_off > max_unbalanced :
					hotspot_windows[chr].append( [ chr , int(start) , int(stop )] )
		hotspot_hap2_file = haplodup_dir + "/" + options.out + ".hotspots.windows.hap2.txt"
		hotspot_hap2 = open(hotspot_hap2_file , "w+")
		for chr in sorted(hotspot_windows.keys()) :
			for element in sorted(hotspot_windows[chr]) :
				print >> hotspot_hap2 , "\t".join([ str(x) for x in element ])
		hotspot_hap2.close()


	if options.rejected :
		# Read rejection relationships
		print >> sys.stdout , '[' + str(datetime.datetime.now()) + '] = Running unused sequences QC analysis'
		print >> sys.stderr , '# Running unused sequences QC analysis'

		structure_comparison_dir = options.out + ".structure_comparison"
		mkdir(structure_comparison_dir)

		all_unused_by_chr = {}
		all_unused_to_chr = {}
		#for seq_id in all_unused_by_seq_id.keys() :
		#	chr_id , orientation = all_unused_by_seq_id[seq_id]
		#	if chr_id not in all_unused_by_chr :
		#		all_unused_by_chr[chr_id] = []
		#	all_unused_by_chr[chr_id].append(seq_id+"|"+orientation)
		for line in open(options.rejected , 'r') :
			seq_id , hapx_id , orientation = line.rstrip().split("\t")
			if hapx_id in hap1_to_chr :
				chr_id = hap1_to_chr[hapx_id]
			elif hapx_id in hap2_to_chr :
				chr_id = hap2_to_chr[hapx_id]
			elif hapx_id in ref_to_chr :
				chr_id = hap2_to_ref[hapx_id]
			else :
				print >> sys.stderr, "[WARNING] Sequence " + seq_id + " is associated to " + hapx_id + " pseudomolecule that is not present in correspondence file. Analysis cannot be performed"
				continue

			if chr_id not in all_unused_by_chr :
				all_unused_by_chr[chr_id] = []
			all_unused_by_chr[chr_id].append(seq_id+"|"+orientation)
			all_unused_to_chr[seq_id] = chr_id
			all_unused_to_chr[seq_id+"|"+orientation] = chr_id

		if options.marker_map and options.markers_hits :
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
						marker_map_by_seq[seq_id] = []
					marker_map_by_seq[seq_id].append([ int(pos) , marker_id ] )
					marker_map_by_id[marker_id] = [ seq_id , int(pos) , marker_id ]
		else :
			marker_map_by_seq = ""
			print >> sys.stderr, "[WARNING] Marker information missing"

		if options.legacy_groups :
			associated_legacy_seqid_file = structure_comparison_dir + "/legacy_components.association.tsv"
			# all_agp_db = legacy_agp + old_legacy_agp >> all sequences related to the components >> grouping on components
			associated_input_seqid_file , seq_group_db = make_seq_pair_from_groups( associated_legacy_seqid_file , options.legacy_groups , legacy_agp , hap1_ids.split(",") , hap2_ids.split(",") , fasta_dict.keys() , {} , "legacy" )
			# associated_input_seqid_file >> [ ... , [ "hap1_legacy" , Tid , Tstart , Tstop , "Query_legacy" , Qid , Qstart , Qstop , group_id] , ... ]
			#
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
					associated_input_seqid_file , seq_group_db = make_seq_pair_from_groups( associated_input_seqid_file ,options.input_group , agp_db , hap1_ids.split(",") , hap2_ids.split(",") , fasta_dict.keys() , fasta_len_dict , "input" )
					# associated_input_seqid_file >> [ ... , [ "hap1_HS" , Tid , Tstart , Tstop , "Query_HS" , Qid , Qstart , Qstop , group_id] , ... ]
			#else :
			#	# parse known_groups and unwanted_pairs
			#	associated_input_seqid_file = make_seq_pair_from_constrains(associated_input_seqid_file, known_groups, unwanted_pairs, alternative_sequences , agp_db, hap1_ids.split(","), hap2_ids.split(","), query_fasta_db.keys(), query_len, "input")
			#	# associated_input_seqid_file >> [ ... , [ "hap1_HS" , Tid , Tstart , Tstop , "Query_HS" , Qid , Qstart , Qstop , group_id] , ... ]

		reported_hits_on_seq = report_marker_usage( options.markers_hits , marker_map_by_seq , marker_map_by_id , agp_db , legacy_agp , hap1_to_chr , hap2_to_chr , all_unused_to_chr , structure_comparison_dir )
		# 	reported_hits_on_seq[seq_id]["id"] : seq_id
		# 	reported_hits_on_seq[seq_id]["chr"] : chr_id
		# 	reported_hits_on_seq[seq_id]["markers"] : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ]
		# 	reported_hits_on_seq[seq_id]["range"] : [marker_pos_min , marker_pos_max] ,
		# 	reported_hits_on_seq[seq_id]["orientation"] : ["+" or "-" or "."] }
		print >> sys.stdout, '[' + str(datetime.datetime.now()) + "] = Reading markers hits"
		structure_plot_db = {"Rejected" : {}}
		for chr_id in sorted(chr_ids) :
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
					structure_plot_db["Rejected"][chr_id][seq_id] = rejected_QC( structure_comparison_dir , seq_id , fasta_dict , chr_id , fasta_db_1 , fasta_db_2 , chr_to_hap1 , chr_to_hap2 , coord_tables["Hap2_vs_Hap1"][0] , associated_input_seqid_file , associated_legacy_seqid_file , agp_db , legacy_agp , "" , seq_group_db , markers_db , reported_hits_on_seq , marker_map_by_seq , options.cores , paths)

		rejected_index_file_name = "index.rejected_sequences.html"
		rejected_index_file_full_path = structure_comparison_dir + "/" + rejected_index_file_name
		rejected_index_file_full_path = make_index_from_report_db(rejected_index_file_name , "." , structure_comparison_dir ,  structure_plot_db  )




	print >> sys.stdout , "------------------------------"
	print >> sys.stdout , "- Done"
	print >> sys.stdout , "------------------------------"
	print >> sys.stderr , "##############################"
	print >> sys.stderr , "# Done"
	print >> sys.stderr , "##############################"



if __name__ == '__main__':
	main()