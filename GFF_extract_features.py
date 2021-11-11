#!/usr/bin/env python

from lib_files.GFF_lib import *
from lib_files.FASTA_lib import *

gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


def main():
	#### Main

	parser = OptionParser()
	parser.add_option("-g", "--genome", dest="genome",
					help="Genome sequence in FASTA format", metavar="genome.fasta [Required]")
	parser.add_option("-a", "--annot", dest="gff",
					help="GFF3 annotation to check", metavar="annotation.gff3 [Required]")
	parser.add_option("-p", "--prefix", dest="prefix",  default="out",
					help="Prefix to output file names", metavar="out")
	parser.add_option("-o", "--out-gff3", dest="print_gff3", default=True, action="store_false",
					help="Avoid printing an updated gff3 file")
	parser.add_option("-l", "--lengths", dest="print_len", default=True, action="store_false",
					help="Avoid printing features length lists")
	parser.add_option("-s", "--sequences", dest="print_seq", default=True, action="store_false",
					help="Avoid printing feature sequences in FASTA format")
	parser.add_option("-c", "--counts", dest="print_counts", default=True, action="store_false",
					help="Avoid printing feature counts per mRNA/locus")
	parser.add_option("-d", "--descriptive-stats", dest="print_stats", default=True, action="store_false",
					help="Avoid printing feature descriptive statistics")
	parser.add_option("-m", "--mono-and-multi", dest="print_mems", default=True, action="store_false",
					help="Avoid printing separate information/files for mono and multi-exonic mRNAs")
	parser.add_option("-i", "--intermediate", dest="temp", default=True, action="store_false",
					help="Avoid printing intermediate gff3 files")
	parser.add_option("-n", "--plot", dest="plot_distributions", default=True, action="store_false",
					help="Avoid plotting length and count distribution (pdf file)")
	parser.add_option("--nopseudo", dest="pseudo", default=True, action="store_false",
					help="Avoid writing pseudogene information files")

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)

	(options, args) = parser.parse_args()

	if not options.genome :
		print >> sys.stderr , "[ERROR] Genome FASTA file missing"
		sys.exit(1)
	if not options.gff :
		print >> sys.stderr , "[ERROR] Genome annotation GFF3 file missing"
		sys.exit(1)

	## Read reference
	reference = dict((seq.id, seq.upper()) for seq in SeqIO.parse(open(options.genome), "fasta"))
	reference_len = get_fasta_lengths_from_file(options.genome)
	genes = {}
	mRNA = {}

	gene_region_length = []
	mRNA_region_length = []

	exon_region_length = []
	UTRp5_region_length = []
	CDS_region_length = []
	UTRp3_region_length = []
	intron_region_length = []
	intergenic_region_length = []
	mRNA_length = []
	UTRp5_length = []
	UTRp3_length = []
	CDS_length = []
	protein_length = []

	gene_region_length_singleExon = []
	mRNA_region_length_singleExon = []

	exon_region_length_singleExon = []
	UTRp5_region_length_singleExon = []
	CDS_region_length_singleExon = []
	UTRp3_region_length_singleExon = []

	mRNA_length_singleExon = []
	UTRp5_length_singleExon = []
	UTRp3_length_singleExon = []
	CDS_length_singleExon = []
	protein_length_singleExon = []


	gene_region_length_multiExon = []
	mRNA_region_length_multiExon = []

	exon_region_length_multiExon = []
	UTRp5_region_length_multiExon = []
	CDS_region_length_multiExon = []
	UTRp3_region_length_multiExon = []
	intron_region_length_multiExon = []

	mRNA_length_multiExon = []
	UTRp5_length_multiExon = []
	UTRp3_length_multiExon = []
	CDS_length_multiExon = []
	protein_length_multiExon = []

	mRNA_per_gene = []
	mRNA_per_gene_multiExon = []
	mRNA_per_gene_singleExon = []

	exons_per_mRNA = []
	UTRp5_exons_per_mRNA = []
	UTRp3_exons_per_mRNA = []
	CDS_exons_per_mRNA = []

	exons_per_mRNA_multiExon  = []
	UTRp5_exons_per_mRNA_multiExon  = []
	UTRp3_exons_per_mRNA_multiExon  = []
	CDS_exons_per_mRNA_multiExon  = []



	mRNA_fasta = {}
	CDS_fasta = {}
	protein_fasta = {}
	pseudo_fasta = {}

	mRNA_singleExon_fasta = {}
	CDS_singleExon_fasta = {}
	protein_singleExon_fasta = {}

	mRNA_multiExon_fasta = {}
	CDS_multiExon_fasta = {}
	protein_multiExon_fasta = {}


	mRNA_feat_order = {"exon": 1, "intron":1, "five_prime_UTR":2 , "CDS":3, "three_prime_UTR":4}

	genes , mRNA = read_gff3(options.gff)
	### genes structure:
	### genes[gene_id] = {}

	if options.temp :
		random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(8))
		tmp_gff3 = options.prefix + ".tmp" + random_string + ".gff3"
		write_gff3( genes , tmp_gff3 , reference_len )
		print >> sys.stderr , "tmp gff3 file generated"

	intergenic_id = 1
	prev_chr = ""
	intergenic_start = ""

	error_exon_file = open(options.prefix+".error_exon_file.txt", 'w')

	### Add introns, check UTRs, check CDS phase, create sequences
	## For each gene, sorted by chr : start position

	Dropped_gene_file = open(options.prefix+".dropped.genes.txt",'w')
	Dropped_gene_list = []
	Dropped_mRNA_file = open(options.prefix+".dropped.mRNAs.txt",'w')

	no_cds_db = {}

	if options.pseudo :
		Dropped_models_file = open(options.prefix+".pseudogenes.gff3",'w')
		Dropped_models_fasta = options.prefix+".pseudogenes.mRNA.fasta"
		print >> Dropped_models_file, "##gff-version 3"

	for g_key, g_value in sorted(genes.iteritems(), key=lambda (k, v): itemgetter(1, 2)(v)) :
		chr_printed = ""
		gene_printed = ""
		try :
			seqname, source, feature, start, end, score, strand, frame, attribute = genes[g_key][0].rstrip().split("\t")
		except:
			print >> sys.stderr , "[ERROR] line for gene " + g_key + ": not declared explicitly, only subfeatures are present."
			#print >> sys.stderr ,  genes[g_key]
			#print >> sys.stderr ,  genes[g_key][0].rstrip()
			sys.exit(2)
		#add intergenic regions
		print >> sys.stderr , "-----------  Gene: " + g_key
		#print >> sys.stderr , str(genes[g_key][3].keys())
		#print >> sys.stderr , "----------- " + prev_chr



		# For each transcript
		loc_mRNA_count = 0

		exons_per_mRNA_local = []

		new_gene_start = ""
		new_gene_end = ""

		for m_key, m_value in sorted(genes[g_key][3].iteritems(), key=lambda (k, v): itemgetter(1)(v)) :
			print >> sys.stderr, "------ mRNA: " + m_key

			if not 3 in genes[g_key][3][m_key][2] :
				# No CDS annotated, drop the transcript
				print >> sys.stderr, "---- " + m_key + " -- dropped -- no CDS annotation provided"
				print >> Dropped_mRNA_file, m_key
				if g_key not in no_cds_db :
					no_cds_db[g_key] = genes[g_key][:]
					no_cds_db[g_key][3] = {}
				try :
					no_cds_db[g_key][3][m_key] = genes[g_key][3][m_key][:]
				except :
					print >> sys.stderr, genes[g_key][3].keys()
					sys.exit(2)

				del genes[g_key][3][m_key]
				continue

			DROPPED = False

			# m_key -> mRNA name
			# feat_key = 1=exon | intron->1 | 2=five_prime_UTR | 3=CDS | 4=three_prime_UTR"
			cds_count = 0
			### Correct CDSs phase and extract whole CDS sequence
			offset = 0
			CDS_id = 0
			strand=genes[g_key][3][m_key][2][3][0][0].split("\t")[6]
			phase = ""
			CDS_seq= Seq("")
			CDS_seq.id = ""


			# TO TEST : coordinates -> correct start and stop of extreme exons and CDSs to fit inside the sequence, drop models with CDS outside boundaries
			# Test CDS, if one has negative start or stops outside of the sequence, drop the mRNA
			for CDS_feat in genes[g_key][3][m_key][2][3] :
				seqname, source, feature, start, end, score, strand, cds_phase, attribute = CDS_feat[0].split("\t")
				seq_length = len(reference[seqname])
				if int(start) < 1 or int(end) > seq_length :
					DROPPED = True
					print >> sys.stderr, "---- " + m_key + " -- dropped -- CDS extending beyond sequence boundaries"
					print >> Dropped_mRNA_file, m_key
					del genes[g_key][3][m_key]
					break

			if DROPPED :
				continue
			else:
				new_exon_list = []
				# Test exon features, if external to the sequence, trim it to stay inside or drop the feature
				for exon_feat in genes[g_key][3][m_key][2][1] :
					seqname, source, feature, start, end, score, strand, cds_phase, attribute = exon_feat[0].split("\t")
					seq_length = len(reference[seqname])
					new_start = max(int(start) , 0)
					new_end = min(int(end), (seq_length+1))

					if (new_start == 0) or (new_end == (seq_length+1) ):
						if (new_start == 1) and (new_end == seq_length) :
							# Drop exon
							print >> sys.stderr, "----- " + m_key + " -- exon dropped -- outside of sequence boundaries"
						else :
							# Edit coords
							new_exon_line = "\t".join([str(x) for x in [seqname, source, feature, new_start, new_end, score, strand, cds_phase, attribute] ])
							new_exon_feat = [ new_exon_line , int(new_start) , int(new_end) ]
							print >> sys.stderr, "----- " + m_key + " -- exon edited -- extending beyond sequence boundaries"
							new_exon_list.append(new_exon_feat)
					else :
						# Good exon
						new_exon_list.append(exon_feat)
				genes[g_key][3][m_key][2][1] = new_exon_list

			# check if phase info is missing, in case that assume the first CDS has phase=0 (frame starts from first base) and uodate other CDS accordingly
			phase_missing = False

			for el in sorted(genes[g_key][3][m_key][2][3], key=lambda x: x[1]) :
				seqname, source, feature, start, end, score, strand, cds_phase, attribute = el[0].split("\t")
				try:
					value = int(cds_phase)
				except ValueError:
					phase_missing = True

			if phase_missing :
				print >> sys.stderr, "---- Some CDS missing phase, correcting "
				actual_seq_length = 0
				new_cds = []

				# phase is missing for at least one CDS, recostruct the whole phasing information

				if strand == "+" :
					for el in sorted(genes[g_key][3][m_key][2][3], key=lambda x: x[1]) :
						seqname, source, feature, start, end, score, strand, tmp_phase, attribute = el[0].split("\t")
						#print >> sys.stderr, "-- Old line " + el[0]
						#print >> sys.stderr, "-- Length of of the CDS so far: " + str(actual_seq_length)

						tmp_phase = 3 - actual_seq_length%3
						if tmp_phase == 3 :
							tmp_phase = 0

						new_line = "\t".join([seqname, source, feature, start, end, score, strand, str(tmp_phase), attribute])
						#print >> sys.stderr, "-- New line " + new_line
						new_cds.append([ new_line , int(start), int(end) ])
						actual_seq_length += int(end) - int(start) + 1
				else :
					for el in sorted(genes[g_key][3][m_key][2][3], key=lambda x: x[1], reverse=True) :
						seqname, source, feature, start, end, score, strand, tmp_phase, attribute = el[0].split("\t")
						#print >> sys.stderr, "-- Old line " + el[0]
						#print >> sys.stderr, "-- Length of of the CDS so far: " + str(actual_seq_length)

						tmp_phase = 3 - actual_seq_length%3
						if tmp_phase == 3 :
							tmp_phase = 0

						new_line = "\t".join([seqname, source, feature, start, end, score, strand, str(tmp_phase), attribute])
						#print >> sys.stderr, "-- Old line " + new_line
						new_cds.append([ new_line , int(start), int(end) ])
						actual_seq_length += int(end) - int(start) + 1


				genes[g_key][3][m_key][2][3]=new_cds


			# Check phase correctness

			for el in sorted(genes[g_key][3][m_key][2][3], key=lambda x: x[1]) :
				cds_count += 1
				#print >> sys.stderr, "-- CDS: " + el[0]
				seqname, source, feature, start, end, score, strand, old_phase, attribute = el[0].split("\t")
				#print >> sys.stderr, el[0]
				if not strand=="-":
					CDS_seq += reference[seqname][int(start)-1:int(end)]
					if phase == "" : phase = int(old_phase)
				else :
					CDS_seq = reference[seqname][int(start)-1:int(end)].reverse_complement() + CDS_seq
					phase = int(old_phase)
				#print >> sys.stderr, reference[seqname][int(start)-1:int(end)].seq
				#print >> sys.stderr, CDS_seq.seq

			loc_CDS_exons_per_mRNA=cds_count

			for new_phase in (0,1,2,3) :
				new_seq = CDS_seq[new_phase:]
				first_codon=new_seq[0:3].seq
				#print >> sys.stderr , "Phase " + str(new_phase) +  ": First codon: " + first_codon
				#print >> sys.stderr , str(new_seq.seq)
				prot_sequence=Seq(str(new_seq.seq)).translate()
				#print >> sys.stderr , " Protein: "+ prot_sequence

				if first_codon=="ATG" :
					# Found a start codon with new_phase, check for in-frame stop codons
					pieces=len(prot_sequence.rstrip("*").split("*"))
					#print >> sys.stderr , pieces

					if pieces==1 :
						print >> sys.stderr , "---- Good first CDS phase: " + str(new_phase) + ";  no stop codons found in frame"
						break

			if new_phase==3 :
				#print >> sys.stdout, m_key + " dropped"
				print >> sys.stderr, "---- " + m_key + "-- dropped"
				print >> sys.stderr, "-- " + m_key + " : " + CDS_seq.seq
				print >> sys.stderr, "-- " + m_key + " : " + prot_sequence
				DROPPED=True
			else :
				offset = new_phase - phase
				print >> sys.stderr , "---- Phase offset detected: " +str(offset)
				#new_phase = old_pahse + offset

			new_CDS_list = []
			loc_CDS_region_length = []


			if DROPPED :
				# Delete mRNA with faulty CDS name in dropped mRNAs list
				print >> Dropped_mRNA_file, m_key
				# Print the mRNA structure as pseudogene

				mRNA_seq = ""
				actual_chr=genes[g_key][1]
				pseudo_name=g_key + ".pseudo"
				#print >> sys.stderr, g_key
				#print >> sys.stderr, genes[g_key][0]
				#print >> sys.stderr, genes[g_key][1]
				#print >> sys.stderr, genes[g_key][2]
				if chr_printed != actual_chr :
					if options.pseudo :
						print >> Dropped_models_file, "##sequence-region " + actual_chr + " 1 " + str(len(reference[actual_chr]))
					chr_printed = actual_chr
				if not g_key == gene_printed :
					if options.pseudo :
						print >> Dropped_models_file, "### " + pseudo_name
					gene_printed = g_key

					seqname, source, feature, start, end, score, strand, old_phase, attribute = genes[g_key][0].split("\t")
					attribute = "ID=" + pseudo_name + ";pseudogene=\"unknown\";pseudo=true"
					if options.pseudo :
						print >> Dropped_models_file, "\t".join([ str(x) for x in [seqname, source, feature, start, end, score, strand, old_phase, attribute] ])

				seqname, source, feature, start, end, score, strand, old_phase, attribute = genes[g_key][3][m_key][0].split("\t")
				attribute = "ID=" + m_key + ";Parent="+ pseudo_name + ";pseudogene=\"unknown\";pseudo=true"
				feature = "transcript"
				if options.pseudo :
					print >> Dropped_models_file, "\t".join([ str(x) for x in [seqname, source, feature, start, end, score, strand, old_phase, attribute] ])

				for feat_key in sorted(genes[g_key][3][m_key][2].keys()):
					for el in sorted(genes[g_key][3][m_key][2][feat_key], key=lambda x: x[1]) :
						seqname, source, feature, start, end, score, strand, old_phase, attribute = el[0].split("\t")
						attribute += ";pseudogene=\"unknown\";pseudo=true"
						if feature=="exon" :
							if not strand=="-":
								mRNA_seq += reference[seqname][int(start)-1:int(end)]
							else :
								mRNA_seq  = reference[seqname][int(start)-1:int(end)].reverse_complement() + mRNA_seq
							if options.pseudo :
								print >> Dropped_models_file, "\t".join([ str(x) for x in [seqname, source, feature, start, end, score, strand, old_phase, attribute] ])

				pseudo_fasta[m_key]=str(mRNA_seq)
				if options.pseudo :
					print >> Dropped_models_file, "###"
				del genes[g_key][3][m_key]
				continue
			else :
				loc_CDS_seq = Seq("")
				loc_CDS_seq.id = ""
				phase = ""
				for el in sorted(genes[g_key][3][m_key][2][3], key=lambda x: x[1]) :
					#print >> sys.stderr, "### OLD line: " + el[0]
					seqname, source, feature, start, end, score, strand, old_phase, attribute = el[0].split("\t")

					if int(old_phase) + offset == 3 : new_phase=0
					elif int(old_phase) + offset == 4 : new_phase=1
					elif int(old_phase) + offset == -1 : new_phase=2
					elif int(old_phase) + offset == -2 : new_phase=1
					else : new_phase = int(old_phase) + offset
					CDS_new_line= "\t".join([seqname, source, feature, start, end, score, strand, str(new_phase), attribute])
					new_CDS_list.append([CDS_new_line,int(start),int(end)])
					loc_CDS_region_length.append(int(end) - int(start) + 1)
					#print >> sys.stderr, CDS_new_line
					if not strand=="-":
						if phase == "" : phase=new_phase
						#print >> sys.stderr , reference[seqname][int(start)-1:int(end)].seq
						loc_CDS_seq += reference[seqname][int(start)-1:int(end)]
						#print >> sys.stderr , loc_CDS_seq.seq
						#print >> sys.stderr , ""
					else :
						phase=new_phase
						#print >> sys.stderr , reference[seqname][int(start)-1:int(end)].reverse_complement().seq
						loc_CDS_seq = reference[seqname][int(start)-1:int(end)].reverse_complement() + loc_CDS_seq
						#print >> sys.stderr , loc_CDS_seq.seq
						#print >> sys.stderr , ""

				loc_CDS_seq = loc_CDS_seq[phase:]
				#print >> sys.stderr , loc_CDS_seq.seq

				genes[g_key][3][m_key][2][3]=new_CDS_list

				loc_prot_sequence=Seq(str(loc_CDS_seq.seq)).translate()

				### STEP: Add missing exons, introns, extract correct mRNA annotation and extract sequence sequence
				# For each exon:
				mRNA_seq = ""
				loc_intron_region_length = []
				loc_exon_region_length = []
				#print >> sys.stderr , genes[g_key][3][m_key][2]
				new_mRNA_start = ""
				new_mRNA_end = ""
				new_exons_list = []

				# If exons are missing in the annotation, create them from CDSs and UTRs
				if not 1 in genes[g_key][3][m_key][2] :
					print >> sys.stderr , "---- Adding exons to " + m_key
					# Exon annotation missing, create it from the CDSs and, eventually, UTRs
					# 1) Add regions from CDS
					for feat in sorted(genes[g_key][3][m_key][2][3], key=lambda x: x[1]) :
						seqname, source, feature, start, end, score, strand, frame, attribute = feat[0].rstrip().split("\t")
						new_exon_line="\t".join([seqname, source, "exon", start, end, score, strand, "." , "name" ])
						new_exons_list.append([new_exon_line, int(start), int(end)])

					# 2) Compare to 5p UTRs, if exist
					if 2 in genes[g_key][3][m_key][2] :
						for feat in sorted(genes[g_key][3][m_key][2][2], key=lambda x: x[1]) :
							for i in range( 0 , len(new_exons_list) ) :
								new_exon = new_exons_list[i]
								if new_exon[2] == feat[1]-1 :
									#UTR end joint to CDS start (UTR left of CDS)
									new_start = str(feat[1])
									seqname, source, feature, start, end, score, strand, frame, attribute = new_exon[0].rstrip().split("\t")
									new_line = "\t".join([seqname, source, feature, new_start, end, score, strand, frame, "name"])
									new_exons_list[i] = [new_line,int(new_start),int(end)]
								elif new_exon[1] == feat[1]+1 :
									# UTR start joint to CDS ens (UTR right of CDS)
									new_end = str(feat[2])
									seqname, source, feature, start, end, score, strand, frame, attribute = new_exon[0].rstrip().split("\t")
									new_line = "\t".join([seqname, source, feature, start, new_end, score, strand, frame, "name"])
									new_exons_list[i] = [new_line,int(new_start),int(end)]
								else :
									#UTR in an independent exon
									seqname, source, feature, start, end, score, strand, frame, attribute = feat[0].rstrip().split("\t")
									new_line = "\t".join([seqname, source, "exon", start, end, score, strand, frame, "name"])
									new_exons_list.append([new_line,int(start),int(end)])

					# 2) Compare to 3p UTRs, if exist
					if 4 in genes[g_key][3][m_key][2] :
						for feat in sorted(genes[g_key][3][m_key][2][4], key=lambda x: x[1]) :
							for i in range( 0 , len(new_exons_list) ) :
								new_exon = new_exons_list[i]
								if new_exon[2] == feat[1]-1 :
									# UTR end joint to CDS start (UTR left of CDS)
									new_start = str(feat[1])
									seqname, source, feature, start, end, score, strand, frame, attribute = new_exon[0].rstrip().split("\t")
									new_line = "\t".join([seqname, source, feature, new_start, end, score, strand, frame, "name"])
									new_exons_list[i] = [new_line,int(new_start),int(end)]
								elif new_exon[1] == feat[1]+1 :
									# UTR start joint to CDS ens (UTR right of CDS)
									new_end = str(feat[2])
									seqname, source, feature, start, end, score, strand, frame, attribute = new_exon[0].rstrip().split("\t")
									new_line = "\t".join([seqname, source, feature, start, new_end, score, strand, frame, "name"])
									new_exons_list[i] = [new_line,int(new_start),int(end)]
								else :
									# UTR in an independent exon
									seqname, source, feature, start, end, score, strand, frame, attribute = feat[0].rstrip().split("\t")
									new_line = "\t".join([seqname, source, "exon", start, end, score, strand, frame, "name"])
									new_exons_list.append([new_line,int(start),int(end)])

				else :
					for feat in sorted(genes[g_key][3][m_key][2][1], key=lambda x: x[1]) :
						try :
							seqname, source, feature, start, end, score, strand, frame, attribute = feat[0].rstrip().split("\t")
						except :
							print >> sys.stderr, feat
							exit(1)
						new_exon_line = "\t".join([seqname, source, "exon", start, end, score, strand, "." , attribute])
						new_exons_list.append([new_exon_line, int(start), int(end)])

				# Update exon and mRNA annotation for m_key
				genes[g_key][3][m_key][2][1]=[]
				exon_count=1
				for new_exon in sorted(new_exons_list, key=lambda x: x[1]) :
					#print >> sys.stderr , "------ " + new_exon[0]
					exon_seq = str( reference[seqname][ int(new_exon[1])-1 : int(new_exon[2]) ].seq )
					#print >> sys.stderr , "------ Exon seq: " + exon_seq
					seqname, source, feature, start, end, score, strand, frame, attribute = new_exon[0].rstrip().split("\t")
					attribute = "ID=" + m_key + ".exon_" + str(exon_count) + ";Parent=" + m_key
					# Test if exon is ending in gap and eventually
					exon_len = len(exon_seq)
					#print >> sys.stderr , "------ Original length: " + str(exon_len) + " | rstrip: " + str(len(exon_seq.rstrip("N"))) + " | lstrip: " + str(len(exon_seq.lstrip("N")))

					if len(exon_seq.rstrip("N")) < exon_len :
						# Falling in a gap on the right
						# correct end coordinates
						corr_end = int(new_exon[2]) - ( exon_len - len(exon_seq.rstrip("N")) )
						exon_seq = exon_seq.rstrip("N")
						exon_len = len(exon_seq)
						print >> sys.stderr , "---- Editing exons end coordinate (gap) ID=" + m_key + ".exon_" + str(exon_count)
						#print >> sys.stderr , exon_seq == str( reference[seqname][ int(new_exon[1])-1 : int(corr_end) ].seq )
					else :
						corr_end = int(new_exon[2])

					if len(exon_seq.lstrip("N")) < exon_len :
						# correct coordinates
						corr_start = int(new_exon[1]) + ( exon_len - len(exon_seq.lstrip("N")) - 1 )
						print >> sys.stderr , "---- Editing exons start coordinate (gap) ID=" + m_key + ".exon_" + str(exon_count)
						#print >> sys.stderr , exon_seq.lstrip("N") == str( reference[seqname][ int(corr_start) : int(corr_end) ].seq )
					else:
						corr_start = int(new_exon[1])

					if new_mRNA_start == "" or new_mRNA_start > int(corr_start) :
						new_mRNA_start = int(corr_start)
					if new_mRNA_end == "" or new_mRNA_end < int(corr_end) :
						new_mRNA_end = int(corr_end)

					new_exon_line = "\t".join([seqname, source, feature, str(corr_start), str(corr_end), score, strand, frame, attribute])
					genes[g_key][3][m_key][2][1].append([ new_exon_line , corr_start , corr_end ])
					exon_count+=1

				#print >> sys.stderr , genes[g_key][3][m_key][0]
				#print >> sys.stderr , genes[g_key][3][m_key][2]

				exon_count = 0
				prev_stop = ""
				intron_count = 0

				for feat in sorted(genes[g_key][3][m_key][2][1], key=lambda x: x[1]) :
					# Create intron if needed
					# recreate (to eventually correct) the mRNA coordinates
					#print >> sys.stderr, feat[0]
					exon_count += 1
					try:
						seqname, source, feature, start, end, score, strand, frame, attribute = feat[0].rstrip().split("\t")
					except :
						print >> sys.stderr, feat
					loc_exon_region_length.append(feat[2]-feat[1]+1)

					if not strand=="-":
						mRNA_seq += reference[seqname][feat[1]-1:feat[2]]
						#print >> sys.stderr, reference[seqname][feat[1]-1:feat[2]].seq
						#print >> sys.stderr, mRNA_seq.seq
					else :
						mRNA_seq  = reference[seqname][feat[1]-1:feat[2]].reverse_complement() + mRNA_seq
						#print >> sys.stderr, reference[seqname][feat[1]-1:feat[2]].reverse_complement().seq
						#print >> sys.stderr, mRNA_seq.seq

					if not prev_stop=="" :
						# 2nd exon or more
						intron_line = "\t".join([seqname, source, "intron", str(prev_stop+1), str(feat[1]-1), ".", strand, ".", "ID="+m_key+".intron."+str(intron_count)+";Parent="+m_key])
						genes[g_key][3][m_key][2][1].append([intron_line, int(prev_stop)+1, int(feat[1])-1])
						#print >> sys.stderr, intron_line
						loc_intron_region_length.append(int(feat[1])-int(prev_stop)-1)
					prev_stop=feat[2]
					intron_count+=1
				#print >> sys.stderr, loc_intron_region_length
				#print >> sys.stderr, mRNA_seq.seq
				loc_exons_per_mRNA = exon_count

				seqname, source, feature, new_mrna_start, new_mrna_end, score, strand, frame, attribute = genes[g_key][3][m_key][0].split("\t")

				loc_mRNA_count += 1
				loc_mRNA_region_length = int(new_mrna_end)-int(new_mrna_start) + 1

				if ( new_gene_start == "" ) or ( new_gene_start > int(new_mrna_start) ) :
					new_gene_start = int(new_mrna_start)
				if (new_gene_end == "" ) or ( new_gene_end < int(new_mrna_end)) :
					new_gene_end = int(new_mrna_end)

				### STEP: Create UTRs (drop if existing)
				#print >> sys.stderr, "mrna - " + mRNA_seq.seq
				#print >> sys.stderr, "CDS -  " + loc_CDS_seq.seq
				#print >> sys.stderr, str(mRNA_seq.seq).split(str(loc_CDS_seq.seq))

				loc_UTRp5_region_length = []
				loc_UTRp3_region_length = []
				loc_UTRp5_length = []
				loc_UTRp3_length = []
				loc_UTRp5_exons_per_mRNA = []
				loc_UTRp3_exons_per_mRNA = []

				if mRNA_seq.seq == loc_CDS_seq.seq :
					print >> sys.stderr, "---- No UTRs"
					# NO UTR regions in the annotation mRNAs == CDSs
				else :
					#UTRp5_seq , UTRp3_seq = mRNA_seq.seq.split(loc_CDS_seq.seq)
					try :
						UTRp5_seq , UTRp3_seq = str(mRNA_seq.seq).split(str(loc_CDS_seq.seq))
					except ValueError :
						print >> error_exon_file, m_key
						print >> sys.stderr , "---- " + m_key + " dropped: exon structure incompatible with CDS"
						print >> Dropped_mRNA_file, m_key
						del genes[g_key][3][m_key]
						continue
					#print >> sys.stderr, "---- Updating UTR information "
					#print >> sys.stderr, UTRp5_seq
					#print >> sys.stderr, len(UTRp5_seq)
					#print >> sys.stderr, UTRp3_seq
					#print >> sys.stderr, len(UTRp3_seq)
					if len(UTRp5_seq) > 0 : loc_UTRp5_length=[len(UTRp5_seq)]
					if len(UTRp3_seq) > 0 : loc_UTRp3_length=[len(UTRp3_seq)]
					#print >> sys.stderr, loc_UTRp5_length
					#print >> sys.stderr, loc_UTRp3_length

					# if 5pUTRs and 3pUTR are not annotated, add annotation
					# For each exon
					# if len(UTRp5_seq) > 0 -> I have to add 5pUTR exons
					# if len(UTRp3_seq) > 0 -> I have to add 3pUTR exons
					# If strand == "+" -> starts with 5pUTR; else start with 3pUTR

					if 2 in genes[g_key][3][m_key][2]: del genes[g_key][3][m_key][2][2]
					if 4 in genes[g_key][3][m_key][2]: del genes[g_key][3][m_key][2][4]
					utr_count = 1

					if strand == "+" :
						prev_utr = "five_prime_UTR"
					else :
						prev_utr = "three_prime_UTR"

					CDS_found=False

					for exon_feat in sorted(genes[g_key][3][m_key][2][1], key=lambda x: x[1]) :
						# Check if any CDS falls in the exon
						seqname, source, feature, start, end, score, strand, old_phase, attribute = exon_feat[0].split("\t")

						if feature == "intron" : continue # create UTRs only where exons are

						#print >> sys.stderr , exon_feat[0]
						exon_range = "_".join(str(c) for c in range(exon_feat[1],exon_feat[2]+1))
						UTR_region = ""

						for CDS_feat in sorted(genes[g_key][3][m_key][2][3], key=lambda x: x[1]) :
							#print >> sys.stderr , CDS_feat[0]
							CDS_range = "_".join(str(c) for c in range(CDS_feat[1],CDS_feat[2]+1))

							region_difference = [ el for el in exon_range.split(CDS_range) if el != '' ]

							#print >> sys.stderr , exon_range.split(CDS_range)
							#print >> sys.stderr , region_difference

							#print >> sys.stderr , len(region_difference)

							if region_difference != [exon_range] :
								# CDS found
								#print >> sys.stderr , "CDS found !"
								#print >> sys.stderr , CDS_feat[0]
								break


						if region_difference == [exon_range] :
							# No CDS found overlapping, whole UTR
							feature = prev_utr
							print >> sys.stderr , "-- Add " + feature + " feature form a whole exon region"
							UTR_line = "\t".join([seqname, source, feature, start, end, "." , strand, old_phase, "ID="+m_key+"."+feature+"."+str(utr_count)+";Parent="+m_key])
							#print >> sys.stderr , UTR_line
							feat_code = mRNA_feat_order[feature]
							if not feat_code in genes[g_key][3][m_key][2] : genes[g_key][3][m_key][2][feat_code] = []
							genes[g_key][3][m_key][2][feat_code].append([UTR_line,int(start),int(end)])

						else :
							# exon has CDS on inside
							if region_difference == [] :
								#print >> sys.stderr , " exon == CDS "
								# CDS is the whole exon, switch feature UTRp5->UTRp3 or UTRp3->UTRp5
								utr_count = 0
								if not CDS_found :
									CDS_found=True
									if strand == "+" : prev_utr = "three_prime_UTR"
									else : prev_utr = "five_prime_UTR"

							else :
								if len(region_difference)==1 :
									# Exon partially UTR on one side only, create one UTR region
									# Identify if 5p or 3p
									UTR_start = region_difference[0].rstrip("_").lstrip("_").split("_")[0]
									UTR_end = region_difference[0].rstrip("_").lstrip("_").split("_")[-1]
									# Identify if 5p or 3p
									# exon and UTR share start coordinate -> UTR on leftmost part
									# exon and UTR share end coordinate -> UTR on rightmost part
									if UTR_start == start and strand == "+" :
										feature = "five_prime_UTR"
									elif UTR_end == end and strand == "-" :
										utr_count = 1
										feature = "five_prime_UTR"
									elif UTR_start == start and strand == "-" :
										feature = "three_prime_UTR"
									else :
										utr_count = 1
										feature = "three_prime_UTR"

									print >> sys.stderr , "-- Add " + feature + " feature form part of an exon "
									UTR_line = "\t".join([seqname, source, feature, UTR_start, UTR_end, "." , strand, old_phase, "ID="+m_key+"."+feature+"."+str(utr_count)+";Parent="+m_key])
									#print >> sys.stderr , UTR_line
									feat_code = mRNA_feat_order[feature]
									if not feat_code in genes[g_key][3][m_key][2] : genes[g_key][3][m_key][2][feat_code] = []
									genes[g_key][3][m_key][2][feat_code].append([UTR_line,int(UTR_start),int(UTR_end)])

									if strand == "+" : prev_utr = "three_prime_UTR"
									else : prev_utr = "five_prime_UTR"

								else :
									# Exon partially UTR on both sides, extract 2 UTR regions
									# 1st UTR
									UTR_start = region_difference[0].rstrip("_").lstrip("_").split("_")[0]
									UTR_end = region_difference[0].rstrip("_").lstrip("_").split("_")[-1]

									if strand == "+" : feature = "five_prime_UTR"
									else : feature = "three_prime_UTR"

									print >> sys.stderr , "-- Add " + feature + " feature form the upstream part of an exon "
									UTR_line = "\t".join([seqname, source, feature, UTR_start, UTR_end, "." , strand, old_phase, "ID="+m_key+"."+feature+"."+str(utr_count)+";Parent="+m_key])
									#print >> sys.stderr , UTR_line
									feat_code = mRNA_feat_order[feature]
									if not feat_code in genes[g_key][3][m_key][2] : genes[g_key][3][m_key][2][feat_code] = []
									genes[g_key][3][m_key][2][feat_code].append([UTR_line,int(UTR_start),int(UTR_end)])
									# 2nd UTR, after CDS
									CDS_found=True
									if feature == "five_prime_UTR" : feature = "three_prime_UTR"
									else : feature = "five_prime_UTR"
									utr_count = 1
									print >> sys.stderr , "-- Add " + feature + " feature form the downstream part of an exon "
									UTR_start = region_difference[1].rstrip("_").lstrip("_").split("_")[0]
									UTR_end = region_difference[1].rstrip("_").lstrip("_").split("_")[-1]
									UTR_line = "\t".join([seqname, source, feature, UTR_start, UTR_end, "." , strand, old_phase, "ID="+m_key+"."+feature+"."+str(utr_count)+";Parent="+m_key])
									#print >> sys.stderr , UTR_line
									feat_code = mRNA_feat_order[feature]
									if not feat_code in genes[g_key][3][m_key][2] : genes[g_key][3][m_key][2][feat_code] = []
									genes[g_key][3][m_key][2][feat_code].append([UTR_line,int(UTR_start),int(UTR_end)])

									prev_utr = feature

						## Calculate
						utr_count += 1

						#print >> sys.stderr, genes[g_key][3][m_key][2].keys()

					# Calculate utr lengths
					if 2 in genes[g_key][3][m_key][2] :
						# if there is any 5pUTR, check the length
						utr5_len_sum = 0
						utr5_len_list = []
						for feat in sorted(genes[g_key][3][m_key][2][2], key=lambda x: x[1]) :
							#build list
							utr5_len_list.append(feat[2]-feat[1]+1)
							utr5_len_sum += feat[2]-feat[1]+1
							#print >> sys.stderr , feat[0]

						#print >> sys.stderr , "5pUTR length from annotated regions: " + str(utr5_len_sum)
						#print >> sys.stderr , "expected legth fomr mRNA - CDS : " + str(len(UTRp5_seq))
						#print >> sys.stderr , utr5_len_list

						if utr5_len_sum == len(UTRp5_seq) :
							loc_UTRp5_region_length += utr5_len_list

						elif utr5_len_sum != 0 :
							# Correct for CDS phase
							offset = len(UTRp5_seq) - utr5_len_sum
							#print >> sys.stderr , "Offset: " + str(offset)
							if strand == "+" :
								utr5_len_list[-1]+=offset
							else :
								utr5_len_list[0]+=offset
							utr5_len_sum += offset
							loc_UTRp5_region_length+=utr5_len_list

						loc_UTRp5_exons_per_mRNA=[len(utr5_len_list)]

					if 4 in genes[g_key][3][m_key][2] :
						# if 3pUTRs are present, check the length
						utr3_len_list = []
						utr3_len_sum = 0
						for feat in sorted(genes[g_key][3][m_key][2][4], key=lambda x: x[1]) :
							#build list
							utr3_len_list.append(feat[2]-feat[1] +1 )
							utr3_len_sum += feat[2]-feat[1] +1
							#print >> sys.stderr , feat[0]

						#print >> sys.stderr , "3pUTR length from annotated regions: " + str(utr3_len_sum)
						#print >> sys.stderr , "expected legth fomr mRNA - CDS : " + str(len(UTRp3_seq))
						#print >> sys.stderr , utr3_len_list

						if utr3_len_sum == len(UTRp3_seq) :
							loc_UTRp3_region_length += utr3_len_list

						elif utr3_len_sum != 0 :
							# Correct for CDS phase
							offset = len(UTRp3_seq) - utr3_len_sum
							#print >> sys.stderr , "Offset: " + str(offset)
							if strand == "+" :
								utr3_len_list[0]+=offset
							else :
								utr3_len_list[-1]+=offset
							utr3_len_sum += offset
							loc_UTRp3_region_length += utr3_len_list

						loc_UTRp3_exons_per_mRNA=[len(utr3_len_list)]

				## mRNA level, add stats
				mRNA_region_length.append(loc_mRNA_region_length)

				exon_region_length += loc_exon_region_length
				UTRp5_region_length += loc_UTRp5_region_length
				CDS_region_length += loc_CDS_region_length
				UTRp3_region_length += loc_UTRp3_region_length
				intron_region_length += loc_intron_region_length

				mRNA_length.append(len(mRNA_seq))
				UTRp5_length += loc_UTRp5_length
				UTRp3_length += loc_UTRp3_length
				CDS_length.append(len(loc_CDS_seq.seq))
				protein_length.append(len(loc_prot_sequence))

				exons_per_mRNA.append(loc_exons_per_mRNA)
				exons_per_mRNA_local.append(loc_exons_per_mRNA)
				UTRp5_exons_per_mRNA += loc_UTRp5_exons_per_mRNA
				UTRp3_exons_per_mRNA += loc_UTRp3_exons_per_mRNA
				CDS_exons_per_mRNA.append(loc_CDS_exons_per_mRNA)

				mRNA_fasta[m_key] = mRNA_seq.seq
				CDS_fasta[m_key] = loc_CDS_seq.seq
				protein_fasta[m_key]= loc_prot_sequence


				if loc_exons_per_mRNA == 1 :
					mRNA_region_length_singleExon.append(loc_mRNA_region_length)
					exon_region_length_singleExon += loc_exon_region_length
					UTRp5_region_length_singleExon += loc_UTRp5_region_length
					CDS_region_length_singleExon += loc_CDS_region_length
					UTRp3_region_length_singleExon += loc_UTRp3_region_length

					mRNA_length_singleExon.append(len(mRNA_seq))
					UTRp5_length_singleExon += loc_UTRp5_length
					UTRp3_length_singleExon += loc_UTRp3_length
					CDS_length_singleExon.append(len(loc_CDS_seq))
					protein_length_singleExon.append(len(loc_prot_sequence))

					mRNA_singleExon_fasta[m_key] = mRNA_seq.seq
					CDS_singleExon_fasta[m_key] = loc_CDS_seq.seq
					protein_singleExon_fasta[m_key]=loc_prot_sequence

				else :
					mRNA_region_length_multiExon.append(loc_mRNA_region_length)
					exon_region_length_multiExon += loc_exon_region_length
					UTRp5_region_length_multiExon += loc_UTRp5_region_length
					CDS_region_length_multiExon += loc_CDS_region_length
					UTRp3_region_length_multiExon += loc_UTRp3_region_length

					exons_per_mRNA_multiExon.append(loc_exons_per_mRNA)
					UTRp5_exons_per_mRNA_multiExon += loc_UTRp5_exons_per_mRNA
					UTRp3_exons_per_mRNA_multiExon += loc_UTRp3_exons_per_mRNA
					CDS_exons_per_mRNA_multiExon.append(loc_CDS_exons_per_mRNA)

					mRNA_length_multiExon.append(len(mRNA_seq))
					UTRp5_length_multiExon += loc_UTRp5_length
					UTRp3_length_multiExon += loc_UTRp3_length
					CDS_length_multiExon.append(len(loc_CDS_seq))
					protein_length_multiExon.append(len(loc_prot_sequence))

					mRNA_multiExon_fasta[m_key] = mRNA_seq.seq
					CDS_multiExon_fasta[m_key] = loc_CDS_seq.seq
					protein_multiExon_fasta[m_key]= loc_prot_sequence


		## Gene level add stats

		#edit gene start and stop and calculate locus length if any mRNA is still annotated on it, otherwise drop the annotation

		if not new_gene_start == "" :
			seqname, source, feature, start, end, score, strand, frame, attribute = genes[g_key][0].rstrip().split("\t")

			genes[g_key][2] = new_gene_start
			new_gene_line = "\t".join([seqname, source, feature, str(new_gene_start), str(new_gene_end), score, strand, frame, attribute])
			genes[g_key][0] = new_gene_line

			gene_region_length.append(new_gene_end - new_gene_start +1)
			mRNA_per_gene.append(loc_mRNA_count)
			#print >> sys.stderr , exons_per_mRNA_local

			if exons_per_mRNA_local.count(1) == len(exons_per_mRNA_local) :
				gene_region_length_singleExon.append(new_gene_end - new_gene_start + 1)
				mRNA_per_gene_singleExon.append(loc_mRNA_count)
			else:
				gene_region_length_multiExon.append(new_gene_end - new_gene_start + 1)
				mRNA_per_gene_multiExon.append(loc_mRNA_count)

			# Add intergenic region feature, if appliable
			if genes[g_key][1] == prev_chr :
				print >> sys.stderr , "---- Add intergenic region"
				#add Intergenic feature to genes{}
				intergenic_name = "intergenic_" + str(intergenic_id)
				int_line = "\t".join([seqname, "GFF_feature", "intergenic", str(int(intergenic_start)+1), str(int(new_gene_start)-1), ".", ".", ".", "ID="+intergenic_name])
				#print >> sys.stderr , "-------- " + int_line
				genes[intergenic_name] = [int_line, seqname, str(int(intergenic_start)+1), {} ]
				intergenic_region_length.append(int(new_gene_start)-int(intergenic_start)-2)
				intergenic_id += 1

			prev_chr = seqname
			intergenic_start = int(new_gene_end)

		else :
			print >> sys.stderr , "------- " + g_key + " locus dropped, no coding mRNA on it "
			print >> Dropped_gene_file, g_key
			del genes[g_key]



	#print >> sys.stderr , UTRp5_exons_per_mRNA


	error_exon_file.close()
	Dropped_gene_file.close()
	Dropped_mRNA_file.close()

	if not no_cds_db == {} :
		no_cds_file = options.prefix+".no_cds.gff3"
		write_gff3( no_cds_db , no_cds_file , reference_len )

	if options.pseudo :
		Dropped_models_file.close()
		write_fasta_from_db(pseudo_fasta, Dropped_models_fasta , False)



	list_feat_length = [ gene_region_length , mRNA_region_length , exon_region_length , UTRp5_region_length , UTRp3_region_length , CDS_region_length , intron_region_length , intergenic_region_length , mRNA_length , UTRp5_length , UTRp3_length , CDS_length , protein_length , gene_region_length_singleExon , mRNA_region_length_singleExon , exon_region_length_singleExon , UTRp5_region_length_singleExon , UTRp3_region_length_singleExon , CDS_region_length_singleExon , mRNA_length_singleExon , UTRp5_length_singleExon , UTRp3_length_singleExon , protein_length_singleExon , CDS_length_singleExon , gene_region_length_multiExon , mRNA_region_length_multiExon , exon_region_length_multiExon , UTRp5_region_length_multiExon , UTRp3_region_length_multiExon , CDS_region_length_multiExon , mRNA_length_multiExon , UTRp5_length_multiExon , UTRp3_length_multiExon , protein_length_multiExon , CDS_length_multiExon]
	list_nested_subfeat_count = [ mRNA_per_gene , exons_per_mRNA , UTRp5_exons_per_mRNA , UTRp3_exons_per_mRNA , CDS_exons_per_mRNA , mRNA_per_gene_singleExon , mRNA_per_gene_multiExon , exons_per_mRNA_multiExon , CDS_exons_per_mRNA_multiExon , UTRp5_exons_per_mRNA_multiExon , UTRp3_exons_per_mRNA_multiExon , UTRp3_exons_per_mRNA_multiExon ]


	### Print GFF3
	if options.print_gff3 :
		write_gff3( genes, options.prefix+".gff3" , reference_len)

	### Print lengths lists
	if options.print_len :
		filename_gene_region_length   = open(options.prefix+".locus.lengths.txt",'w')
		filename_mRNA_region_length   = open(options.prefix+".mRNA_region.lengths.txt",'w')
		filename_exon_region_length   = open(options.prefix+".exon_region.lengths.txt",'w')
		filename_UTRp5_region_length  = open(options.prefix+".5pUTR_exon.lengths.txt",'w')
		filename_UTRp3_region_length  = open(options.prefix+".3pUTR_exon.lengths.txt",'w')
		filename_CDS_region_length    = open(options.prefix+".CDS_exon.lengths.txt",'w')
		filename_intron_region_length = open(options.prefix+".intron.lengths.txt",'w')
		filename_intergenic_region_length    = open(options.prefix+".intergenic.lengths.txt",'w')

		filename_mRNA_length          = open(options.prefix+".mRNA_sequence.lengths.txt",'w')
		filename_UTRp5_length         = open(options.prefix+".5pUTR_sequence.lengths.txt",'w')
		filename_UTRp3_length         = open(options.prefix+".3pUTR_sequence.lengths.txt",'w')
		filename_protein_length       = open(options.prefix+".protein_sequence.lengths.txt",'w')
		filename_CDS_length           = open(options.prefix+".CDS_sequence.lengths.txt",'w')


		print >> filename_gene_region_length, "\n".join(str(a) for a in gene_region_length)
		print >> filename_mRNA_region_length, "\n".join(str(a) for a in mRNA_region_length)
		print >> filename_exon_region_length, "\n".join(str(a) for a in exon_region_length)
		print >> filename_UTRp5_region_length, "\n".join(str(a) for a in UTRp5_region_length)
		print >> filename_UTRp3_region_length, "\n".join(str(a) for a in UTRp3_region_length)
		print >> filename_CDS_region_length, "\n".join(str(a) for a in CDS_region_length)
		print >> filename_intron_region_length, "\n".join(str(a) for a in intron_region_length)
		print >> filename_intergenic_region_length, "\n".join(str(a) for a in intergenic_region_length)

		print >> filename_mRNA_length, "\n".join(str(a) for a in mRNA_length)
		print >> filename_UTRp5_length, "\n".join(str(a) for a in UTRp5_length)
		print >> filename_UTRp3_length, "\n".join(str(a) for a in UTRp3_length)
		print >> filename_protein_length, "\n".join(str(a) for a in CDS_length)
		print >> filename_CDS_length, "\n".join(str(a) for a in protein_length)

		filename_gene_region_length.close()
		filename_mRNA_region_length.close()
		filename_exon_region_length.close()
		filename_UTRp5_region_length.close()
		filename_UTRp3_region_length.close()
		filename_CDS_region_length.close()
		filename_intron_region_length.close()
		filename_intergenic_region_length.close()
		filename_mRNA_length.close()
		filename_UTRp5_length.close()
		filename_UTRp3_length.close()
		filename_protein_length.close()
		filename_CDS_length.close()

		if options.print_mems :

			filename_gene_region_length_single   = open(options.prefix+".locus.monoexonic.lengths.txt",'w')
			filename_mRNA_region_length_single   = open(options.prefix+".mRNA_region.monoexonic.lengths.txt",'w')
			filename_exon_region_length_single   = open(options.prefix+".exon_region.monoexonic.lengths.txt",'w')
			filename_UTRp5_region_length_single  = open(options.prefix+".5pUTR_exon.monoexonic.lengths.txt",'w')
			filename_UTRp3_region_length_single  = open(options.prefix+".3pUTR_exon.monoexonic.lengths.txt",'w')
			filename_CDS_region_length_single    = open(options.prefix+".CDS_exon.monoexonic.lengths.txt",'w')

			filename_mRNA_length_single          = open(options.prefix+".mRNA_sequence.monoexonic.lengths.txt",'w')
			filename_UTRp5_length_single         = open(options.prefix+".5pUTR_sequence.monoexonic.lengths.txt",'w')
			filename_UTRp3_length_single         = open(options.prefix+".3pUTR_sequence.monoexonic.lengths.txt",'w')
			filename_protein_length_single       = open(options.prefix+".protein_sequence.monoexonic.lengths.txt",'w')
			filename_CDS_length_single           = open(options.prefix+".CDS_sequence.monoexonic.lengths.txt",'w')

			print >> filename_gene_region_length_single , "\n".join(str(a) for a in gene_region_length_singleExon)
			print >> filename_mRNA_region_length_single , "\n".join(str(a) for a in mRNA_region_length_singleExon)
			print >> filename_exon_region_length_single , "\n".join(str(a) for a in exon_region_length_singleExon)
			print >> filename_UTRp5_region_length_single , "\n".join(str(a) for a in UTRp5_region_length_singleExon)
			print >> filename_UTRp3_region_length_single , "\n".join(str(a) for a in UTRp3_region_length_singleExon)
			print >> filename_CDS_region_length_single , "\n".join(str(a) for a in CDS_region_length_singleExon)

			print >> filename_mRNA_length_single , "\n".join(str(a) for a in mRNA_length_singleExon)
			print >> filename_UTRp5_length_single , "\n".join(str(a) for a in UTRp5_length_singleExon)
			print >> filename_UTRp3_length_single , "\n".join(str(a) for a in UTRp3_length_singleExon)
			print >> filename_protein_length_single , "\n".join(str(a) for a in protein_length_singleExon)
			print >> filename_CDS_length_single , "\n".join(str(a) for a in CDS_length_singleExon)

			filename_gene_region_length_single.close()
			filename_mRNA_region_length_single.close()
			filename_exon_region_length_single.close()
			filename_UTRp5_region_length_single.close()
			filename_UTRp3_region_length_single.close()
			filename_CDS_region_length_single.close()

			filename_mRNA_length_single.close()
			filename_UTRp5_length_single.close()
			filename_UTRp3_length_single.close()
			filename_protein_length_single.close()
			filename_CDS_length_single.close()



			filename_gene_region_length_multi   = open(options.prefix+".locus.multiexonic.lengths.txt",'w')
			filename_mRNA_region_length_multi   = open(options.prefix+".mRNA_region.multiexonic.lengths.txt",'w')
			filename_exon_region_length_multi   = open(options.prefix+".exon_region.multiexonic.lengths.txt",'w')
			filename_UTRp5_region_length_multi  = open(options.prefix+".5pUTR_exon.multiexonic.lengths.txt",'w')
			filename_UTRp3_region_length_multi  = open(options.prefix+".3pUTR_exon.multiexonic.lengths.txt",'w')
			filename_CDS_region_length_multi    = open(options.prefix+".CDS_exon.multiexonic.lengths.txt",'w')

			filename_mRNA_length_multi          = open(options.prefix+".mRNA_sequence.multiexonic.lengths.txt",'w')
			filename_UTRp5_length_multi         = open(options.prefix+".5pUTR_sequence.multiexonic.lengths.txt",'w')
			filename_UTRp3_length_multi         = open(options.prefix+".3pUTR_sequence.multiexonic.lengths.txt",'w')
			filename_protein_length_multi       = open(options.prefix+".protein_sequence.multiexonic.lengths.txt",'w')
			filename_CDS_length_multi           = open(options.prefix+".CDS_sequence.multiexonic.lengths.txt",'w')

			print >> filename_gene_region_length_multi , "\n".join(str(a) for a in gene_region_length_multiExon)
			print >> filename_mRNA_region_length_multi , "\n".join(str(a) for a in mRNA_region_length_multiExon)
			print >> filename_exon_region_length_multi , "\n".join(str(a) for a in exon_region_length_multiExon)
			print >> filename_UTRp5_region_length_multi , "\n".join(str(a) for a in UTRp5_region_length_multiExon)
			print >> filename_UTRp3_region_length_multi , "\n".join(str(a) for a in UTRp3_region_length_multiExon)
			print >> filename_CDS_region_length_multi , "\n".join(str(a) for a in CDS_region_length_multiExon)

			print >> filename_mRNA_length_multi , "\n".join(str(a) for a in mRNA_length_multiExon)
			print >> filename_UTRp5_length_multi , "\n".join(str(a) for a in UTRp5_length_multiExon)
			print >> filename_UTRp3_length_multi , "\n".join(str(a) for a in UTRp3_length_multiExon)
			print >> filename_protein_length_multi , "\n".join(str(a) for a in protein_length_multiExon)
			print >> filename_CDS_length_multi , "\n".join(str(a) for a in CDS_length_multiExon)

			filename_gene_region_length_multi.close()
			filename_mRNA_region_length_multi.close()
			filename_exon_region_length_multi.close()
			filename_UTRp5_region_length_multi.close()
			filename_UTRp3_region_length_multi.close()
			filename_CDS_region_length_multi.close()

			filename_mRNA_length_multi.close()
			filename_UTRp5_length_multi.close()
			filename_UTRp3_length_multi.close()
			filename_protein_length_multi.close()
			filename_CDS_length_multi.close()


	### Print counts lists
	if options.print_counts:
		filename_mRNA_per_gene        = open(options.prefix+".mRNA_per_locus.counts.txt",'w')
		filename_exons_per_mRNA       = open(options.prefix+".exon_per_mRNA.counts.txt",'w')
		filename_UTRp5_exons_per_mRNA = open(options.prefix+".5pUTR_per_mRNA.counts.txt",'w')
		filename_UTRp3_exons_per_mRNA = open(options.prefix+".3pUTR_per_mRNA.counts.txt",'w')
		filename_CDS_exons_per_mRNA   = open(options.prefix+".CDS_per_mRNA.counts.txt",'w')

		print >> filename_mRNA_per_gene, "\n".join(str(a) for a in mRNA_per_gene)
		print >> filename_exons_per_mRNA, "\n".join(str(a) for a in exons_per_mRNA)
		print >> filename_UTRp5_exons_per_mRNA, "\n".join(str(a) for a in UTRp5_exons_per_mRNA)
		print >> filename_UTRp3_exons_per_mRNA, "\n".join(str(a) for a in UTRp3_exons_per_mRNA)
		print >> filename_CDS_exons_per_mRNA, "\n".join(str(a) for a in CDS_exons_per_mRNA)

		filename_mRNA_per_gene.close()
		filename_exons_per_mRNA.close()
		filename_UTRp5_exons_per_mRNA.close()
		filename_UTRp3_exons_per_mRNA.close()
		filename_CDS_exons_per_mRNA.close()

		if options.print_mems :
			filename_mRNA_per_gene_single = open(options.prefix+".mRNA_per_locus.monoexonic.counts.txt",'w')
			filename_mRNA_per_gene_multi  = open(options.prefix+".mRNA_per_locus.multiexonic.counts.txt",'w')
			filename_exons_per_mRNA_multi = open(options.prefix+".exon_per_mRNA.multiexonic.counts.txt",'w')
			filename_UTRp5_exons_per_mRNA_multi = open(options.prefix+".5pUTR_per_mRNA.multiexonic.counts.txt",'w')
			filename_UTRp3_exons_per_mRNA_multi = open(options.prefix+".3pUTR_per_mRNA.multiexonic.counts.txt",'w')
			filename_CDS_exons_per_mRNA_multi= open(options.prefix+".CDS_per_mRNA.multiexonic.counts.txt",'w')

			print >> filename_mRNA_per_gene_single, "\n".join(str(a) for a in mRNA_per_gene_singleExon)
			print >> filename_mRNA_per_gene_multi, "\n".join(str(a) for a in mRNA_per_gene_multiExon)
			print >> filename_exons_per_mRNA_multi, "\n".join(str(a) for a in exons_per_mRNA_multiExon)
			print >> filename_UTRp5_exons_per_mRNA_multi, "\n".join(str(a) for a in UTRp5_exons_per_mRNA_multiExon)
			print >> filename_UTRp3_exons_per_mRNA_multi, "\n".join(str(a) for a in UTRp3_exons_per_mRNA_multiExon)
			print >> filename_CDS_exons_per_mRNA_multi, "\n".join(str(a) for a in CDS_exons_per_mRNA_multiExon)

			filename_mRNA_per_gene_single.close()
			filename_mRNA_per_gene_multi.close()
			filename_exons_per_mRNA_multi.close()
			filename_UTRp5_exons_per_mRNA_multi.close()
			filename_UTRp3_exons_per_mRNA_multi.close()
			filename_CDS_exons_per_mRNA_multi.close()

	### Print stats
	if options.print_stats:
		### Stats comupted for each feature
		# Count
		# Global Length
		# Average per Locus
		# Avereage per mRNA
		# Max Length
		# Average Length
		# Length St. Dev.
		# Length 25% quantile
		# Length 50% quantile / Median
		# Length 75% quantile

		colnames = [
					"Count",
					"Cumulative Length",
					"Average per locus",
					"Average per mRNA",
					"Max",
					"Average",
					"St. Dev.",
					"25th quantile",
					"50th quantile / Median",
					"75th quantile"
					]

		formatting_db = {
					'Count': '{:,.0f}',
					'Cumulative Length': '{:,.0f}',
					'Average per locus': '{:,.2f}',
					'Average per mRNA': '{:,.2f}',
					'Max': '{:,.0f}',
					'Average': '{:,.2f}',
					'St. Dev.': '{:,.2f}',
					'25th quantile': '{:,.0f}',
					'50th quantile / Median': '{:,.0f}',
					'75th quantile': '{:,.0f}',
				}

		formatter_obj = {k: v.format for k, v in formatting_db.items()}

		filename_stats = open(options.prefix+".stats.txt", 'w')

		print >> filename_stats, "#########################################################"
		print >> filename_stats, "#### Annotation statistics ##############################"
		print >> filename_stats, "Annotation file\t"+options.gff
		print >> filename_stats, "Genome file\t"+options.genome

		genome_seq_lenths = sorted([ len(reference[x].seq) for x in reference ] ,reverse=True)

		#print >> sys.stderr, genome_seq_lenths

		print >> filename_stats , "\n\n"
		print >> filename_stats , "#########################################################"
		print >> filename_stats , "#### Genome Statistics ##################################"
		print >> filename_stats , "Genome length:" + "\t" + str('%.0f' % sum(genome_seq_lenths))
		print >> filename_stats , "Number of sequences:" + "\t" + str('%.0f' % len(genome_seq_lenths))
		print >> filename_stats , "Average sequence length:" + "\t" + str('%.2f' % np.mean(genome_seq_lenths))
		print >> filename_stats , "Median_seqence_length" + "\t" + str('%.0f' % np.median(genome_seq_lenths))
		print >> filename_stats , "Minimum sequence length:" + "\t" + str('%.0f' % np.min(genome_seq_lenths))
		print >> filename_stats , "Maximum sequence length:" + "\t" + str('%.0f' % np.max(genome_seq_lenths))

		N , L = Nvalue(25 ,genome_seq_lenths)
		print >> filename_stats , "N25 length:\t" + str(N) + "\tIndex:\t" + str(L)
		N , L = Nvalue(50 ,genome_seq_lenths)
		print >> filename_stats , "N50 length:\t" + str(N) + "\tIndex:\t" + str(L)
		N , L = Nvalue(75 ,genome_seq_lenths)
		print >> filename_stats , "N75 length:\t" + str(N) + "\tIndex:\t" + str(L)
		N , L = Nvalue(90 ,genome_seq_lenths)
		print >> filename_stats , "N90 length:\t" + str(N) + "\tIndex:\t" + str(L)

		print >> filename_stats , "Sequences > 100b - Count:\t"  +  str('%.0f' % len([x for x in genome_seq_lenths if x > 100 ])) + "\tCumulative length:\t" + str('%.0f' % sum([x for x in genome_seq_lenths if x > 100 ]))
		print >> filename_stats , "Sequences > 500b - Count:\t"  +  str('%.0f' % len([x for x in genome_seq_lenths if x > 500 ])) + "\tCumulative length:\t" + str('%.0f' % sum([x for x in genome_seq_lenths if x > 500 ]))
		print >> filename_stats , "Sequences > 1Kb - Count:\t"  +  str('%.0f' % len([x for x in genome_seq_lenths if x > 1000 ])) + "\tCumulative length:\t" + str('%.0f' % sum([x for x in genome_seq_lenths if x > 1000 ]))
		print >> filename_stats , "Sequences > 5Kb - Count:\t"  +  str('%.0f' % len([x for x in genome_seq_lenths if x > 5000 ])) + "\tCumulative length:\t" + str('%.0f' % sum([x for x in genome_seq_lenths if x > 5000 ]))
		print >> filename_stats , "Sequences > 10Kb - Count:\t"  +  str('%.0f' % len([x for x in genome_seq_lenths if x > 10000 ])) + "\tCumulative length:\t" + str('%.0f' % sum([x for x in genome_seq_lenths if x > 10000 ]))
		print >> filename_stats , "Sequences > 50Kb - Count:\t"  +  str('%.0f' % len([x for x in genome_seq_lenths if x > 50000 ])) + "\tCumulative length:\t" + str('%.0f' % sum([x for x in genome_seq_lenths if x > 50000 ]))
		print >> filename_stats , "Sequences > 100Kb - Count:\t"  +  str('%.0f' % len([x for x in genome_seq_lenths if x > 100000 ])) + "\tCumulative length:\t" + str('%.0f' % sum([x for x in genome_seq_lenths if x > 100000 ]))
		print >> filename_stats , "Sequences > 500Kb - Count:\t"  +  str('%.0f' % len([x for x in genome_seq_lenths if x > 500000 ])) + "\tCumulative length:\t" + str('%.0f' % sum([x for x in genome_seq_lenths if x > 500000 ]))
		print >> filename_stats , "Sequences > 1Mb - Count:\t"  +  str('%.0f' % len([x for x in genome_seq_lenths if x > 1000000 ])) + "\tCumulative length:\t" + str('%.0f' % sum([x for x in genome_seq_lenths if x > 1000000 ]))
		print >> filename_stats , "Sequences > 5Mb - Count:\t"  +  str('%.0f' % len([x for x in genome_seq_lenths if x > 5000000 ])) + "\tCumulative length:\t" + str('%.0f' % sum([x for x in genome_seq_lenths if x > 5000000 ]))

		print >> filename_stats, "\n\n"
		print >> filename_stats, "#########################################################"
		print >> filename_stats, "#### Overall ############################################"

		num_loci = len(gene_region_length)
		num_mRNAs = len(mRNA_length)

		overall = np.array(feat_len_stats(gene_region_length,num_loci,0))
		overall = np.vstack([overall,feat_len_stats(mRNA_region_length,num_loci,num_mRNAs)])
		overall = np.vstack([overall,count_stats(mRNA_per_gene)])
		overall = np.vstack([overall,feat_len_stats(intergenic_region_length,0,0)])
		overall = np.vstack([overall,feat_len_stats(exon_region_length,num_loci,num_mRNAs)])
		overall = np.vstack([overall,count_stats(exons_per_mRNA)])
		overall = np.vstack([overall,feat_len_stats(intron_region_length,num_loci,num_mRNAs)])
		overall = np.vstack([overall,feat_len_stats(UTRp5_region_length,num_loci,num_mRNAs)])
		overall = np.vstack([overall,count_stats(UTRp5_exons_per_mRNA)])
		overall = np.vstack([overall,feat_len_stats(UTRp3_region_length,num_loci,num_mRNAs)])
		overall = np.vstack([overall,count_stats(UTRp3_exons_per_mRNA)])
		overall = np.vstack([overall,feat_len_stats(CDS_region_length,num_loci,num_mRNAs)])
		overall = np.vstack([overall,count_stats(CDS_exons_per_mRNA)])
		overall = np.vstack([overall,feat_len_stats(mRNA_length,num_loci,num_mRNAs)])
		overall = np.vstack([overall,feat_len_stats(CDS_length,num_loci,num_mRNAs)])
		overall = np.vstack([overall,feat_len_stats(UTRp5_length,num_loci,num_mRNAs)])
		overall = np.vstack([overall,feat_len_stats(UTRp3_length,num_loci,num_mRNAs)])
		overall = np.vstack([overall,feat_len_stats(protein_length,num_loci,num_mRNAs)])

		rownames_overall = [
					"Loci",
					"mRNAs",
					"mRNAs count/locus",
					"Intergenic regions lengths",
					"Exons lengths",
					"Exons count/mRNA",
					"Introns lengths",
					"5' UTRs exons lengths",
					"5' UTRs exons count/mRNA",
					"3' UTRs exons lengths",
					"3' UTRs exons count/mRNA",
					"CDSs exons lengths",
					"CDSs exons count/mRNA",
					"mRNAs sequences lengths",
					"CDSs sequences lengths",
					"5' UTRs sequences lengths",
					"3' UTRs sequences lengths",
					"Protein sequences lengths"
					]

		print >> filename_stats, pd.DataFrame(overall,columns=colnames, index=rownames_overall).to_string(formatters=formatter_obj)

		if options.print_mems :

			print >> filename_stats, "\n\n"
			print >> filename_stats, "#########################################################"
			print >> filename_stats, "#### Mono-exonic ########################################"

			num_loci = len(gene_region_length_singleExon)
			num_mRNAs = len(mRNA_region_length_singleExon)

			mono = np.array(feat_len_stats(gene_region_length_singleExon,num_loci,0))
			mono = np.vstack([mono,feat_len_stats(mRNA_region_length_singleExon,num_loci,num_mRNAs)])
			mono = np.vstack([mono,count_stats(mRNA_per_gene_singleExon)])
			mono = np.vstack([mono,feat_len_stats(exon_region_length_singleExon,num_loci,num_mRNAs)])
			mono = np.vstack([mono,feat_len_stats(UTRp5_region_length_singleExon,num_loci,num_mRNAs)])
			mono = np.vstack([mono,feat_len_stats(UTRp3_region_length_singleExon,num_loci,num_mRNAs)])
			mono = np.vstack([mono,feat_len_stats(CDS_region_length_singleExon,num_loci,num_mRNAs)])
			mono = np.vstack([mono,feat_len_stats(mRNA_length_singleExon,num_loci,num_mRNAs)])
			mono = np.vstack([mono,feat_len_stats(CDS_length_singleExon,num_loci,num_mRNAs)])
			mono = np.vstack([mono,feat_len_stats(UTRp5_length_singleExon,num_loci,num_mRNAs)])
			mono = np.vstack([mono,feat_len_stats(UTRp3_length_singleExon,num_loci,num_mRNAs)])
			mono = np.vstack([mono,feat_len_stats(protein_length_singleExon,num_loci,num_mRNAs)])

			rownames_mono = [
						"Loci lengths",
						"mRNAs lengths",
						"mRNAs count/gene",
						"Exons lengths",
						"5' UTRs exons lengths",
						"3' UTRs exons lengths",
						"CDSs exons lengths",
						"mRNAs sequences lengths",
						"CDSs sequences lengths",
						"5' UTRs sequences lengths",
						"3' UTRs sequences lengths",
						"Protein sequences lengths"
						]

			print >> filename_stats, pd.DataFrame(mono,columns=colnames, index=rownames_mono).to_string(formatters=formatter_obj)

			print >> filename_stats, "\n\n"
			print >> filename_stats, "#########################################################"
			print >> filename_stats, "#### Multi-exonic #######################################"

			num_loci = len(gene_region_length_multiExon)
			num_mRNAs = len(mRNA_region_length_multiExon)

			multi = np.array(feat_len_stats(gene_region_length_multiExon,num_loci,0))
			multi = np.vstack([multi,feat_len_stats(mRNA_region_length_multiExon,num_loci,num_mRNAs)])
			multi = np.vstack([multi,count_stats(mRNA_per_gene_multiExon)])
			multi = np.vstack([multi,feat_len_stats(exon_region_length_multiExon,num_loci,num_mRNAs)])
			multi = np.vstack([multi,count_stats(exons_per_mRNA_multiExon)])
			multi = np.vstack([multi,feat_len_stats(UTRp5_region_length_multiExon,num_loci,num_mRNAs)])
			multi = np.vstack([multi,count_stats(UTRp5_exons_per_mRNA_multiExon)])
			multi = np.vstack([multi,feat_len_stats(UTRp3_region_length_multiExon,num_loci,num_mRNAs)])
			multi = np.vstack([multi,count_stats(UTRp3_exons_per_mRNA_multiExon)])
			multi = np.vstack([multi,feat_len_stats(CDS_region_length_multiExon,num_loci,num_mRNAs)])
			multi = np.vstack([multi,count_stats(CDS_exons_per_mRNA_multiExon)])
			multi = np.vstack([multi,feat_len_stats(mRNA_length_multiExon,num_loci,num_mRNAs)])
			multi = np.vstack([multi,feat_len_stats(CDS_length_multiExon,num_loci,num_mRNAs)])
			multi = np.vstack([multi,feat_len_stats(UTRp5_length_multiExon,num_loci,num_mRNAs)])
			multi = np.vstack([multi,feat_len_stats(UTRp3_length_multiExon,num_loci,num_mRNAs)])
			multi = np.vstack([multi,feat_len_stats(protein_length_multiExon,num_loci,num_mRNAs)])

			rownames_multi = [
						"Loci lengths",
						"mRNAs lengths",
						"mRNAs count/gene",
						"Exons lengths",
						"Exons count/mRNA",
						"5' UTRs exons lengths",
						"5' UTRs exons count/mRNA",
						"3' UTRs exons lengths",
						"3' UTRs exons count/mRNA",
						"CDSs exons lengths",
						"CDSs exons count/mRNA",
						"mRNAs sequences lengths",
						"CDSs sequences lengths",
						"5' UTRs sequences lengths",
						"3' UTRs sequences lengths",
						"Protein sequences lengths"
						]

			print >> filename_stats, pd.DataFrame(multi,columns=colnames, index=rownames_multi).to_string(formatters=formatter_obj)


	if options.plot_distributions :
		plot_distributions( reference , list_feat_length , list_nested_subfeat_count  , options.prefix , options.print_mems )


	### Print fasta files
	if options.print_seq :

		filename_mRNA_fasta           = options.prefix+".mRNA.fasta"
		filename_CDS_fasta            = options.prefix+".CDS.fasta"
		filename_protein_fasta        = options.prefix+".protein.fasta"

		filename_mRNA_fasta_single    = options.prefix+".mRNA.monoexonic.fasta"
		filename_CDS_fasta_single     = options.prefix+".CDS.monoexonic.fasta"
		filename_protein_fasta_single = options.prefix+".protein.monoexonic.fasta"

		filename_mRNA_fasta_multi     = options.prefix+".mRNA.multiexonic.fasta"
		filename_CDS_fasta_multi      = options.prefix+".CDS.multiexonic.fasta"
		filename_protein_fasta_multi  = options.prefix+".protein.multiexonic.fasta"


		if options.print_mems :

			process_list = [ (filename_mRNA_fasta, mRNA_fasta ),
				(filename_CDS_fasta, CDS_fasta),
				(filename_protein_fasta, protein_fasta),
				(filename_mRNA_fasta_single, mRNA_singleExon_fasta),
				(filename_CDS_fasta_single, CDS_singleExon_fasta),
				(filename_protein_fasta_single, protein_singleExon_fasta),
				(filename_mRNA_fasta_multi , mRNA_multiExon_fasta),
				(filename_CDS_fasta_multi , CDS_multiExon_fasta),
				(filename_protein_fasta_multi , protein_multiExon_fasta) ]

		else :
			process_list = [ (filename_mRNA_fasta, mRNA_fasta ), (filename_CDS_fasta, CDS_fasta), (filename_protein_fasta, protein_fasta) ]

		for process in process_list :
			write_fasta_from_db( process[1] , process[0] , False)

	print >> sys.stderr, "-------------------------"
	print >> sys.stderr, "----    Completed    ----"
	print >> sys.stderr, "-------------------------"

	print >> sys.stdout, "-------------------------"
	print >> sys.stdout, "----    Completed    ----"
	print >> sys.stdout, "-------------------------"

if __name__ == '__main__':
	main()
