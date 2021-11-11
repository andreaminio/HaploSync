#!/usr/bin/env python

import sys
import string
import random
from operator import itemgetter
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from optparse import OptionParser
import textwrap
from multiprocessing import Pool
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.font_manager as font_manager
import gc
import os
from collections import defaultdict
from HaploFunct import *
from AGP_lib import *
from FASTA_lib import *

#### Functions

#### genes db structure:
#### genes[gene_id] = { gff_line , chr , int(start) , mRNA_dict }
####		mRNA_dict[mRNA_id] = { gff_line , int(start), feat_dict{} }
####			feat_dict[mRNA_order_num] = [ [ line , int(start) , int(end) ]  ]
####		 		with mRNA_order_num meaning:	1 = "exon"
####												2 = "five_prime_UTR"
####												3 = "CDS"
####												4 = "three_prime_UTR"


def feat_len_stats(list , n_loci , n_mRNAs ):

	if list == [] : stats_array = [0,0,0,0,0,0,0,0,0,0]
	else :
		stats_array = []

		stats_array.append(len(list)) # Counts
		stats_array.append(sum(list)) # Global length
		if n_loci >0 :
			stats_array.append(float(len(list))/float(n_loci)) # Average per Locus
		else :
			stats_array.append(0)
		if n_mRNAs>0 :
			stats_array.append(float(len(list))/float(n_mRNAs)) # Avereage per mRNA
		else :
			stats_array.append(0)
		stats_array.append(np.max(list)) # Max Length
		stats_array.append(np.mean(list)) # Average Length
		stats_array.append(np.std(list)) # Length St. Dev.
		stats_array += np.percentile(list,[25,50,75]).tolist()# Length 25% quantile , Length 50% quantile / Median,  Length 75% quantile

	return(stats_array)


def count_stats(list):

	if list == [] : stats_array = [0,0,0,0,0,0,0,0,0,0]
	else :
		stats_array = []

		stats_array.append(0) # Counts
		stats_array.append(0) # Global length
		stats_array.append(0) # Average per Locus
		stats_array.append(0) # Avereage per mRNA
		stats_array.append(np.max(list)) # Max Length
		stats_array.append(np.mean(list)) # Average Length
		stats_array.append(np.std(list)) # Length St. Dev.
		stats_array += np.percentile(list,[25,50,75]).tolist()# Length 25% quantile , Length 50% quantile / Median,  Length 75% quantile

	return(stats_array)


def read_gff3(gff3_file) :
	mRNA_feat_order = {"exon": 1, "intron":1, "five_prime_UTR":2 , "CDS":3, "three_prime_UTR":4}
	genes = {}
	mRNA = {}

	## Read original gff
	for line in open(gff3_file):
		#print >> sys.stdout, "#\t" + line.rstrip()
		if line[0] == "#" or line.rstrip() == "" :
			continue  # print >> sys.stderr, line.rstrip()
		else:
			attributes_dict = {}
			seqname, source, feature, start, end, score, strand, frame, attribute = line.rstrip().split("\t")
			for chunk in attribute.rstrip(";").split(";"):
				att = chunk.split("=")
				try : attributes_dict[att[0]] = att[1]
				except : print >> sys.stderr , "error in attributes at line: " + line

			if feature == "gene":
				try:
					gene_name = attributes_dict["ID"]
				except :
					print >> sys.stderr, "[ERROR] Gene feature missing ID @ line: " + line
					sys.exit(1)
				if gene_name in genes :
					genes[gene_name][0] = line.rstrip()
					genes[gene_name][1] = seqname
					genes[gene_name][2] = int(start)
				else:
					genes[gene_name] = [line.rstrip(), seqname, int(start), {}]

			elif feature == "mRNA":
				try :
					gene_name = attributes_dict["Parent"]
				except :
					print >> sys.stderr, "[ERROR] mRNA missing parent information @ line: " + line
					sys.exit(1)

				try:
					mRNA_name = attributes_dict["ID"]
				except :
					print >> sys.stderr, "[ERROR] mRNA feature missing ID @ line: " + line
					sys.exit(1)

				if not attributes_dict["ID"] in mRNA :
					mRNA[attributes_dict["ID"]] = gene_name
				else :
					print >> sys.stderr, "[ERROR] Duplicated mRNA ID: " + attributes_dict["ID"]
					sys.exit(1)

				if gene_name in genes :
					genes[gene_name][3][attributes_dict["ID"]] = [line.rstrip(), int(start), {}]
				else :
					genes[gene_name] = ["", "", "", {}]
					genes[gene_name][3][attributes_dict["ID"]] = [line.rstrip(), int(start), {}]

			elif feature == "intergenic" :
				continue

			elif feature == "intron" :
				continue

			else :
				try :
					mRNA_ID = attributes_dict["Parent"]
				except :
					print >> sys.stderr, "[ERROR] Feature missing parent information @ line: " + line
					sys.exit(1)

				if len(mRNA_ID.split(","))>1 :
					for single_id in mRNA_ID.split(",") :
						if single_id in mRNA :
							gene_name = mRNA[single_id]
							attributes_dict["Parent"] = single_id
							attribute = ""
							for key in reversed(attributes_dict.keys()) : attribute = attribute + key + "=" + attributes_dict[key] + ";"
							if mRNA_feat_order[feature] not in genes[gene_name][3][single_id][2] :
								genes[gene_name][3][single_id][2][mRNA_feat_order[feature]]=[]
							new_line = "\t".join([seqname, source, feature, start, end, score, strand, frame, attribute])
							genes[gene_name][3][single_id][2][mRNA_feat_order[feature]].append([new_line, int(start),int(end)])
				else:
					if mRNA_ID in mRNA :
						gene_name = mRNA[mRNA_ID]
						if mRNA_feat_order[feature] not in genes[gene_name][3][mRNA_ID][2] :
							genes[gene_name][3][mRNA_ID][2][mRNA_feat_order[feature]]=[]
						genes[gene_name][3][mRNA_ID][2][mRNA_feat_order[feature]].append([line.rstrip(), int(start),int(end)])
					else :
						print >> sys.stderr, "[ERROR] Feature declared before parent mRNA @ line: " + line
						sys.exit(1)
	return genes , mRNA


def write_gff3( genes , out_name , reference_seq_len ) :
	filename_gff3 = open( out_name ,"w")
	print >> filename_gff3, "##gff-version 3"
	chr_printed=""
	for g_key, g_value in sorted(genes.items(), key=lambda gene: gene[1:2]) :
		actual_chr=genes[g_key][1]
		if chr_printed != actual_chr :
			print >> filename_gff3, "##sequence-region " + str(actual_chr) + " 1 " + str( reference_seq_len[actual_chr] )
			chr_printed = actual_chr
		print >> filename_gff3, "### " + g_key
		print >> filename_gff3, genes[g_key][0]
		# Iterate on mRNAs
		for m_key in sorted(genes[g_key][3], key=lambda mRNA: mRNA[1]) :
			print >> filename_gff3, genes[g_key][3][m_key][0]
			for feat_key in sorted(genes[g_key][3][m_key][2].keys()):
				for el in sorted( genes[g_key][3][m_key][2][feat_key] , key=lambda feat : feat[1]) :
					print >> filename_gff3, el[0]


def print_len( prefix , lengths_list , subfeat_count_list , reference , print_mems=True) :
	prefix = str(prefix)
	gene_region_length , mRNA_region_length , exon_region_length , UTRp5_region_length , UTRp3_region_length , CDS_region_length , intron_region_length , intergenic_region_length , mRNA_length , UTRp5_length , UTRp3_length , CDS_length , protein_length , gene_region_length_singleExon , mRNA_region_length_singleExon , exon_region_length_singleExon , UTRp5_region_length_singleExon , UTRp3_region_length_singleExon , CDS_region_length_singleExon , mRNA_length_singleExon , UTRp5_length_singleExon , UTRp3_length_singleExon , protein_length_singleExon , CDS_length_singleExon , gene_region_length_multiExon , mRNA_region_length_multiExon , exon_region_length_multiExon , UTRp5_region_length_multiExon , UTRp3_region_length_multiExon , CDS_region_length_multiExon , mRNA_length_multiExon , UTRp5_length_multiExon , UTRp3_length_multiExon , protein_length_multiExon , CDS_length_multiExon = lengths_list
	mRNA_per_gene , exons_per_mRNA , UTRp5_exons_per_mRNA , UTRp3_exons_per_mRNA , CDS_exons_per_mRNA , mRNA_per_gene_singleExon , mRNA_per_gene_multiExon , exons_per_mRNA_multiExon , CDS_exons_per_mRNA_multiExon , UTRp5_exons_per_mRNA_multiExon , UTRp3_exons_per_mRNA_multiExon , UTRp3_exons_per_mRNA_multiExon = subfeat_count_list
	genome_seq_lenths = sorted([ len(reference[x].seq) for x in reference ] ,reverse=True)

	filename_gene_region_length   = open(prefix+".locus.lengths.txt",'w')
	filename_mRNA_region_length   = open(prefix+".mRNA_region.lengths.txt",'w')
	filename_exon_region_length   = open(prefix+".exon_region.lengths.txt",'w')
	filename_UTRp5_region_length  = open(prefix+".5pUTR_exon.lengths.txt",'w')
	filename_UTRp3_region_length  = open(prefix+".3pUTR_exon.lengths.txt",'w')
	filename_CDS_region_length    = open(prefix+".CDS_exon.lengths.txt",'w')
	filename_intron_region_length = open(prefix+".intron.lengths.txt",'w')
	filename_intergenic_region_length    = open(prefix+".intergenic.lengths.txt",'w')

	filename_mRNA_length          = open(prefix+".mRNA_sequence.lengths.txt",'w')
	filename_UTRp5_length         = open(prefix+".5pUTR_sequence.lengths.txt",'w')
	filename_UTRp3_length         = open(prefix+".3pUTR_sequence.lengths.txt",'w')
	filename_protein_length       = open(prefix+".protein_sequence.lengths.txt",'w')
	filename_CDS_length           = open(prefix+".CDS_sequence.lengths.txt",'w')


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

	if print_mems :

		filename_gene_region_length_single   = open(prefix+".locus.monoexonic.lengths.txt",'w')
		filename_mRNA_region_length_single   = open(prefix+".mRNA_region.monoexonic.lengths.txt",'w')
		filename_exon_region_length_single   = open(prefix+".exon_region.monoexonic.lengths.txt",'w')
		filename_UTRp5_region_length_single  = open(prefix+".5pUTR_exon.monoexonic.lengths.txt",'w')
		filename_UTRp3_region_length_single  = open(prefix+".3pUTR_exon.monoexonic.lengths.txt",'w')
		filename_CDS_region_length_single    = open(prefix+".CDS_exon.monoexonic.lengths.txt",'w')

		filename_mRNA_length_single          = open(prefix+".mRNA_sequence.monoexonic.lengths.txt",'w')
		filename_UTRp5_length_single         = open(prefix+".5pUTR_sequence.monoexonic.lengths.txt",'w')
		filename_UTRp3_length_single         = open(prefix+".3pUTR_sequence.monoexonic.lengths.txt",'w')
		filename_protein_length_single       = open(prefix+".protein_sequence.monoexonic.lengths.txt",'w')
		filename_CDS_length_single           = open(prefix+".CDS_sequence.monoexonic.lengths.txt",'w')

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



		filename_gene_region_length_multi   = open(prefix+".locus.multiexonic.lengths.txt",'w')
		filename_mRNA_region_length_multi   = open(prefix+".mRNA_region.multiexonic.lengths.txt",'w')
		filename_exon_region_length_multi   = open(prefix+".exon_region.multiexonic.lengths.txt",'w')
		filename_UTRp5_region_length_multi  = open(prefix+".5pUTR_exon.multiexonic.lengths.txt",'w')
		filename_UTRp3_region_length_multi  = open(prefix+".3pUTR_exon.multiexonic.lengths.txt",'w')
		filename_CDS_region_length_multi    = open(prefix+".CDS_exon.multiexonic.lengths.txt",'w')

		filename_mRNA_length_multi          = open(prefix+".mRNA_sequence.multiexonic.lengths.txt",'w')
		filename_UTRp5_length_multi         = open(prefix+".5pUTR_sequence.multiexonic.lengths.txt",'w')
		filename_UTRp3_length_multi         = open(prefix+".3pUTR_sequence.multiexonic.lengths.txt",'w')
		filename_protein_length_multi       = open(prefix+".protein_sequence.multiexonic.lengths.txt",'w')
		filename_CDS_length_multi           = open(prefix+".CDS_sequence.multiexonic.lengths.txt",'w')

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


def print_counts (prefix , lengths_list , subfeat_count_list , reference , print_mems=True):
	prefix = str(prefix)
	gene_region_length , mRNA_region_length , exon_region_length , UTRp5_region_length , UTRp3_region_length , CDS_region_length , intron_region_length , intergenic_region_length , mRNA_length , UTRp5_length , UTRp3_length , CDS_length , protein_length , gene_region_length_singleExon , mRNA_region_length_singleExon , exon_region_length_singleExon , UTRp5_region_length_singleExon , UTRp3_region_length_singleExon , CDS_region_length_singleExon , mRNA_length_singleExon , UTRp5_length_singleExon , UTRp3_length_singleExon , protein_length_singleExon , CDS_length_singleExon , gene_region_length_multiExon , mRNA_region_length_multiExon , exon_region_length_multiExon , UTRp5_region_length_multiExon , UTRp3_region_length_multiExon , CDS_region_length_multiExon , mRNA_length_multiExon , UTRp5_length_multiExon , UTRp3_length_multiExon , protein_length_multiExon , CDS_length_multiExon = lengths_list
	mRNA_per_gene , exons_per_mRNA , UTRp5_exons_per_mRNA , UTRp3_exons_per_mRNA , CDS_exons_per_mRNA , mRNA_per_gene_singleExon , mRNA_per_gene_multiExon , exons_per_mRNA_multiExon , CDS_exons_per_mRNA_multiExon , UTRp5_exons_per_mRNA_multiExon , UTRp3_exons_per_mRNA_multiExon , UTRp3_exons_per_mRNA_multiExon = subfeat_count_list
	genome_seq_lenths = sorted([ len(reference[x].seq) for x in reference ] ,reverse=True)

	filename_mRNA_per_gene        = open(prefix+".mRNA_per_locus.counts.txt",'w')
	filename_exons_per_mRNA       = open(prefix+".exon_per_mRNA.counts.txt",'w')
	filename_UTRp5_exons_per_mRNA = open(prefix+".5pUTR_per_mRNA.counts.txt",'w')
	filename_UTRp3_exons_per_mRNA = open(prefix+".3pUTR_per_mRNA.counts.txt",'w')
	filename_CDS_exons_per_mRNA   = open(prefix+".CDS_per_mRNA.counts.txt",'w')

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

	if print_mems :
		filename_mRNA_per_gene_single = open(prefix+".mRNA_per_locus.monoexonic.counts.txt",'w')
		filename_mRNA_per_gene_multi  = open(prefix+".mRNA_per_locus.multiexonic.counts.txt",'w')
		filename_exons_per_mRNA_multi = open(prefix+".exon_per_mRNA.multiexonic.counts.txt",'w')
		filename_UTRp5_exons_per_mRNA_multi = open(prefix+".5pUTR_per_mRNA.multiexonic.counts.txt",'w')
		filename_UTRp3_exons_per_mRNA_multi = open(prefix+".3pUTR_per_mRNA.multiexonic.counts.txt",'w')
		filename_CDS_exons_per_mRNA_multi= open(prefix+".CDS_per_mRNA.multiexonic.counts.txt",'w')

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


def print_stats( prefix , genome , reference , gff , lengths_list , subfeat_count_list , print_mems=True ):

	prefix = str(prefix)
	gene_region_length , mRNA_region_length , exon_region_length , UTRp5_region_length , UTRp3_region_length , CDS_region_length , intron_region_length , intergenic_region_length , mRNA_length , UTRp5_length , UTRp3_length , CDS_length , protein_length , gene_region_length_singleExon , mRNA_region_length_singleExon , exon_region_length_singleExon , UTRp5_region_length_singleExon , UTRp3_region_length_singleExon , CDS_region_length_singleExon , mRNA_length_singleExon , UTRp5_length_singleExon , UTRp3_length_singleExon , protein_length_singleExon , CDS_length_singleExon , gene_region_length_multiExon , mRNA_region_length_multiExon , exon_region_length_multiExon , UTRp5_region_length_multiExon , UTRp3_region_length_multiExon , CDS_region_length_multiExon , mRNA_length_multiExon , UTRp5_length_multiExon , UTRp3_length_multiExon , protein_length_multiExon , CDS_length_multiExon = lengths_list
	mRNA_per_gene , exons_per_mRNA , UTRp5_exons_per_mRNA , UTRp3_exons_per_mRNA , CDS_exons_per_mRNA , mRNA_per_gene_singleExon , mRNA_per_gene_multiExon , exons_per_mRNA_multiExon , CDS_exons_per_mRNA_multiExon , UTRp5_exons_per_mRNA_multiExon , UTRp3_exons_per_mRNA_multiExon , UTRp3_exons_per_mRNA_multiExon = subfeat_count_list
	genome_seq_lenths = sorted([ len(reference[x].seq) for x in reference ] ,reverse=True)


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

	filename_stats = open(prefix+".stats.txt", 'w')

	print >> filename_stats, "#########################################################"
	print >> filename_stats, "#### Annotation statistics ##############################"
	print >> filename_stats, "Annotation file\t" + gff
	print >> filename_stats, "Genome file\t" + genome

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

	overall = np.array(len_stats(gene_region_length,num_loci,0))
	overall = np.vstack([overall,len_stats(mRNA_region_length,num_loci,num_mRNAs)])
	overall = np.vstack([overall,count_stats(mRNA_per_gene)])
	overall = np.vstack([overall,len_stats(intergenic_region_length,0,0)])
	overall = np.vstack([overall,len_stats(exon_region_length,num_loci,num_mRNAs)])
	overall = np.vstack([overall,count_stats(exons_per_mRNA)])
	overall = np.vstack([overall,len_stats(intron_region_length,num_loci,num_mRNAs)])
	overall = np.vstack([overall,len_stats(UTRp5_region_length,num_loci,num_mRNAs)])
	overall = np.vstack([overall,count_stats(UTRp5_exons_per_mRNA)])
	overall = np.vstack([overall,len_stats(UTRp3_region_length,num_loci,num_mRNAs)])
	overall = np.vstack([overall,count_stats(UTRp3_exons_per_mRNA)])
	overall = np.vstack([overall,len_stats(CDS_region_length,num_loci,num_mRNAs)])
	overall = np.vstack([overall,count_stats(CDS_exons_per_mRNA)])
	overall = np.vstack([overall,len_stats(mRNA_length,num_loci,num_mRNAs)])
	overall = np.vstack([overall,len_stats(CDS_length,num_loci,num_mRNAs)])
	overall = np.vstack([overall,len_stats(UTRp5_length,num_loci,num_mRNAs)])
	overall = np.vstack([overall,len_stats(UTRp3_length,num_loci,num_mRNAs)])
	overall = np.vstack([overall,len_stats(protein_length,num_loci,num_mRNAs)])

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

	if print_mems :

		print >> filename_stats, "\n\n"
		print >> filename_stats, "#########################################################"
		print >> filename_stats, "#### Mono-exonic ########################################"

		num_loci = len(gene_region_length_singleExon)
		num_mRNAs = len(mRNA_region_length_singleExon)

		mono = np.array(len_stats(gene_region_length_singleExon,num_loci,0))
		mono = np.vstack([mono,len_stats(mRNA_region_length_singleExon,num_loci,num_mRNAs)])
		mono = np.vstack([mono,count_stats(mRNA_per_gene_singleExon)])
		mono = np.vstack([mono,len_stats(exon_region_length_singleExon,num_loci,num_mRNAs)])
		mono = np.vstack([mono,len_stats(UTRp5_region_length_singleExon,num_loci,num_mRNAs)])
		mono = np.vstack([mono,len_stats(UTRp3_region_length_singleExon,num_loci,num_mRNAs)])
		mono = np.vstack([mono,len_stats(CDS_region_length_singleExon,num_loci,num_mRNAs)])
		mono = np.vstack([mono,len_stats(mRNA_length_singleExon,num_loci,num_mRNAs)])
		mono = np.vstack([mono,len_stats(CDS_length_singleExon,num_loci,num_mRNAs)])
		mono = np.vstack([mono,len_stats(UTRp5_length_singleExon,num_loci,num_mRNAs)])
		mono = np.vstack([mono,len_stats(UTRp3_length_singleExon,num_loci,num_mRNAs)])
		mono = np.vstack([mono,len_stats(protein_length_singleExon,num_loci,num_mRNAs)])

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

		multi = np.array(len_stats(gene_region_length_multiExon,num_loci,0))
		multi = np.vstack([multi,len_stats(mRNA_region_length_multiExon,num_loci,num_mRNAs)])
		multi = np.vstack([multi,count_stats(mRNA_per_gene_multiExon)])
		multi = np.vstack([multi,len_stats(exon_region_length_multiExon,num_loci,num_mRNAs)])
		multi = np.vstack([multi,count_stats(exons_per_mRNA_multiExon)])
		multi = np.vstack([multi,len_stats(UTRp5_region_length_multiExon,num_loci,num_mRNAs)])
		multi = np.vstack([multi,count_stats(UTRp5_exons_per_mRNA_multiExon)])
		multi = np.vstack([multi,len_stats(UTRp3_region_length_multiExon,num_loci,num_mRNAs)])
		multi = np.vstack([multi,count_stats(UTRp3_exons_per_mRNA_multiExon)])
		multi = np.vstack([multi,len_stats(CDS_region_length_multiExon,num_loci,num_mRNAs)])
		multi = np.vstack([multi,count_stats(CDS_exons_per_mRNA_multiExon)])
		multi = np.vstack([multi,len_stats(mRNA_length_multiExon,num_loci,num_mRNAs)])
		multi = np.vstack([multi,len_stats(CDS_length_multiExon,num_loci,num_mRNAs)])
		multi = np.vstack([multi,len_stats(UTRp5_length_multiExon,num_loci,num_mRNAs)])
		multi = np.vstack([multi,len_stats(UTRp3_length_multiExon,num_loci,num_mRNAs)])
		multi = np.vstack([multi,len_stats(protein_length_multiExon,num_loci,num_mRNAs)])

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


def plot_distributions( reference , length_lists , subfeat_count_list , prefix , print_mems=True ) :

	prefix = str(prefix)
	gene_region_length , mRNA_region_length , exon_region_length , UTRp5_region_length , UTRp3_region_length , CDS_region_length , intron_region_length , intergenic_region_length , mRNA_length , UTRp5_length , UTRp3_length , CDS_length , protein_length , gene_region_length_singleExon , mRNA_region_length_singleExon , exon_region_length_singleExon , UTRp5_region_length_singleExon , UTRp3_region_length_singleExon , CDS_region_length_singleExon , mRNA_length_singleExon , UTRp5_length_singleExon , UTRp3_length_singleExon , protein_length_singleExon , CDS_length_singleExon , gene_region_length_multiExon , mRNA_region_length_multiExon , exon_region_length_multiExon , UTRp5_region_length_multiExon , UTRp3_region_length_multiExon , CDS_region_length_multiExon , mRNA_length_multiExon , UTRp5_length_multiExon , UTRp3_length_multiExon , protein_length_multiExon , CDS_length_multiExon = length_lists
	mRNA_per_gene , exons_per_mRNA , UTRp5_exons_per_mRNA , UTRp3_exons_per_mRNA , CDS_exons_per_mRNA , mRNA_per_gene_singleExon , mRNA_per_gene_multiExon , exons_per_mRNA_multiExon , CDS_exons_per_mRNA_multiExon , UTRp5_exons_per_mRNA_multiExon , UTRp3_exons_per_mRNA_multiExon , UTRp3_exons_per_mRNA_multiExon = subfeat_count_list

	genome_seq_lenths = sorted([ len(reference[x].seq) for x in reference ] ,reverse=True)

	with PdfPages( prefix + ".distribution.pdf") as pdf:
		plt.style.use("ggplot")

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist(genome_seq_lenths , 200, cumulative=True, normed=True, edgecolor='#444444', color='#444444')
		plt.title("Length of genomic sequences", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist(gene_region_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("Loci length", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( mRNA_region_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("mRNA lengths", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( mRNA_per_gene , max( [ 1 , max(mRNA_per_gene) - min(mRNA_per_gene) ] )  , edgecolor='#E6E6E6', color='#444444')
		plt.title("mRNAs per locus", fontsize=8)
		plt.xlabel("Count" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( intergenic_region_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("Intergenic regions length", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( exon_region_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("Exon length", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( exons_per_mRNA , max( [ 1 , max(exons_per_mRNA)-min(exons_per_mRNA)] ) , edgecolor='#E6E6E6', color='#444444')
		plt.title("Exon per mRNA", fontsize=8)
		plt.xlabel("Count" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( intron_region_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("Intron length", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( UTRp5_region_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("5' UTRs length", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( UTRp5_exons_per_mRNA , max( [ 1 , max(UTRp5_exons_per_mRNA)-min(UTRp5_exons_per_mRNA) ] )  , edgecolor='#E6E6E6', color='#444444')
		plt.title("5' UTRs per mRNA", fontsize=8)
		plt.xlabel("Count" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( UTRp3_region_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("3' UTRs length", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( UTRp3_exons_per_mRNA , max( [ 1 , max(UTRp3_exons_per_mRNA)-min(UTRp3_exons_per_mRNA) ] ) , edgecolor='#E6E6E6', color='#444444')
		plt.title("3' UTRs per mRNA", fontsize=8)
		plt.xlabel("Count" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( CDS_region_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("CDS length", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( CDS_exons_per_mRNA , max( [ 1 , max(CDS_exons_per_mRNA)-min(CDS_exons_per_mRNA) ] ) , edgecolor='#E6E6E6', color='#444444')
		plt.title("CDS per mRNA", fontsize=8)
		plt.xlabel("Count" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( mRNA_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("mRNA sequence length", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( CDS_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("Coding sequence length", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( UTRp5_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("5' UTR sequence length", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( UTRp3_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("3' UTR sequence length", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()

		plt.figure()
		plt.xticks(fontsize=6,rotation=45)
		plt.yticks(fontsize=6)
		n, bins, patches = plt.hist( protein_length , 50, edgecolor='#E6E6E6', color='#444444')
		plt.title("Protein length", fontsize=8)
		plt.xlabel("Length (bp)" , fontsize=7 )
		plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
		plt.close()


	if print_mems :
		with PdfPages(prefix+".distribution.monoexonic.pdf") as pdf:
			plt.style.use("ggplot")

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( gene_region_length_singleExon , 50, edgecolor='#444444', color='#444444')
			plt.title("Loci length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( mRNA_region_length_singleExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("mRNA lengths", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( mRNA_per_gene_singleExon , max( [ 1 , max(mRNA_per_gene_singleExon)-min(mRNA_per_gene_singleExon) ] ) , edgecolor='#E6E6E6', color='#444444')
			plt.title("mRNAs per locus", fontsize=8)
			plt.xlabel("Count" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( exon_region_length_singleExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("Exon length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( UTRp5_region_length_singleExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("5' UTRs length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( UTRp3_region_length_singleExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("3' UTRs length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( CDS_region_length_singleExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("CDS length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( mRNA_length_singleExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("mRNA sequence length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( CDS_length_singleExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("Coding sequence length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( UTRp5_length_singleExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("5' UTR sequence length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( UTRp3_length_singleExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("3' UTR sequence length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( protein_length_singleExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("Protein length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

		with PdfPages(prefix+".distribution.multiexonic.pdf") as pdf:
			plt.style.use("ggplot")

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist(gene_region_length_multiExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("Loci length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( mRNA_region_length_multiExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("mRNA lengths", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( mRNA_per_gene_multiExon , max( [ 1 , max(mRNA_per_gene_multiExon)-min(mRNA_per_gene_multiExon) ] ) , edgecolor='#E6E6E6', color='#444444')
			plt.title("mRNAs per locus", fontsize=8)
			plt.xlabel("Count" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( exon_region_length_multiExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("Exon length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( exons_per_mRNA_multiExon , max( [ 1 , max(exons_per_mRNA_multiExon)-min(exons_per_mRNA_multiExon) ] ) , edgecolor='#E6E6E6', color='#444444')
			plt.title("Exon per mRNA", fontsize=8)
			plt.xlabel("Count" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( UTRp5_region_length_multiExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("5' UTRs length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( UTRp5_exons_per_mRNA_multiExon , max( [ 1 , max(UTRp5_exons_per_mRNA_multiExon)-min(UTRp5_exons_per_mRNA_multiExon) ] ) , edgecolor='#E6E6E6', color='#444444')
			plt.title("5' UTRs per mRNA", fontsize=8)
			plt.xlabel("Count" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( UTRp3_region_length_multiExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("3' UTRs length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( UTRp3_exons_per_mRNA_multiExon , max( [ 1 , max(UTRp3_exons_per_mRNA_multiExon)-min(UTRp3_exons_per_mRNA_multiExon) ] ) , edgecolor='#E6E6E6', color='#444444')
			plt.title("3' UTRs per mRNA", fontsize=8)
			plt.xlabel("Count" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( CDS_region_length_multiExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("CDS length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( CDS_exons_per_mRNA_multiExon , max( [ 1 , max(CDS_exons_per_mRNA_multiExon)-min(CDS_exons_per_mRNA_multiExon) ] )  , edgecolor='#E6E6E6', color='#444444')
			plt.title("CDS per mRNA", fontsize=8)
			plt.xlabel("Count" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( mRNA_length_multiExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("mRNA sequence length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( CDS_length_multiExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("Coding sequence length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( UTRp5_length_multiExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("5' UTR sequence length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( UTRp3_length_multiExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("3' UTR sequence length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()

			plt.figure()
			plt.xticks(fontsize=6,rotation=45)
			plt.yticks(fontsize=6)
			n, bins, patches = plt.hist( protein_length_multiExon , 50, edgecolor='#E6E6E6', color='#444444')
			plt.title("Protein length", fontsize=8)
			plt.xlabel("Length (bp)" , fontsize=7 )
			plt.tight_layout() ; pdf.savefig()  # saves the current figure into a pdf page
			plt.close()


def new_CDS( mRNA_id , count , old_CDS_entry ) :

	seqname, source, feature, start, end, score, strand, frame, attribute = old_CDS_entry[0].split("\t")
	attribute = "ID=" + mRNA_id + ".cds" + str(count) + ";Parent=" + mRNA_id
	new_line = "\t".join([seqname, source, feature, str(start), str(end), score, strand, frame, attribute])

	return [ new_line , int(start), int(end) ]


def new_exon( coords , strand , count , mRNA_id , seqname, source ) :
	start , end = coords
	attribute = "ID=" + mRNA_id + ".exon" + str(count) + ";Parent=" + mRNA_id
	new_line = "\t".join([seqname, source, "exon", str(start), str(end), ".", strand, ".", attribute])

	return [ new_line , int(start), int(end) ]


def new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ) :
	attribute = "ID=" + mRNA_id + "." + feature + str(count) + ";Parent=" + mRNA_id
	new_line = "\t".join([seqname, source, feature, str(start), str(end), ".", strand, ".", attribute])

	return [ new_line , int(start) , int(end) ]


def make_UTRs( seqname, source , strand , exon_coord_list , CDS_coord_list , mRNA_id ) :
	five_prime_UTR = []
	three_prime_UTR = []
	CDS_coords = sorted(CDS_coord_list[:])
	CDS_range = [ CDS_coords[0][0] , CDS_coords[-1][1] ]
	count = 0
	if strand == "+":
		# 5' -> 3'
		for exon in exon_coord_list:
			found_cds = False
			count += 1

			if exon[1] < CDS_range[0] :
				# Exon starts and ends before the CDS region
				feature = "five_prime_UTR"
				start , end = exon
				five_prime_UTR.append(new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ))

			elif exon[0] < CDS_range[0] and exon[1] < CDS_range[1] :
				# Exon starts before CDS and ends within.
				# 5' utr extends up to CDS begin
				feature = "five_prime_UTR"
				start = exon[0]
				end = CDS_range[0] - 1
				five_prime_UTR.append(new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ))
				found_cds = True

			elif exon[0] < CDS_range[0] and exon[1] > CDS_range[1] :
				# Exon starts before and ends after
				# Add UTRs on both sides
				# 5'
				feature = "five_prime_UTR"
				start = exon[0]
				end = CDS_range[0] - 1
				five_prime_UTR.append(new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ))
				# 3'
				feature = "three_prime_UTR"
				count = 1
				end = exon[1]
				start = CDS_range[1] + 1
				three_prime_UTR.append(new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ))
				found_cds = True

			elif exon[0] > CDS_range[0] and exon[1] < CDS_range[1] :
				# Exon is within CDS
				# Do nothing
				feature = "CDS"
				found_cds = True
				count = -1

			elif exon[0] < CDS_range[1] and exon[1] > CDS_range[1] :
				# Exon start inside CDS and ends after
				# 3' utr starts where CDS ends
				feature = "three_prime_UTR"
				count = 1
				end = exon[1]
				start = CDS_range[1] + 1
				three_prime_UTR.append(new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ))
				found_cds = True

			elif exon[0] > CDS_range[1]  :
				# Exon is downstream CDS region
				if found_cds and count == 0 :
					count=1
				feature = "three_prime_UTR"
				start , end = exon
				three_prime_UTR.append(new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ))

	else:
		# 3' -> 5'
		for exon in exon_coord_list:
			found_cds = False
			count += 1

			if exon[1] < CDS_range[0] :
				# Exon starts and ends before the CDS region
				feature = "three_prime_UTR"
				start , end = exon
				three_prime_UTR.append(new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ))

			elif exon[0] < CDS_range[0] and exon[1] < CDS_range[1] :
				# Exon starts before CDS and ends within.
				# 5' utr extends up to CDS begin
				feature = "three_prime_UTR"
				start = exon[0]
				end = CDS_range[0] - 1
				three_prime_UTR.append(new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ))
				found_cds = True

			elif exon[0] < CDS_range[0] and exon[1] > CDS_range[1] :
				# Exon starts before and ends after
				# Add UTRs on both sides
				# 3'
				feature = "three_prime_UTR"
				start = exon[0]
				end = CDS_range[0] - 1
				three_prime_UTR.append(new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ))
				# 5'
				feature = "five_prime_UTR"
				count = 1
				end = exon[1]
				start = CDS_range[1] + 1
				five_prime_UTR.append(new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ))
				found_cds = True

			elif exon[0] > CDS_range[0] and exon[1] < CDS_range[1] :
				# Exon is within CDS
				# Do nothing
				feature = "CDS"
				found_cds = True
				count = -1

			elif exon[0] < CDS_range[1] and exon[1] > CDS_range[1] :
				# Exon start inside CDS and ends after
				# 3' utr starts where CDS ends
				feature = "five_prime_UTR"
				count = 1
				end = exon[1]
				start = CDS_range[1] + 1
				five_prime_UTR.append(new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ))
				found_cds = True

			elif exon[0] > CDS_range[1]  :
				# Exon is downstream CDS region
				if found_cds and count == 0 :
					count=1
				feature = "five_prime_UTR"
				start , end = exon
				five_prime_UTR.append(new_utr( start , end , feature , mRNA_id , seqname, source , strand , count ))



	return five_prime_UTR , three_prime_UTR


def mRNA_consensus( new_name , matching_ids , matching_structures , CDS , gene_id , seqname, source, strand ):

	mRNA_feat_order = {"exon": 1, "intron":1, "five_prime_UTR":2 , "CDS":3, "three_prime_UTR":4}
	new_mRNA_id = new_name
	consensus_model = [ "" , "" , {} ]
	consensus_splice_sites = []
	mRNA_seqname = seqname
	mRNA_source = source


	# Make exons consensus
	consensus_model[2][1] = []
	all_exons = []
	#print >> sys.stderr , "matching_structures"
	#print >> sys.stderr , matching_structures
	for coords in matching_structures :
		coords_list = [ ( coords[i*2] , coords[i*2 +1] ) for i in range((len(coords) + 1)/2)]
		#print >> sys.stderr , coords_list
		all_exons += coords_list

	all_exons = sorted(set(all_exons))
	#print >> sys.stderr , "all_exons"
	#print >> sys.stderr , all_exons
	actual_exon = ""
	for interval in all_exons :
		#print >> sys.stderr , interval
		if actual_exon == "" :
			actual_exon = list(interval)
		else :
			if interval[0] < actual_exon[1] :
				# new exon overlaps the previous exon
				# Update exon range
				actual_exon = [ min(interval[0] , actual_exon[0]) , max(interval[1], actual_exon[1]) ]
			else :
				consensus_splice_sites.append(actual_exon)
				actual_exon = list(interval)
	consensus_splice_sites.append(actual_exon)
	consensus_splice_sites = sorted(consensus_splice_sites)
	#print >> sys.stderr , "consensus_splice_sites"
	#print >> sys.stderr , consensus_splice_sites

	counter = 0
	if strand == "+":
		for interval in consensus_splice_sites:
			counter+=1
			consensus_model[2][1].append(new_exon( interval , strand , counter , new_mRNA_id , mRNA_seqname, mRNA_source))
	else :
		for interval in consensus_splice_sites[::-1]:
			counter+=1
			consensus_model[2][1].append(new_exon( interval , strand , counter , new_mRNA_id , mRNA_seqname, mRNA_source ))

	mRNA_start = int(consensus_splice_sites[0][0])
	mRNA_stop = int(consensus_splice_sites[-1][-1])
	#print >> sys.stderr, [mRNA_start , mRNA_stop]
	# Add CDSs from mRNA_1. CDS are identical in the two models
	consensus_model[2][3] = []
	counter=0
	cds_ranges = []
	if strand == "+":
		for entry in sorted(CDS , key=lambda feat : feat[1]) :
			counter+=1
			new_cds_entry = new_CDS(new_mRNA_id , counter , entry)
			cds_ranges.append([new_cds_entry[1],new_cds_entry[2]])
			consensus_model[2][3].append(new_cds_entry)
	else :
		for entry in sorted(CDS , key=lambda feat : feat[1], reverse=True ) :
			counter+=1
			new_cds_entry = new_CDS(new_mRNA_id , counter , entry)
			cds_ranges.append([new_cds_entry[1],new_cds_entry[2]])
			consensus_model[2][3].append(new_cds_entry)

	# Make mRNA consensus
	mRNA_attribute = "ID=" + new_mRNA_id + ";Parent=" + gene_id
	consensus_model[1] = int(mRNA_start)
	consensus_model[0] = "\t".join([mRNA_seqname, mRNA_source, "mRNA", str(mRNA_start), str(mRNA_stop), ".", strand, ".", mRNA_attribute])

	# Extract UTR information
	consensus_model[2][2] , consensus_model[2][4]  = make_UTRs( mRNA_seqname, mRNA_source , strand , consensus_splice_sites , cds_ranges , new_mRNA_id)  # five_prime_UTR , three_prime_UTR

	return consensus_model


def update_coords_gff_search( line , offset_dict ) :
	seqname, source, feature, start, end, score, strand, frame, attribute = line.split("\t")
	gene_offset = []
	# Search for right set of new_chr, offset and direction from the ranges available in offset_list
	for element in sorted(offset_dict.keys()) :
		if int(start) >= int(element[0]) and int(end) <= int(element[1]) :
			new_chr , offset , direction = offset_dict[element]

			if direction == "+" :
				new_start = min( ( int(start) + offset ) , ( int(end) + offset ) )
				new_end = max( ( int(start) + offset ) , ( int(end) + offset ) )
				new_strand = strand
				new_line = '\t'.join([str(x) for x in [new_chr, source, feature, new_start, new_end, score, new_strand, frame, attribute]])
				gene_offset.append( [ new_line , new_chr , offset , direction , new_start , new_end ] )
			elif direction == "-" :
				new_start = min( (-int(start) + offset) , (-int(end) + offset ) )
				new_end = max( (-int(start) + offset) , (-int(end) + offset ) )
				if strand == "-" :
					new_strand = "+"
				else :
					new_strand = "-"
				new_line = '\t'.join([str(x) for x in [new_chr, source, feature, new_start, new_end, score, new_strand, frame, attribute]])
				gene_offset.append( [ new_line , new_chr , offset , direction , new_start , new_end ] )

	if gene_offset == [] :
		return ""
	else :
		return gene_offset


def update_coords_gff_from_offset( line , new_chr , offset , direction ) :
	seqname, source, feature, start, end, score, strand, frame, attribute = line.split("\t")

	if direction == "+" :
		new_start = min( ( int(start) + offset ) , ( int(end) + offset ) )
		new_end = max( ( int(start) + offset ) , ( int(end) + offset ) )
		new_strand = strand
		new_line = '\t'.join([str(x) for x in [new_chr, source, feature, new_start, new_end, score, new_strand, frame, attribute]])
		return new_line , new_chr , new_start , new_end
	elif direction == "-" :
		new_start = min( (-int(start) + offset) , (-int(end) + offset ) )
		new_end = max( (-int(start) + offset) , (-int(end) + offset ) )
		if strand == "-" :
			new_strand = "+"
		else :
			new_strand = "-"
		new_line = '\t'.join([str(x) for x in [new_chr, source, feature, new_start, new_end, score, new_strand, frame, attribute]])
		return new_line , new_chr , new_start , new_end

	else :  # direction = "" , no compliant range found -> untranslatable -> broken
		return ""


def translate_gff3( gff3_db , coordinate_offset_db , broken_file , multiple_copies_file="gene_ids.multiple_copy.txt" ) :
	# Coordinate offset dict direct format:
	# coordinate_offset_db[old_seq_id][(old_seq_start , old_seq_stop)] = [ new_seq_id , offset , strand ]
	new_gff3 = {}
	multiplied_list = {}
	multiple_copies = open( multiple_copies_file , "w" )
	# Update each element in gff3_db
	### Fields to update :
	#### GFF line - chr = el[0] , start = el[3] , stop = el[4])
	#### db attributes: [chr] , int(start)
	### genes db structure:
	#### genes[gene_id] = [ gff_line , chr , int(start) , mRNA_dict{} ]
	####		mRNA_dict[mRNA_id] = [ gff_line , int(start), feat_dict{} ]
	####			feat_dict[mRNA_order_num] = [ [ line , int(start) , int(end) ]  ]
	####		 		with mRNA_order_num meaning:	1 = "exon"
	####												2 = "five_prime_UTR"
	####												3 = "CDS"
	####												4 = "three_prime_UTR"
	#### line (unsplitted on tab) = seqname, source, feature, start, end, score, strand, frame, attribute

	drop_file = open(broken_file, "w")

	for gene_id in gff3_db :
		suffix_id = -1
		# update gene feat
		gff_line , chr , start , mRNA_dict = gff3_db[gene_id]
		if chr in coordinate_offset_db :
			updated = update_coords_gff_search( gff_line , coordinate_offset_db[chr] )
			# Returns a list of gene loci on the destination.
		else :
			updated = ""

		if updated == "" :
			# update_coords_gff returned an empty value because the gene was not entirely falling in any good AGP range
			# The locus may have been broken between 2 regions or falling in an unplaced region
			# Locus must be drop (along with the rest of the annotation)
			print >> sys.stderr , "[Warning] Gene locus dropped: " + gene_id
			print >> drop_file , gene_id
			continue
		else :
			# New region found for the locus
			# Save the offset info and pass them to the subfeatures
			# NEW: handles duplicated ids -> uniquify names by adding a literal suffix
			# 	updated may be a list of list of multiple destination loci
			# 	handle them separately and update names
			if len(updated) == 1 :
				# Only one destination
				new_gff_line , new_chr , offset , direction , new_start , new_stop = updated[0]
				new_gff3[gene_id] = [ new_gff_line , new_chr , int(new_start) , {} ]
				# update gene subfeatures
				for mRNA_id in mRNA_dict.keys():
					# update mrna feat
					mRNA_gff_line , mRNA_start, feat_dict = mRNA_dict[mRNA_id]
					new_mRNA_gff_line , new_chr , new_mRNA_start , new_mRNA_end = update_coords_gff_from_offset( mRNA_gff_line , new_chr , offset , direction )
					new_gff3[gene_id][3][mRNA_id] = [ new_mRNA_gff_line , int(new_mRNA_start) , {} ]
					# update mRNA subfeat
					for feat_id in feat_dict.keys():
						# for each mRNA subfeat type:
						new_gff3[gene_id][3][mRNA_id][2][feat_id] = []
						feat_list = feat_dict[feat_id]
						for element in feat_list :
							# for each subfeat of mRNA subfeat category
							feat_line , feat_start , feat_end = element
							new_feat_line , new_chr , new_feat_start , new_feat_end = update_coords_gff_from_offset( feat_line , new_chr , offset , direction )
							new_gff3[gene_id][3][mRNA_id][2][feat_id].append( [ new_feat_line , int(new_feat_start) , int(new_feat_end) ] )
			else :
				# Multiple destinations, make multiple copies of the gene with different ids
				print >> sys.stderr , "[DEBUG] Gene: " + gene_id + " has " + str(len(updated)) + " destinations "
				alphabet = "abcdefghijklmnopqrstuvwxyz"
				multiplied_list[gene_id] = []
				for gene_destination in updated :
					suffix_id += 1
					try :
						suffix = alphabet[suffix_id]
					except :
						print >> sys.stderr , "[ERROR] Suffix index out of range of alphabet for a gene:"
						print >> sys.stderr , "[ERROR] Suffix: " + str(suffix_id)
						print >> sys.stderr , "[ERROR] Gene: " + gene_id
						sys.exit(1)
					renamed_gene_id = gene_id + "." + suffix
					multiplied_list[gene_id].append(renamed_gene_id)
					new_gff_line , new_chr , offset , direction , new_start , new_stop = gene_destination
					renamed_gff_line = new_gff_line.replace( gene_id , renamed_gene_id )
					new_gff3[renamed_gene_id] = [ renamed_gff_line , new_chr , int(new_start) , {} ]
					for mRNA_id in mRNA_dict.keys():
						# update mrna feat
						mRNA_gff_line , mRNA_start, feat_dict = mRNA_dict[mRNA_id]
						new_mRNA_gff_line , new_chr , new_mRNA_start , new_mRNA_end = update_coords_gff_from_offset( mRNA_gff_line , new_chr , offset , direction )
						old_mRNA_id = gff_feature_id(new_mRNA_gff_line)
						renamed_mRNA_gff_line = new_mRNA_gff_line.replace( gene_id , renamed_gene_id )
						renamed_mRNA_id = gff_feature_id(renamed_mRNA_gff_line)
						renamed_mRNA = True
						if renamed_mRNA_id == old_mRNA_id :
							renamed_mRNA = False
							renamed_mRNA_id = old_mRNA_id + "." + suffix
							actual_gff_line = renamed_mRNA_gff_line
							renamed_mRNA_gff_line = actual_gff_line.replace( old_mRNA_id , renamed_mRNA_id )
						new_gff3[renamed_gene_id][3][renamed_mRNA_id] = [ renamed_mRNA_gff_line , int(new_mRNA_start) , {} ]
						# update mRNA subfeat

						for feat_id in feat_dict.keys():
							# for each mRNA subfeat type:
							new_gff3[renamed_gene_id][3][renamed_mRNA_id][2][feat_id] = []
							feat_list = feat_dict[feat_id]
							for element in feat_list :
								# for each subfeat of mRNA subfeat category
								feat_line , feat_start , feat_end = element
								new_feat_line , new_chr , new_feat_start , new_feat_end = update_coords_gff_from_offset( feat_line , new_chr , offset , direction )
								renamed_feat_line = new_feat_line.replace( gene_id , renamed_gene_id )
								if not renamed_mRNA :
									actual_gff_line = renamed_feat_line
									renamed_feat_line = actual_gff_line.replace( old_mRNA_id , renamed_mRNA_id )
								# Check if the subfeature has an id that needs to be uniquified
								old_feat_id = gff_feature_id(new_feat_line)
								renamed_feat_id = gff_feature_id(renamed_feat_line)
								if (old_feat_id == renamed_feat_id) and not ( old_feat_id == "" ) :
									renamed_feat_id = old_feat_id + "." + suffix
									actual_gff_line = renamed_feat_line
									renamed_feat_line = actual_gff_line.replace( old_feat_id , renamed_feat_id )

								new_gff3[renamed_gene_id][3][renamed_mRNA_id][2][feat_id].append( [ renamed_feat_line , int(new_feat_start) , int(new_feat_end) ] )

	for id in multiplied_list :
		print >> multiple_copies, id + "\t" + ",".join( [ str(x) for x in multiplied_list[id] ] )
	drop_file.close()
	multiple_copies.close()
	return new_gff3


def gff_feature_id( gff_line ) :
	attributes_dict = {}
	seqname, source, feature, start, end, score, strand, frame, attribute = gff_line.rstrip().split("\t")
	for chunk in attribute.rstrip(";").split(";"):
		att = chunk.split("=")
		try :
			attributes_dict[str(att[0]).upper()] = att[1]
		except :
			print >> sys.stderr , "[ERROR] GFF attributes structure are not compliant at line: " + gff_line
			exit(1)
	if "ID" in attributes_dict :
		return attributes_dict["ID"]
	else :
		return ""


def get_sequence( gff3_db , fasta_db , filename_prefix , feat = "CDS" ) :
	filename = filename_prefix + "." + feat + ".fasta"
	sequence_fasta_file = open(filename , "w")
	sequence_fasta_db = {}

	if feat == "CDS" or feat == "protein"  :
		feat_type = 3
	elif feat == "five_prime_UTR" :
		feat_type = 2
	elif feat == "three_prime_UTR" :
		feat_type = 4
	else :
		feat_type = 1

	for gene_id in gff3_db.keys() :
		gff_line , chr , start , mRNA_dict = gff3_db[gene_id]
		for mRNA_id in mRNA_dict.keys():
			mRNA_seq = ""
			mRNA_gff_line , mRNA_start, feat_dict = mRNA_dict[mRNA_id]
			# select feature kind

			if feat_type not in feat_dict :
				print >> sys.stderr , "[WARNING] The selected feature is missing for the transcript " + mRNA_id + ". It will be skipped. (" +  mRNA_gff_line + ")"
				continue

			for feat in sorted(feat_dict[feat_type] , key=lambda x: x[1]) :
				feat_line , feat_start , feat_end = feat
				#print >> sys.stderr, feat_line
				seqname, source, feature, start, end, score, strand, old_phase, attribute = feat_line.split("\t")
				if not strand=="-":
					mRNA_seq += fasta_db[seqname][int(start)-1:int(end)]
				else :
					mRNA_seq = str( Seq(fasta_db[seqname][int(start)-1:int(end)]).reverse_complement() ) + mRNA_seq

			if feat == "protein" :
				sequence_fasta_db[mRNA_id] = str(Seq(str(mRNA_seq.seq)).translate())
			else :
				sequence_fasta_db[mRNA_id] = mRNA_seq

	for mRNA_id in sorted(sequence_fasta_db.keys()) :
		print >> sequence_fasta_file, ">" + mRNA_id
		print >> sequence_fasta_file, str(sequence_fasta_db[mRNA_id]).rstrip()

	sequence_fasta_file.close()

	return filename


def get_gene2mRNA( gff3 ) :
	mRNA_db = {}

	for line in open(gff3, 'r') :
		if line.rstrip() == "" or line[0] == "#" :
			continue
		else :
			seqname, source, feature, start, end, score, strand, frame, attribute = line.rstrip().split("\t")
		if not feature=="mRNA" :
			# Skip all lines not of features different from mRNA
			continue
		else :
			attributes_dict = {}
			for chunk in attribute.rstrip(";").split(";"):
					att = chunk.split("=")
					attributes_dict[att[0]] = att[1]

			mRNA_db[attributes_dict["ID"]] = attributes_dict["Parent"]

	return mRNA_db


def get_gene2mRNA_from_db( gff_db ) :
	mRNA_db = {}
	for gene_id in gff_db.keys() :
		gff_line , chr , start , mRNA_dict = gff_db[gene_id]
		for mRNA_id in mRNA_dict.keys():
			mRNA_db[mRNA_id] = gene_id
	return mRNA_db


def gff3_filter2table_Hap(gff3_db, selected_feature, chr_name_pattern) :

	filtered_table = []

	#### GFF line - chr = el[0] , start = el[3] , stop = el[4])
	#### db attributes: [chr] , int(start)
	### genes db structure:
	#### genes[gene_id] = [ gff_line , chr , int(start) , mRNA_dict{} ]
	####		mRNA_dict[mRNA_id] = [ gff_line , int(start), feat_dict{} ]
	####			feat_dict[mRNA_order_num] = [ [ line , int(start) , int(end) ]  ]
	####		 		with mRNA_order_num meaning:	1 = "exon"
	####												2 = "five_prime_UTR"
	####												3 = "CDS"
	####												4 = "three_prime_UTR"
	#### line (unsplitted on tab) = seqname, source, feature, start, end, score, strand, frame, attribute

	for gene_id in gff3_db.keys() :
		gene_line , chr , start , mRNA_dict = gff3_db[gene_id]
		if not chr_name_pattern in chr :
			continue
		else :
			if selected_feature == "gene" :
				seqname, source, feature, start, end, score, strand, frame, attribute = gene_line.rstrip().split("\t")
				filtered_table.append( [chr, int(start) , int(end) , gene_id] )
				continue

			else :
				for mRNA_id in mRNA_dict.keys() :
					mRNA_line , start, feat_dict = mRNA_dict[mRNA_id]

					if selected_feature == "mRNA" :
						seqname, source, feature, start, end, score, strand, frame, attribute = gene_line.rstrip().split("\t")
						filtered_table.append( [chr, int(start) , int(end) , gene_id] )
						continue

					else :
						if selected_feature == "CDS" : subfeat=3
						elif selected_feature == "five_prime_UTR" : subfeat=2
						elif selected_feature == "three_prime_UTR" : subfeat=4
						else : subfeat=1 #exon
						for element in feat_dict[subfeat] :
							subfeat_line , start , end  = element
							seqname, source, feature, start, end, score, strand, frame, attribute = subfeat_line.rstrip().split("\t")
							attributes_dict = {}
							for chunk in attribute.rstrip(";").split(";"):
								att = chunk.split("=")
								attributes_dict[att[0]] = att[1]
							filtered_table.append( [chr, int(start) , int(end) , attributes_dict["ID"]] )
	return filtered_table


def gff3_filter2table(gff3_db, selected_feature, column , pairs) :
	filtered_table = []

	#### GFF line - chr = el[0] , start = el[3] , stop = el[4])
	#### db attributes: [chr] , int(start)
	### genes db structure:
	#### genes[gene_id] = [ gff_line , chr , int(start) , mRNA_dict{} ]
	####		mRNA_dict[mRNA_id] = [ gff_line , int(start), feat_dict{} ]
	####			feat_dict[mRNA_order_num] = [ [ line , int(start) , int(end) ]  ]
	####		 		with mRNA_order_num meaning:	1 = "exon"
	####												2 = "five_prime_UTR"
	####												3 = "CDS"
	####												4 = "three_prime_UTR"
	#### line (unsplitted on tab) = seqname, source, feature, start, end, score, strand, frame, attribute

	selected_sequences = []
	for pair in pairs:
		try :
			selected_id = pair[ int(column) - 1 ]
		except :
			print >> sys.stderr , "[ERROR] Unexpected selection of haplotype: " + column
			sys.exit(1)
		else :
			selected_sequences.append(selected_id)

	for gene_id in gff3_db.keys() :
		gene_line , chr , start , mRNA_dict = gff3_db[gene_id]
		if chr not in selected_sequences :
			continue
		else :
			if selected_feature == "gene" :
				seqname, source, feature, start, end, score, strand, frame, attribute = gene_line.rstrip().split("\t")
				filtered_table.append( [chr, int(start) , int(end) , gene_id] )
				continue

			else :
				for mRNA_id in mRNA_dict.keys() :
					mRNA_line , start, feat_dict = mRNA_dict[mRNA_id]

					if selected_feature == "mRNA" :
						seqname, source, feature, start, end, score, strand, frame, attribute = gene_line.rstrip().split("\t")
						filtered_table.append( [chr, int(start) , int(end) , gene_id] )
						continue

					else :
						if selected_feature == "CDS" : subfeat=3
						elif selected_feature == "five_prime_UTR" : subfeat=2
						elif selected_feature == "three_prime_UTR" : subfeat=4
						else : subfeat=1 #exon
						for element in feat_dict[subfeat] :
							subfeat_line , start , end  = element
							seqname, source, feature, start, end, score, strand, frame, attribute = subfeat_line.rstrip().split("\t")
							attributes_dict = {}
							for chunk in attribute.rstrip(";").split(";"):
								att = chunk.split("=")
								attributes_dict[att[0]] = att[1]
							filtered_table.append( [chr, int(start) , int(end) , attributes_dict["ID"]] )
	return filtered_table


def feature_ranges( gff3_db , selected_feature ) :

	filtered_table = {}

	#### GFF line - chr = el[0] , start = el[3] , stop = el[4])
	#### db attributes: [chr] , int(start)
	### genes db structure:
	#### genes[gene_id] = [ gff_line , chr , int(start) , mRNA_dict{} ]
	####		mRNA_dict[mRNA_id] = [ gff_line , int(start), feat_dict{} ]
	####			feat_dict[mRNA_order_num] = [ [ line , int(start) , int(end) ]  ]
	####		 		with mRNA_order_num meaning:	1 = "exon"
	####												2 = "five_prime_UTR"
	####												3 = "CDS"
	####												4 = "three_prime_UTR"
	#### line (unsplitted on tab) = seqname, source, feature, start, end, score, strand, frame, attribute

	for gene_id in gff3_db.keys() :
		gene_line , chr , start , mRNA_dict = gff3_db[gene_id]
		if chr not in filtered_table: 
			filtered_table[chr] = []
			
		if selected_feature == "gene" :
			seqname, source, feature, start, end, score, strand, frame, attribute = gene_line.rstrip().split("\t")
			filtered_table[chr].append( [chr, int(start) , int(end) , gene_id] )
			continue

		else :
			for mRNA_id in mRNA_dict.keys() :
				mRNA_line , start, feat_dict = mRNA_dict[mRNA_id]

				if selected_feature == "mRNA" :
					seqname, source, feature, start, end, score, strand, frame, attribute = gene_line.rstrip().split("\t")
					filtered_table[chr].append( [chr, int(start) , int(end) , gene_id] )
					continue

				else :
					if selected_feature == "CDS" : subfeat = 3
					elif selected_feature == "five_prime_UTR" : subfeat = 2
					elif selected_feature == "three_prime_UTR" : subfeat = 4
					else : subfeat = 1  # exon
					for element in feat_dict[subfeat] :
						subfeat_line , start , end = element
						seqname, source, feature, start, end, score, strand, frame, attribute = subfeat_line.rstrip().split("\t")
						attributes_dict = {}
						for chunk in attribute.rstrip(";").split(";"):
							att = chunk.split("=")
							attributes_dict[att[0]] = att[1]
						filtered_table[chr].append( [chr, int(start) , int(end) , attributes_dict["ID"]] )

	return filtered_table


