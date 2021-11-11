import os
import sys
import datetime
from collections import defaultdict
import json
import subprocess
from FASTA_lib import *



def read_paf( paf_file ) :
	original_hits = {}
	id=0
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
		original_hits[(Tid,Qid)].append([ id , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ])
	return original_hits


def map_minimap_unplaced_on_target(region_db, cores , chunk_db, workdir):
	unplaced_on_target_paf = workdir + "/" + region_db["id"] + ".map.unplaced.paf"
	target_fasta = region_db["sequence_file"]
	original_dir = os.path.dirname(chunk_db["sequences"][ chunk_db["sequences"].keys()[0] ]["folder"])
	unplaced_for_file = original_dir + "/tmp.unplaced.for"
	unplaced_rev_file = original_dir + "/tmp.unplaced.rev"
	unplaced_len_file = original_dir + "/tmp.unplaced.len.json"
	if not ( os.path.exists(unplaced_for_file) and os.path.exists(unplaced_rev_file) ) :
		unplaced_len = {}
		print >> sys.stderr, "### original_dir: " + original_dir
		print >> sys.stderr, "#### Temporary unplaced sequences FASTA file: " + unplaced_for_file
		print >> sys.stderr, "#### Temporary unplaced sequences reversed FASTA file: " + unplaced_rev_file
		print >> sys.stderr, "#### Temporary unplaced sequences length database: " + unplaced_len_file
		#make fasta files
		unplaced_for = open(unplaced_for_file, 'w')
		unplaced_rev = open(unplaced_rev_file, 'w')
		unplaced_fasta = read_fasta( chunk_db["inputs"]["U"] )
		for seq in sorted( unplaced_fasta.keys() ) :
			print >> unplaced_for , ">" + seq + "|+"
			print >> unplaced_for , str(Seq(unplaced_fasta[seq])).upper()
			print >> unplaced_rev , ">" + seq + "|-"
			print >> unplaced_rev , str(Seq(unplaced_fasta[seq]).reverse_complement()).upper()
			unplaced_len[ seq+ "|-" ] = len(unplaced_fasta[seq])
			unplaced_len[ seq+ "|+" ] = len(unplaced_fasta[seq])

		unplaced_for.close()
		unplaced_rev.close()

		json.dump(unplaced_len , open(unplaced_len_file, 'w') , indent=4 )
	else :
		unplaced_len = json.load(open(unplaced_len_file))

	query_fasta = [ unplaced_for_file , unplaced_rev_file ]

	return map_minimap( target_fasta , query_fasta , int(cores) , " -x asm20 --for-only " , unplaced_on_target_paf , chunk_db["tool_paths"]["minimap2"] ) , unplaced_len


def map_minimap( ref_file , query_file , cores , parameters , out_file_name , path , filter="" ):
	if path == "" :
		minimap2_search=subprocess.Popen( "which minimap2" , shell=True, stdout=subprocess.PIPE )
		command_line , error = minimap2_search.communicate()
		command_line = command_line.rstrip()
	else :
		command_line = path + "/minimap2"

	if not os.path.exists(command_line) :
		print >> sys.stderr , "[ERROR] wrong or no path to minimap2"
		sys.exit(1)

	command_line+=" --cs -t " + str(cores) + " " + parameters + " " + ref_file + " " + " ".join(str(x) for x in query_file) + filter
	print >> sys.stderr, "### Running command line: " + command_line
	out_file = open(out_file_name , 'w')
	mapProcess = subprocess.Popen( command_line , shell=True, stdout=out_file)
	output, error = mapProcess.communicate()
	out_file.close()

	return out_file_name


def read_psl( psl_file ) :
	original_hits = {}
	id=0
	for line in open( psl_file ) :
		id+=1
		# subset only the necessary columns 1 -> 17
		matches , mismatch , rep , Ns , Qgap_num , Qgap_len , Tgap_num , Tgap_len , strand , Qid , Qlen , Qstart , Qstop , Tid , Tlen , Tstart , Tstop = line.rstrip().split("\t")[0:17]
		hitLen = int(matches) + int(mismatch) + int(Tgap_len)
		if ( Tid , Qid ) not in original_hits :
			# create a new entry
			original_hits[(Tid,Qid)] = []
		# Add a new range to the db
		original_hits[(Tid,Qid)].append([ id , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ])
	return original_hits


def map_blat(ref_file , query_file , strand_filter = "" , paramters = "" , out_file_name = "tmp" , path="" ) :
	if path == "" :
		blat_search=subprocess.Popen( "which blat" , shell=True, stdout=subprocess.PIPE )
		command_line , error = blat_search.communicate()
		command_line = command_line.rstrip()
		if command_line == "" :
			print >> sys.stderr , '[ERROR] Blat expected to be in $PATH, not found'
			exit(1)
	else :
		command_line = path + "/blat"

	if strand_filter == "+" or strand_filter == "-" :
		map_file = out_file_name + ".raw"
		blat_err_file = open(map_file + ".blat.err"  , "w")
		mappingCommand = command_line + " -noHead " + str(paramters) + " " + ref_file + " " + query_file + " " + map_file
		print >> sys.stderr, "#### Mapping command line: " + mappingCommand
		mapProcess = subprocess.Popen(mappingCommand, shell=True , stdout=blat_err_file )
		output, error = mapProcess.communicate()
		psl_file = open(out_file_name , "w")
		command_line = "awk \'$9==\"" + strand_filter + "\"\' " + map_file
		print >> sys.stderr, "#### Filtering command line: " + command_line
		filterProcess = subprocess.Popen(command_line, shell=True, stdout=psl_file  )
		output, error = filterProcess.communicate()
		psl_file.close()
		blat_err_file.close()

	else :
		map_file = out_file_name
		blat_err_file = open(map_file + ".blat.err"  , "w")
		mappingCommand = command_line + " -noHead " + str(paramters) + " " + ref_file + " " + query_file + " " + map_file
		print >> sys.stderr, "#### Mapping command line: " + mappingCommand
		mapProcess = subprocess.Popen(mappingCommand, shell=True , stdout=blat_err_file )
		output, error = mapProcess.communicate()
		blat_err_file.close()

	return out_file_name


def read_nucmer_coords( coords_file ) :
	# coords file as produced with show-coords + " -l -r -T -H "
	original_hits = {}
	id=0
	for line in open(coords_file) :
		id+=1
		if line[0] == "#" or line.rstrip() == "": continue ;
		if line[0] == "" or line.rstrip() == "": continue ;
		Tstart , Tstop , Qstart , Qstop , hitLen , QhitLen , iden , Tlen , Qlen , Tcov , Qcov , Tid , Qid = line.rstrip().split("\t")[0:13]
		matches = int( float(hitLen) * float(iden)/100 )

		if ( Tid , Qid ) not in original_hits :
			# create a new entry
			original_hits[(Tid,Qid)] = []

		# Add a new range to the db
		original_hits[(Tid,Qid)].append([ id , int(Tstart) , int(Tstop) , int(Qstart) , int(Qstop) , int(matches) , int(hitLen) ])
	return original_hits


def map_nucmer_unplaced_on_target(target_fasta, name , cores , chunk_db, workdir):
	unplaced_on_target_coords = workdir + "/" + name #map.unplaced.coords
	unplaced_file = workdir + "/tmp.unplaced.for_and_rev"
	unplaced_len_file = workdir + "/tmp.unplaced.len.json"
	if not ( os.path.exists(unplaced_file) ) :
		unplaced_len = {}
		print >> sys.stderr, "#### Temporary unplaced sequences FASTA file: " + unplaced_file
		print >> sys.stderr, "#### Temporary unplaced sequences length database: " + unplaced_len_file
		#make fasta files
		unplaced = open(unplaced_file, 'w')
		unplaced_fasta = read_fasta( chunk_db["inputs"]["U"] )
		for seq in sorted( unplaced_fasta.keys() ) :
			print >> unplaced , ">" + seq + "|+"
			print >> unplaced , str(Seq(unplaced_fasta[seq])).upper()
			print >> unplaced , ">" + seq + "|-"
			print >> unplaced , str(Seq(unplaced_fasta[seq]).reverse_complement()).upper()
			unplaced_len[ seq+ "|-" ] = len(unplaced_fasta[seq])
			unplaced_len[ seq+ "|+" ] = len(unplaced_fasta[seq])
		unplaced.close()
		json.dump(unplaced_len , open(unplaced_len_file, 'w') , indent=4 )
	else :
		unplaced_len = json.load(open(unplaced_len_file))

	return map_nucmer( target_fasta , unplaced_file , int(cores) , unplaced_on_target_coords , chunk_db["tool_paths"]["nucmer"] , chunk_db["tool_paths"]["show-coords"] , " --forward " , " -l -r -T -H ") , unplaced_len


def map_nucmer( ref_file , query_file ,  cores ,  out_file_name , nucmer_path , showcoords_path , parameters = "" , filter = "" ) :
	out_file_name_prefix = out_file_name.split(".coords")[0]

	if nucmer_path == "" :
		nucmer_search=subprocess.Popen( "which nucmer" , shell=True, stdout=subprocess.PIPE )
		nucmer_command_line , error = nucmer_search.communicate()
		nucmer_command_line = nucmer_command_line.rstrip()
		if nucmer_command_line == "" :
			print >> sys.stderr , '[ERROR] Nucmer expected to be in $PATH, not found'
			exit(1)
	else :
		nucmer_command_line = nucmer_path + "/nucmer"
	nucmer_command_line += " -p " + out_file_name_prefix + " -t " + str(cores) + " " + parameters + " "
	map_file_err = open( out_file_name_prefix + ".err" , "w")
	print >> sys.stderr, "### Running command line: " + nucmer_command_line + ref_file + " " + query_file + " "
	mapProcess = subprocess.Popen(nucmer_command_line + ref_file + " " + query_file + " ", shell=True, stderr=map_file_err)
	output, error = mapProcess.communicate()
	map_file_err.close()

	if showcoords_path == "" :
		showcoords_search=subprocess.Popen( "which show-coords" , shell=True, stdout=subprocess.PIPE )
		nucmer_command_line , error = showcoords_search.communicate()
		extract_coords_process = nucmer_command_line.rstrip()
		if extract_coords_process == "" :
			print >> sys.stderr , '[ERROR] show-coords expected to be in $PATH, not found'
			exit(1)
	else :
		extract_coords_process = showcoords_path + "/show-coords"

	extract_coords_process += " -c " + filter + " "
	coords_file = open( out_file_name , 'w' )
	input_delta = out_file_name_prefix + ".delta"
	print >> sys.stderr, "### Running command line: " + extract_coords_process + input_delta
	coordsProcess = subprocess.Popen(extract_coords_process + input_delta, shell=True, stdout=coords_file)
	output, error = coordsProcess.communicate()
	coords_file.close()

	return out_file_name


def index_bam( bam_file , samtools_path) :
	#print >> sys.stderr, "#### Indexing BAM file: " + bam_file
	if samtools_path == "" :
		samtools_search=subprocess.Popen( "which samtools" , shell=True, stdout=subprocess.PIPE )
		samtools_command , error = samtools_search.communicate()
		samtools_command = samtools_command.rstrip()
		if samtools_command == "" :
			print >> sys.stderr , '[ERROR] Samtools expected to be in $PATH, not found'
			exit(1)
	else :
		samtools_command = samtools_path + "/samtools"

	if not os.path.exists(samtools_command) :
		print >> sys.stderr , "[ERROR] Wrong or no path to Samtools"
		sys.exit(1)

	command_line = samtools_command + " index " + bam_file
	print >> sys.stderr, "##### Running command line: " + command_line
	index_process = subprocess.Popen( command_line , shell=True, stdout=subprocess.PIPE)
	output, error = index_process.communicate()


def sam2sorted_bam( sam_file , bam_file ,  ref_file , cores , path , clean=True) :

	if path == "" :
		samtools_search=subprocess.Popen( "which samtools" , shell=True, stdout=subprocess.PIPE )
		samtools_command , error = samtools_search.communicate()
		samtools_command = samtools_command.rstrip()
		if samtools_command == "" :
			print >> sys.stderr , '[ERROR] Samtools expected to be in $PATH, not found'
			exit(1)
	else :
		samtools_command = path + "/samtools"

	if not os.path.exists(samtools_command) :
		print >> sys.stderr , "[ERROR] Wrong or no path to Samtools"
		sys.exit(1)

	# sort
	command_line = samtools_command +" sort -l 9 -@ " + str(cores) +  " -m 500M -o " + bam_file + " " + sam_file
	print >> sys.stderr, "##### Running command line: " + command_line
	sort_process = subprocess.Popen( command_line, shell=True , stdout=subprocess.PIPE )
	output, error = sort_process.communicate()
	# Index
	index_bam( bam_file , path)

	if clean :
		os.remove(sam_file)

	return bam_file


def paf2pair( unique_paf_file ,  seq_id_1 , seq_id_2 , chunk_db) :
	pairing_file_1 = chunk_db["sequences"][seq_id_1]["folder"] + "/" + seq_id_1 + ".pairing.txt.gz"
	pairing_file_2 = chunk_db["sequences"][seq_id_2]["folder"] + "/" + seq_id_2 + ".pairing.txt.gz"

	f = open( pairing_file_1 , 'w' )
	g = open( pairing_file_2 , 'w' )
	for line in open(unique_paf_file) :
		if line[0] == "#" or line.rstrip() == "": continue ;
		if line[0] == "" or line.rstrip() == "": continue ;
		Qid , Qlen , Qstart , Qstop , strand , Tid , Tlen , Tstart , Tstop , matches , hitLen = line.rstrip().split("\t")[0:11]
		if Qid == seq_id_1 and Tid == seq_id_2 :
			pairing_line_1 = "\t".join( [ Qid, Qstart , Qstop , Tid , Tstart , Tstop , matches , hitLen ] )
			pairing_line_2 = "\t".join( [ Tid , Tstart , Tstop , Qid, Qstart , Qstop , matches , hitLen ] )
		elif Tid == seq_id_1 and Qid == seq_id_2 :
			pairing_line_1 = "\t".join( [ Tid , Tstart , Tstop , Qid, Qstart , Qstop , matches , hitLen ] )
			pairing_line_2 = "\t".join( [ Qid, Qstart , Qstop , Tid , Tstart , Tstop , matches , hitLen ] )
		else :
			print >> sys.stderr, "#### Unexpected sequence name in alignment file"
			exit(6)
		print >> f, pairing_line_1
		print >> g, pairing_line_2
	f.close()
	g.close()

	return pairing_file_1 , pairing_file_2


def map_minimap_dotplot( ref_prefix , ref_file_name , query_prefix , query_file_name , out_folder , cores , paths,  forward_only = True) :
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

	if forward_only :
		out_file_name_prefix = query_prefix + ".on." + ref_prefix + ".forward"
		mappingCommand = command_line + " --cs -x asm20 -r 1000 --for-only -t " + str(cores) + " "
	else :
		out_file_name_prefix = query_prefix + ".on." + ref_prefix
		mappingCommand = command_line + " --cs -x asm20 -r 1000 -t " + str(cores) + " "

	map_file = open( out_folder + "/" + out_file_name_prefix + ".paf" , "w")
	print >> sys.stderr, "##### Running command line: " + mappingCommand + ref_file_name + " " + query_file_name + " "
	mapProcess = subprocess.Popen(mappingCommand + ref_file_name + " " + query_file_name + " ", shell=True, stdout=map_file)
	output, error = mapProcess.communicate()
	map_file.close()
	return out_file_name_prefix


def map_nucmer_dotplot( ref_prefix , ref_file_name , query_prefix , query_file_name , out_folder , cores , paths , forward_only = False ) :
	nucmer_path = paths["nucmer"]

	if nucmer_path == "" :
		minimap2_search=subprocess.Popen( "which nucmer" , shell=True, stdout=subprocess.PIPE )
		command_line , error = minimap2_search.communicate()
		command_line = command_line.rstrip()
		if command_line == "" :
			print >> sys.stderr , '[ERROR] Nucmer expected to be in $PATH, not found'
			exit(1)
	else :
		command_line = nucmer_path + "/nucmer"

	if not os.path.exists(command_line) :
		print >> sys.stderr , "[ERROR] Wrong or no path to nucmer (Mummer4)"
		sys.exit(1)

	if forward_only :
		out_file_name_prefix = query_prefix + ".on." + ref_prefix + ".forward"
		mappingCommand = "nucmer --forward -p " + out_folder + "/" + out_file_name_prefix + " -t " + str(cores) + " "
	else :
		out_file_name_prefix = query_prefix + ".on." + ref_prefix
		mappingCommand = "nucmer -p " + out_folder + "/" + out_file_name_prefix + " -t " + str(cores) + " "

	map_file = open( out_folder + "/" + out_file_name_prefix + ".err" , "w")
	print >> sys.stderr, "##### Running command line: " + mappingCommand + ref_file_name + " " + query_file_name + " "
	mapProcess = subprocess.Popen(mappingCommand + ref_file_name + " " + query_file_name + " ", shell=True, stderr=map_file)
	output, error = mapProcess.communicate()
	map_file.close()

	showcoords_path = paths["show-coords"]

	if showcoords_path == "" :
		minimap2_search=subprocess.Popen( "which show-coords" , shell=True, stdout=subprocess.PIPE )
		command_line , error = minimap2_search.communicate()
		command_line = command_line.rstrip()
	else :
		command_line = showcoords_path + "/show-coords"

	if not os.path.exists(command_line) :
		print >> sys.stderr , "[ERROR] Wrong or no path to show-coords (Mummer4)"
		sys.exit(1)

	extract_coords_process = command_line + " -c "
	coords_file = open( out_folder + "/" + out_file_name_prefix + ".coords" , 'w' )
	input_delta = out_folder + "/" + out_file_name_prefix + ".delta"
	print >> sys.stderr, "##### Running command line: " + extract_coords_process + input_delta
	coordsProcess = subprocess.Popen(extract_coords_process + input_delta, shell=True, stdout=coords_file)
	output, error = coordsProcess.communicate()
	coords_file.close()

	return out_file_name_prefix


def map_regions( left_seq , right_seq , mapper , cores , tempdir , paths ) :
	map_results = []
	left_seq_file_name = tempdir + "/left.fasta"
	#print >> sys.stderr , left_seq_file_name
	left_seq_file = open(left_seq_file_name , 'w')
	print >> left_seq_file , ">left"
	print >> left_seq_file , left_seq
	left_seq_file.close()
	right_seq_file_name = tempdir + "/right.fasta"
	#print >> sys.stderr , right_seq_file_name
	right_seq_file = open(right_seq_file_name , 'w')
	print >> right_seq_file , ">right"
	print >> right_seq_file , right_seq
	right_seq_file.close()

	if mapper == "blat" :
		blat_path = paths["blat"]
		psl_file = tempdir + "/alignment.psl"
		psl_file = map_blat( left_seq_file_name , right_seq_file_name , "+" , " -minIdentity=80 " , psl_file , blat_path )
		map_results = read_psl( psl_file )
	#elif mapper == "nucmer" :
	#	nucmer_path = paths["nucmer"]
	#	showcoords_path = paths["show-coords"]
	#	coords_file = tempdir + "/alignment.coords"
	#	coords_file = map_nucmer( left_seq_file_name , right_seq_file_name ,  int(cores) ,  coords_file , nucmer_path , showcoords_path , " --forward " , " -l -r -T -H ")
	#	map_results = read_nucmer_coords( coords_file )
	#elif mapper == "minimap" :
	#	# Minimap hits never include extremities (0 or END) of any sequence
	#	minimap_path = paths["minimap2"]
	#	paf_file = tempdir + "/alignment.paf"
	#	paf_file = map_minimap( left_seq_file_name , right_seq_file_name , int(cores) , " -x asm20 --for-only " , paf_file , minimap_path )
	#	map_results = read_paf( paf_file )
	else :
		print >> sys.stderr, "[ERROR] Unknown mapping tool (" + mapper + "). Valid values for \"-m\"|\"--mapper\" flags are \"minimap\" or \"nucmer\""

	return map_results


def find_extremity_overlap( hits_db , left_db , right_db ) :
	# left_db/right_db format:
	# right_block["region_given"] = [seqID , int(start) , int(stop) , strand]
	# right_block["region_corrected"] = [seqID , int(corr_start) , int(corr_stop) , strand]
	# right_block["seq"] = seq
	# right_block["annot_on_fasta"] = [ [ gene_start , gene_stop , gene_id ] , [...] , ... ]
	#
	# hits_db format:
	# hits_db[ ( Target_id , Query_id ) ]
	# hits_db[("left","right")] = [
	#	#	...
	#	#	[ id , int(left_start) , int(left_stop) , int(right_start) , int(right_stop) , int(matches) , int(hitLen) ])
	#	#	...
	#	]
	if hits_db == {} :
		overlap_region = {}
		return overlap_region
	else :
		left_len = len(left_db["seq"])
		right_len = len(right_db["seq"])
		# Sort hits by start on right
		sorted_hits = sorted( hits_db[("left","right")] , key=lambda x: int(x[3]) , reverse=True)
		keep_searching = True
		overlap_region = {}
		overlap_length = 0
		# Good candiadate for overlap:
		#	max 500bp from start of right region
		#	max 500bp form end of left region
		#	length > 1Kb
		#	Do not take all the left or all the right region >> would mean delete it entirely >> not allowed

		# Scan all hits from the beginning of the right region
		while keep_searching  and ( len(sorted_hits) > 0 ) :
			candidate_hit = sorted_hits.pop(0)
			id , left_start , left_stop , right_start , right_stop , matches , range = candidate_hit
			if right_start > 500 :
				# This and the following hits are too far from the beginning fo right region to be good
				keep_searching = False
			else:
				if not left_stop < left_len - 500 :
					# hit is not too far from the end of left region
					if matches > 1000 :
						# hit is not too short
						# good candidate >> if valid save it in overlap_region and stop searching
						if ( left_start < 500 ) or ( right_stop > (right_len - 500) ) :
							# The overlap does cover almost completely the left or the right sequence >> could obliterate it >> not allowed
							print >> sys.stderr, "[WARNING] Overlap found that covers an entire sequence. Discarded mapping to avoid obliteration"
							keep_searching = False
						else :
							if overlap_region == {} or overlap_length < matches :
								overlap_region = [left_start , left_stop , right_start , right_stop]
								overlap_length = matches
								print >> sys.stderr, "#### Overlap found between blocks: " + str(overlap_region)

		# overlap_region =
		# 	if overlap found: [left_start , left_stop , right_start , right_stop]
		#	else : {}
		return overlap_region


def map_sequences( left_id , left_seq , right_id , right_seq , mapper , cores , tempdir , prefix ,  paths ) :
		map_results = []
		left_seq_file_name = tempdir + "/" + left_id + ".fasta"
		left_seq_file = open(left_seq_file_name)
		print >> left_seq_file , ">" + left_id
		print >> left_seq_file , left_seq
		right_seq_file_name = tempdir + "/" + right_id + ".fasta"
		right_seq_file = open(right_seq_file_name)
		print >> right_seq_file , ">" + right_id
		print >> right_seq_file , right_seq
		left_seq_file.close()
		right_seq_file.close()

		if mapper == "blat" :
			blat_path = paths["blat"]
			psl_file = tempdir + "/" + prefix + ".psl"
			psl_file = map_blat( left_seq_file , right_seq_file , "+" , " -minIdentity=80 " , psl_file , blat_path )
			map_results = read_psl( psl_file )
		elif mapper == "minimap" :
			# Minimap hits never include extremities (0 or END) of any sequence
			minimap_path = paths["minimap2"]
			paf_file = tempdir + "/" + prefix + ".paf"
			paf_file = map_minimap( left_seq_file , right_seq_file , int(cores) , " -x asm20 --for-only " , paf_file , minimap_path )
			map_results = read_paf( paf_file )
		elif mapper == "nucmer" :
			nucmer_path = paths["nucmer"]
			showcoords_path = paths["show-coords"]
			coords_file = tempdir + "/" + prefix + ".coords"
			coords_file = map_nucmer( left_seq_file , right_seq_file ,  int(cores) ,  coords_file , nucmer_path , showcoords_path , " --forward " , " -l -r -T -H ")
			map_results = read_nucmer_coords( coords_file )
		else :
			print >> sys.stderr, "[ERROR] Unknown mapping tool (" + mapper + "). Valid values for \"-m\"|\"--mapper\" flags are \"minimap\" or \"nucmer\""

		return map_results


def do_count_hits_Hap(ref_list, hit1, hit2, ref , annotation_dict ) :
	count_db = {}

	for element in sorted(ref_list):
		try :
			chr, start , end , id = element
		except :
			print >> sys.stderr, element
			exit(1)

		if id in annotation_dict :
			feats = annotation_dict[id]
			descriptions = ""
			counts = ""

			for desc in sorted(feats.keys()) :
				if not descriptions == "" :
					descriptions += ";"
					counts += ";"
				descriptions += desc
				counts += str(feats[desc])
		else :
			descriptions = "-"
			counts = "-"

		if ref == "_Hap1_" :
			chr_hap1 = chr
			chr_hap2 = chr.replace("_Hap1_", "_Hap2_")

			if chr_hap1 not in count_db :
				count_db[chr_hap1] = []
			try :
				h1_len = len(hit1[chr_hap1][id])
			except :
				h1_len = 0
			try :
				h2_len = len(hit2[chr_hap2][id])
			except :
				h2_len = 0

			if h1_len == 0 :
				ratio = "inf"
			else :
				ratio = float(h2_len) / float(h1_len)
			count_db[chr_hap1].append([chr_hap1, start , end , id , h1_len , h2_len , ratio , descriptions , counts])

		elif ref == "_Hap2_" :
			chr_hap1 = chr.replace("_Hap2_", "_Hap1_")
			chr_hap2 = chr

			if chr_hap2 not in count_db :
				count_db[chr_hap2] = []
			try :
				h1_len = len(hit1[chr_hap1][id])
			except :
				h1_len = 0
			try :
				h2_len = len(hit2[chr_hap2][id])
			except :
				h2_len = 0

			if h1_len == 0 :
				ratio = "inf"
			else :
				ratio = float(h2_len) / float(h1_len)
			count_db[chr_hap2].append([chr_hap2, start , end , id , h1_len , h2_len , ratio , descriptions , counts])

		else :
			print >> sys.stderr , "[ERROR] Unexpected reference haplotype"
			sys.exit(1)

	return count_db


def print_hit_counts(hit_db, outfile_name):
	outfile = open(outfile_name, 'w')
	print >> outfile , "\t".join(["Chr","Start","Stop","Gene_Id","Hap1_count","Hap2_count","Hap2_to_Hap1_ratio","Description","Gene_count_with_description"])
	for chr in sorted(hit_db.keys()) :
		for hit in sorted(hit_db[chr]) :
			print >> outfile , "\t".join([str(x) for x in hit])

	outfile.close()

	return outfile_name


def do_count_hits( ref_list, hit1, hit2, ref , pairs_list , annotation_dict ) :
	# annotation_dict:
	# annotation_dict[feat_name][feat_description] = feat_count
	count_db = {}
	pair_db = {}
	if ref == 1 :
		for pair in pairs_list:
			pair_db[pair[0]] = pair[1]
	elif ref == 2 :
		for pair in pairs_list:
			pair_db[pair[1]] = pair[0]
	else :
		print >> sys.stderr , "[ERROR] Unexpected reference haplotype"
		sys.exit(1)
	for element in sorted(ref_list):
		try :
			chr, start , end , id = element
		except :
			print >> sys.stderr, "[ERROR] Error in annotation for: " + str(element)
			exit(1)
		else:
			if id in annotation_dict :
				feats = annotation_dict[id]
				descriptions = ""
				counts = ""

				for desc in sorted(feats.keys()) :
					if not descriptions == "" :
						descriptions += ";"
						counts += ";"
					descriptions += desc
					counts += str(feats[desc])
			else :
				descriptions = "-"
				counts = "-"

			if ref == 1 :
				# Use first element in pairs == Hap1
				chr_hap1 = chr
				chr_hap2 = pair_db[chr]
				if chr_hap1 not in count_db :
					count_db[chr_hap1] = []
				try :
					h1_len = len(hit1[chr_hap1][id])
				except :
					h1_len = 0
				try :
					h2_len = len(hit2[chr_hap2][id])
				except :
					h2_len = 0

				if h1_len == 0 :
					ratio = "inf"
				else :
					ratio = float(h2_len) / float(h1_len)
				count_db[chr_hap1].append([chr_hap1, start , end , id , h1_len , h2_len , ratio , descriptions , counts])
			elif ref == 2 :
				# Use second element in pairs == Hap2
				chr_hap1 = pair_db[chr]
				chr_hap2 = chr
				if chr_hap2 not in count_db :
					count_db[chr_hap2] = []
				try :
					h1_len = len(hit1[chr_hap1][id])
				except :
					h1_len = 0
				try :
					h2_len = len(hit2[chr_hap2][id])
				except :
					h2_len = 0

				if h1_len == 0 :
					ratio = "inf"
				else :
					ratio = float(h2_len) / float(h1_len)
				count_db[chr_hap2].append([chr_hap2, start , end , id , h1_len , h2_len , ratio , descriptions , counts])
	return count_db


def read_gmap_results_paired_sequences(gff3, threshold_cov, threshold_iden, group_by, pairs_list, mRNA_to_gene =""):
	hap1_sequences = []
	hap2_sequences = []
	pairs_db = {}
	for pair in pairs_list:
		try :
			hap1 , hap2 = pair
		except ValueError :
			print >> sys.stderr , "[ERROR] Unexpected number of sequences in pair: " + pair
			sys.exit(1)
		else :
			hap1_sequences.append(hap1)
			hap2_sequences.append(hap2)
			pairs_db[hap1] = hap2
			pairs_db[hap2] = hap1

	hit_table = {}
	# Filter hits
	for line in open(gff3, 'r') :
		if line.rstrip() == "" or line[0] == "#" :
			continue
		seqname, source, feature, start, end, score, strand, frame, attribute = line.rstrip().split("\t")
		if not feature == "mRNA" :
			# Skip all lines not of features different from mRNA
			continue
		else :
			if seqname not in hit_table :
				hit_table[seqname] = defaultdict(list)

			attributes_dict = {}
			for chunk in attribute.rstrip(";").split(";"):
				att = chunk.split("=")
				attributes_dict[att[0]] = att[1]

			if "coverage" not in attributes_dict :
				attribute = str(attribute).replace( "coverage" , ";coverage")
				for chunk in attribute.rstrip(";").split(";"):
					att = chunk.split("=")
					attributes_dict[att[0]] = att[1]

			if "identity" not in attributes_dict :
				attribute = str(attribute).replace( "identity" , ";identity")
				for chunk in attribute.rstrip(";").split(";"):
					att = chunk.split("=")
					attributes_dict[att[0]] = att[1]

			if ( "coverage" not in attributes_dict ) or ( "identity" not in attributes_dict ) :
				print >> sys.stdout , "[ERROR] Unexpected syntax for GMAP output, coverage or identity fields missing "
				print >> sys.stderr , "[ERROR] Unexpected syntax for GMAP output, coverage or identity fields missing "
				print >> sys.stdout , "[ERROR] Line: " + line.rstrip()
				print >> sys.stdout , "[ERROR] Attributes: " + attribute.rstrip()
				sys.exit(1)

			if float(attributes_dict["coverage"] > int(threshold_cov) ) and float(attributes_dict["identity"] > int(threshold_iden) ):
				mRNA_id = attributes_dict["Name"]
				if group_by == "by_mRNA" :
					element_id = mRNA_id
				else :
					element_id = mRNA_to_gene[mRNA_id]

				if element_id not in hit_table[seqname] :
					hit_table[seqname][element_id] = []

				hit_table[seqname][element_id].append([int(start) , int(end)])

	# Collapse overlapping calls and split by haplotype
	hit_table_hap1 = defaultdict(list)
	hit_table_hap2 = defaultdict(list)
	for chr in sorted(hit_table.keys()) :
		clean_hits = {}
		for locus in sorted( hit_table[chr].keys() ):
			clean_hits[locus] = []
			stop = ""
			for hit in sorted(hit_table[chr][locus]) :
				if stop == "" :
					clean_hits[locus].append(hit)
					stop = hit[1]
				else :
					if hit[0] > stop :
						clean_hits[locus].append(hit)
						stop = hit[1]

		if chr in hap1_sequences :
			hit_table_hap1[chr] = clean_hits
		elif chr in hap2_sequences :
			hit_table_hap2[chr] = clean_hits
		else :
			print >> sys.stderr, "[WARNING] Hits mapping on unpaired sequence: " + str(chr)
	return hit_table_hap1 , hit_table_hap2


def read_gmap_results_Hap(gff3, threshold_cov , threshold_iden , group_by, mRNA_to_gene =""):
	hit_table = {}
	hit_table_hap1 = defaultdict(list)
	hit_table_hap2 = defaultdict(list)

	# Extract good hits
	for line in open(gff3, 'r') :
		if line.rstrip() == "" or line[0] == "#" :
			continue

		seqname, source, feature, start, end, score, strand, frame, attribute = line.rstrip().split("\t")

		if seqname not in hit_table :
			hit_table[seqname] = defaultdict(list)

		if not feature == "mRNA" :
			# Skip all lines not of features different from mRNA
			continue

		attributes_dict = {}
		for chunk in attribute.rstrip(";").split(";"):
			att = chunk.split("=")
			attributes_dict[att[0]] = att[1]

		if "coverage" not in attributes_dict :
			attribute = str(attribute).replace( "coverage" , ";coverage")
			for chunk in attribute.rstrip(";").split(";"):
				att = chunk.split("=")
				attributes_dict[att[0]] = att[1]

		if "identity" not in attributes_dict :
			attribute = str(attribute).replace( "identity" , ";identity")
			for chunk in attribute.rstrip(";").split(";"):
				att = chunk.split("=")
				attributes_dict[att[0]] = att[1]

		if ( "coverage" not in attributes_dict ) or ( "identity" not in attributes_dict ) :
			print >> sys.stdout , "[ERROR] Unexpected syntax for GMAP output, coverage or identity fields missing "
			print >> sys.stderr , "[ERROR] Unexpected syntax for GMAP output, coverage or identity fields missing "
			print >> sys.stdout , "[ERROR] Line: " + line.rstrip()
			print >> sys.stdout , "[ERROR] Attributes: " + attribute.rstrip()
			sys.exit(1)

		if float(attributes_dict["coverage"] > int(threshold_cov) ) and float(attributes_dict["identity"] > int(threshold_iden) ):
			mRNA_id = attributes_dict["Name"]
			if group_by == "by_mRNA" :
				element_id = mRNA_id
			else :
				element_id = mRNA_to_gene[mRNA_id]

			if element_id not in hit_table[seqname] :
				hit_table[seqname][element_id] = []

			hit_table[seqname][element_id].append([int(start) , int(end)])

	# Collapse overlapping calls and split by haplotype
	for chr in sorted(hit_table.keys()) :
		clean_hits = {}
		for locus in sorted( hit_table[chr].keys() ):
			clean_hits[locus] = []
			stop = ""
			for hit in sorted(hit_table[chr][locus]) :
				if stop == "" :
					clean_hits[locus].append(hit)
					stop = hit[1]
				else :
					if hit[0] > stop :
						clean_hits[locus].append(hit)
						stop = hit[1]

		if "_Hap1_" in chr :
			hit_table_hap1[chr] = clean_hits
		elif "_Hap2_" in chr :
			hit_table_hap2[chr] = clean_hits

	return hit_table_hap1 , hit_table_hap2


def read_gmap_results(annot_gff3_db, map_file , threshold_cov , threshold_iden , group_by):
	#### annot_gff3_db[gene_id] = [ gff_line , chr , int(start) , mRNA_dict{} ]
	####		mRNA_dict[mRNA_id] = [ gff_line , int(start), feat_dict{} ]
	####			feat_dict[mRNA_order_num] = [ [ line , int(start) , int(end) ]  ]
	####		 		with mRNA_order_num meaning:	1 = "exon"
	####												2 = "five_prime_UTR"
	####												3 = "CDS"
	####												4 = "three_prime_UTR"
	#### line (unsplitted on tab) = seqname, source, feature, start, end, score, strand, frame, attribute
	hit_table = {}
	feat_position = {}
	if group_by == "gene" :
		#Make mRNA ID to gene ID correspondence db
		mRNA_to_gene = {}
		for gene_id in annot_gff3_db.keys() :
			gff_line = annot_gff3_db[gene_id][0].split("\t")
			feat_chr = gff_line[0]
			feat_start = gff_line[3]
			feat_stop = gff_line[4]
			feat_position[gene_id] = [feat_chr , feat_start, feat_stop]
			for mrna_id in annot_gff3_db[gene_id][3].keys() :
				mRNA_to_gene[mrna_id] = gene_id
	else :
		for gene_id in annot_gff3_db.keys() :
			for mrna_id in annot_gff3_db[gene_id][3].keys() :
				gff_line = annot_gff3_db[gene_id][3][mrna_id][0].split("\t")
				feat_chr = gff_line[0]
				feat_start = gff_line[3]
				feat_stop = gff_line[4]
				feat_position[mrna_id] = [feat_chr , feat_start, feat_stop]

	# Extract good hits
	for line in open(map_file, 'r') :
		if line.rstrip() == "" or line[0] == "#" :
			continue
		seqname, source, feature, start, end, score, strand, frame, attribute = line.rstrip().split("\t")
		if seqname not in hit_table :
			hit_table[seqname] = defaultdict(list)
		if not feature == "mRNA" :
			# Skip all lines not of features different from mRNA
			continue
		attributes_dict = {}
		for chunk in attribute.rstrip(";").split(";"):
			att = chunk.split("=")
			attributes_dict[att[0]] = att[1]
		if "coverage" not in attributes_dict :
			attribute = str(attribute).replace( "coverage" , ";coverage")
			for chunk in attribute.rstrip(";").split(";"):
				att = chunk.split("=")
				attributes_dict[att[0]] = att[1]
		if "identity" not in attributes_dict :
			attribute = str(attribute).replace( "identity" , ";identity")
			for chunk in attribute.rstrip(";").split(";"):
				att = chunk.split("=")
				attributes_dict[att[0]] = att[1]
		if ( "coverage" not in attributes_dict ) or ( "identity" not in attributes_dict ) :
			print >> sys.stdout , "[ERROR] Unexpected syntax for GMAP output, coverage or identity fields missing "
			print >> sys.stderr , "[ERROR] Unexpected syntax for GMAP output, coverage or identity fields missing "
			print >> sys.stdout , "[ERROR] Line: " + line.rstrip()
			print >> sys.stdout , "[ERROR] Attributes: " + attribute.rstrip()
			sys.exit(1)

		if float(attributes_dict["coverage"] > int(threshold_cov) ) and float(attributes_dict["identity"] > int(threshold_iden) ):
			mRNA_id = attributes_dict["Name"]
			if group_by == "mRNA" :
				element_id = mRNA_id
			elif group_by == "gene" :
				element_id = mRNA_to_gene[mRNA_id]

			if element_id not in hit_table[seqname] :
				hit_table[seqname][element_id] = []
			hit_table[seqname][element_id].append([int(start) , int(end)])

	# remove overlapping calls of same feature
	clean_hits = {}
	for chr in sorted(hit_table.keys()) :
		clean_hits[chr] = {}
		for locus in sorted( hit_table[chr].keys() ):
			clean_hits[chr][locus] = []
			stop = ""
			for hit in sorted(hit_table[chr][locus]) :
				if stop == "" :
					clean_hits[chr][locus].append(hit)
					stop = hit[1]
				else :
					if hit[0] > stop :
						clean_hits[chr][locus].append(hit)
						stop = hit[1]
	# do counts
	filtered_hit_table = []
	for chr in sorted(clean_hits.keys()) :
		for locus in sorted( clean_hits[chr].keys() ):
			feat_chr , feat_start, feat_stop = feat_position[locus]
			hit_count = len(clean_hits[chr][locus])
			filtered_hit_table.append( [feat_chr , feat_start, feat_stop , locus , hit_count]  )

	return filtered_hit_table


def index_gmap( ref_name , ref_id , out_prefix , out_dir , index_dir , paths = "" ) :
	try :
		gmap_path = paths["gmap"]
	except :
		print >> sys.stderr , '[ERROR] Path to gmap not set in configuration file'
		exit(1)
	else :
		if gmap_path == "" :
			gmap_buid_search=subprocess.Popen( "which gmap_build" , shell=True, stdout=subprocess.PIPE )
			gmap_buid_exec , error = gmap_buid_search.communicate()
			gmap_buid_exec = gmap_buid_exec.rstrip()
			if gmap_buid_exec == "" :
				print >> sys.stderr , '[ERROR] gmap_build expected to be in $PATH, not found'
				exit(1)
		else :
			gmap_buid_exec = gmap_path + "/gmap_buid"

	indexing_out_file = open( out_dir + "/" + out_prefix + ".gmap_index.log" ,"w" )
	indexing_err_file = open( out_dir + "/" + out_prefix + ".gmap_index.err" ,"w" )
	indexing_command = gmap_buid_exec + " -D " + index_dir + " -d " + ref_id + " " + ref_name
	indexProcess = subprocess.Popen(indexing_command, shell=True, stdout=indexing_out_file , stderr=indexing_err_file)
	output, error = indexProcess.communicate()
	indexing_out_file.close()
	indexing_err_file.close()
	return index_dir


def map_gmap( index_dir , ref_id , CDS_file_name , out_prefix , out_dir , cores = 4 , paths = "" ) :
	try :
		gmap_path = paths["gmap"]
	except :
		print >> sys.stderr , '[ERROR] Path to gmap not set in configuration file'
		exit(1)
	else :
		if gmap_path == "" :
			gmap_buid_search=subprocess.Popen( "which gmap" , shell=True, stdout=subprocess.PIPE )
			gmap_buid_exec , error = gmap_buid_search.communicate()
			gmap_exec = gmap_buid_exec.rstrip()
			if gmap_exec == "" :
				print >> sys.stderr , '[ERROR] Gmap expected to be in $PATH, not found'
				exit(1)
		else :
			gmap_exec = gmap_path + "/gmap"
	# Gmap CDS on results
	gmap_results = out_dir + "/" + out_prefix + ".on.ref.gmap.gff3"
	gmap_gff3 = open( gmap_results , "w" )
	gmap_err = open( gmap_results + ".err" , "w" )
	gmapCommand = gmap_exec + " -D " + index_dir + " -d " + ref_id + " -f 2 -n 500 -t " + str(cores) + " " + CDS_file_name
	print >> sys.stderr , '### gmap command: ' + gmapCommand
	gmapProcess = subprocess.Popen(gmapCommand, shell=True, stdout=gmap_gff3 , stderr=gmap_err)
	output, error = gmapProcess.communicate()
	gmap_err.close()
	gmap_gff3.close()
	return gmap_results
