import HTSeq
import itertools
import sys
import argparse
import collections
import csv

features = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
features_properties = dict()


def count_reads_by_position_in_features(bam_file,gtf_file,out_file,feature_type,id_attribute,multimapped_mode,stranded,minaqual,quiet):


	gff = HTSeq.GFF_Reader( gtf_file, end_included=True )
	k = 0

	# TODO: Check this to avoid the HACK
	# CIGAR match characters (including alignment match, sequence match, and
	# sequence mismatch
	# Ref https://drive5.com/usearch/manual/cigar.html
	com = ('M', '=', 'X')
	
	try:
		for f in gff:
			
			# feature type (3rd column in GFF file)
			if f.type == feature_type:
				try:
					feature_id = f.attr[id_attribute]
				except KeyError:
					raise ValueError(
							"Feature %s does not contain a '%s' attribute" %
							(f.name, id_attribute))
							
				
				features[ f.iv ] += feature_id
				features_properties[ f.name ] = f.iv.length,f.iv.start,f.iv.end
		

				if stranded != "no" and f.iv.strand == ".":
					raise ValueError(
							"Feature %s at %s does not have strand information" %
							(f.name, f.iv))

				k += 1
				if k % 100000 == 0 and not quiet:
					sys.stderr.write("%d GFF lines processed.\n" % k)
					sys.stderr.flush()
	except:
		sys.stderr.write("Error occured when processing GFF file (%s):\n" %	gff.get_line_number_string())
		raise
	
	if not k:
		sys.stderr.write("Nothing was found when processing GFF file.\nCheck the gff file or feature type")
		raise SystemExit
	 
	if not quiet:
		sys.stderr.write("%d GFF lines processed.\n" % k)
		sys.stderr.flush()

	# GET COUNTS
	bam = HTSeq.BAM_Reader(bam_file)
	counts = collections.Counter( )
	k = 0
	genes_set = set()

	# TODO remove this line
	for align in itertools.islice( bam, 100000):
	# ~ for align in bam:

		
		if k > 0 and k % 100000 == 0 and not quiet:
			sys.stderr.write("%d SAM alignment records processed.\n" % k)
			sys.stderr.flush()
		
		k += 1
		

		# TODO: Remove this
		# ~ for co in align.cigar:
			# ~ # cut long cigar strings, only the best of reads, maximum one variation
			# ~ if co.type not in com and co.size > 0:
				# ~ print("NO",co,co.type,co.size,co.ref_iv)
				# ~ continue
			# ~ else:
				# ~ print("SI",co,co.type,co.size,co.ref_iv)

					# ~ # print ("g_ids",g_ids,iv_f.start)
		
		
		if not align.aligned:
			counts[ "_unmapped" ] += 1
			continue

		# TODO: Remove Hack			
		if align.iv.length > 51:
			# counts[ "_too_long" ] += 1
			continue
		
		if align.aQual < minaqual:
			counts[ "_too_low_aQual" ] += 1
			continue

		gene_ids = set()
		
		# interval, gen
		for iv_f, gene in features[ align.iv ].steps():
			# Add element to the set.
			gene_ids |= gene

		if len(gene_ids) == 1:
			gene_id = list(gene_ids)[0]
			counts[ gene_id, iv_f.start] += 1
			genes_set |= gene
		
		elif len(gene_ids) == 0:
			counts[ "_no_feature" ] += 1

		else:
			if multimapped_mode == 'none':
				counts[ "_ambiguous" ] += 1
					
			elif multimapped_mode == 'all':
				for g_ids in list(gene_ids):
					gene_id = g_ids
					counts[ gene_id, iv_f.start] += 1
					genes_set |= gene
			else:
				sys.exit("Illegal multimap mode.")
	
	# GET RESULTS
	# TODO: Remove HACK
	counts[ "_zzz" ] = 0

	# New dictionary sorted on gene_name and then on nt position on that gene
	sorted_counts = sorted(counts.items(),key = lambda i: (i[0][0], i[0][1]))

	last_gene_id = ''
	additional_data = []
	list_of_reads = []
	
	# ~ my_iter = iter(sorted_counts)
	
	# ~ for i in my_iter:
		# ~ print(i)
		# ~ next(my_iter)
		
	# ~ iterator = iter(sorted_counts)
	# ~ try:
		# ~ while True:
			# ~ item = next(iterator)
			# ~ print(item)
	# ~ except StopIteration:
		# ~ pass
	# ~ finally:
		# ~ del iterator


	for (index, gene_pos_count) in enumerate(sorted_counts[:-1]):
		if index < len(sorted_counts):
			current, next_ = gene_pos_count, sorted_counts[index + 1]
			
			# print("CN",current, next_)			
			# Gene ID and read position
			gene_id_pos = current[0]
			# Read counts for current position
			read_counts = current[1]
			
		# _unmapped _no_feature or _ambiguous
		if isinstance(gene_id_pos, str):
			additional_data.append((gene_id_pos,read_counts))
			
		else:
			gene_id = gene_id_pos[0] 
			position = gene_id_pos[1]
			next_gene_id = next_[0][0]
			
			gene_id_length = features_properties.get(gene_id)[0]
			gene_id_start = features_properties.get(gene_id)[1]
			gene_id_end = features_properties.get(gene_id)[2]
		
		# Check if there is something mapped
		try: 
			gene_id
		except:
			sys.stderr.write("No genes mapped")
			raise SystemExit

		# Print only the Gene name
		if (last_gene_id != gene_id) :
			i = gene_id_start
			
			reads = []
			reads.append(gene_id)
			last_gene_id = gene_id

		# Still working with the same gene
		if (gene_id == next_gene_id):

			while (i < position):
				# TODO: Extend with 0s
				# reads.extend(repeat(0, position - i))
				# i = position
				reads.append(0)
				i+=1
			if i == position:
				reads.append(read_counts)
				i+=1
				
		# Last position(s) of the gene
		else:
			if i == position:
				reads.append(read_counts)
				i+=1

			while (i < gene_id_end):
				# TODO: Extend with 0s
				# reads.extend(repeat(0, gene_id_end - i))
				reads.append(0)
				i+=1

			list_of_reads.append(reads)
		
	# Last position of the file
	# TODO put this one
	gene_id_pos = next_[0]
	read_counts = next_[1]
	
	
	# Write to file
	write_to_file(out_file,additional_data)
	write_to_file(out_file,list_of_reads)

def write_to_file (FILE_TO_WRITE,data_to_write):
	with open(FILE_TO_WRITE,"a") as f:
		wr = csv.writer(f,delimiter=';')
		wr.writerows(data_to_write)
		

def main():
	pa = argparse.ArgumentParser(
		usage="%(prog)s [options] alignment_file gff_file",
		description="This script takes one or more alignment files in SAM/BAM " +
		"format and a feature file in GFF format and calculates for each feature " +
		"the number of reads mapping to each position")
	
	# multiple bams
	# ~ pa.add_argument(
			# ~ "samfilenames", nargs='+', type=str,
			# ~ help="Path to the SAM/BAM files containing the mapped reads. " +
			# ~ "If '-' is selected, read from standard input")
	
	pa.add_argument(
			"samfilenames", type=str,
			help="Path to the SAM/BAM files containing the mapped reads. " +
			"If '-' is selected, read from standard input")

	pa.add_argument(
			"featuresfilename", type=str,
			help="Path to the GTF file containing the features")
	pa.add_argument(
			"-f", "--format", dest="samtype",
			choices=("sam", "bam"), default="sam",
			help="Type of <alignment_file> data, either 'sam' or 'bam' (default: sam)")

	pa.add_argument(
			"-a", "--minaqual", type=int, dest="minaqual",
			default=10,
			help="Skip all reads with MAPQ alignment quality lower than the given " +
			"minimum value (default: 10). MAPQ is the 5th column of a SAM/BAM " +
			"file and its usage depends on the software used to map the reads.")
	
	 # ~ pa.add_argument(
			# ~ "-o", "--samout", type=str, dest="samouts",
			# ~ action='append',
			# ~ default=[],
			# ~ help="Write out all SAM alignment records into " +
			# ~ "SAM files (one per input file needed), annotating each line " +
			# ~ "with its feature assignment (as an optional field with tag 'XF')")
	pa.add_argument(
			"-o", "--samout", type=str, dest="samouts",
			help="Write out all SAM alignment records into " +
			"SAM files (one per input file needed), annotating each line " +
			"with its feature assignment (as an optional field with tag 'XF')")
	pa.add_argument(
		"-t", "--type", type=str, dest="featuretype",
		default="exon",
		help="Feature type (3rd column in GTF file) to be used, " +
		"all features of other type are ignored (default, suitable for Ensembl " +
		"GTF files: exon)")
	
	pa.add_argument(
		"-i", "--idattr", type=str, dest="idattr",
		# ~ default="gene_id",
		default="ID",
		help="GTF attribute to be used as feature ID (default, " +
		"suitable for Ensembl GTF files: gene_id)")
			
		# In my case i will define default all for nonunique - multimapped
	pa.add_argument(
		"--multimapped", dest="multimapped", type=str,
		# TODO set default none
		choices=("none", "all"), default="all",
		help="Whether to score reads that are not uniquely aligned " +
		"or ambiguously assigned to features")
	# TODO implement this
	# ~ pa.add_argument(
		# ~ "--supplementary-alignments", dest="supplementary_alignments", type=str,
		# ~ choices=("score", "ignore"), default="ignore",
		# ~ help="Whether to score supplementary alignments (0x800 flag)")
	


	# parser.addArgument("-m", "--mincount")
				# ~ .help("min TOTAL counts for given transcript, default 5").required(false).setDefault(5).dest("mincount");
	args = pa.parse_args()

	# Clean file
	with open(args.samouts, 'w'):
		pass
	count_reads_by_position_in_features(
		args.samfilenames,
		args.featuresfilename,
		args.samouts,
		args.featuretype,
		args.idattr,
		args.multimapped,
		'yes',
		20,
		False)
	# ~ try:
		# ~ # type idattr stranded quiet
		# ~ count_reads_by_position_in_features(
			# ~ args.samfilenames,
			# ~ args.featuresfilename,
			# ~ 'exon',
			# ~ 'ID',
			# ~ 'yes',
			# ~ 20,
			# ~ False)
        # ~ count_reads_in_features(
                # ~ args.samfilenames,
                # ~ args.featuresfilename,
                # ~ args.samtype,
                # ~ args.order,
                # ~ args.max_buffer_size,
                # ~ args.stranded,
                # ~ args.mode,
                # ~ args.nonunique,
                # ~ args.secondary_alignments,
                # ~ args.supplementary_alignments,
                # ~ args.featuretype,
                # ~ args.idattr,
                # ~ args.additional_attr,
                # ~ args.quiet,
                # ~ args.minaqual,
                # ~ args.samouts)
	# ~ except:
		# ~ print("Error")
		# ~ sys.stderr.write("  %s\n" % str(sys.exc_info()[1]))
		# ~ sys.stderr.write("  [Exception type: %s, raised in %s:%d]\n" %
						# ~ (sys.exc_info()[1].__class__.__name__,
						# ~ os.path.basename(traceback.extract_tb(
							# ~ sys.exc_info()[2])[-1][0]),
							# ~ traceback.extract_tb(sys.exc_info()[2])[-1][1]))
	sys.exit(1)

if __name__ == "__main__":
	main()

