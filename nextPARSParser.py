import HTSeq
import itertools
import sys
import argparse
import collections
import csv
import time


# source activate nextPARSParser
# python parser_test.py

# TODO
# Check for perfect match


# Experienced Python developers will recognize that the for loop can be replaced with a single line using a generator comprehension and the reduce function:
# sorted(set.union(*[val for iv, val in gas[ read_iv ].steps()]))



features = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
# TODO remove this line
# features = HTSeq.GenomicArrayOfSets("auto", stranded != "no")

# TODO, put it in an object
features_properties = dict()


def count_reads_by_position_in_features(bam_file,gtf_file,out_file,feature_type,id_attribute,multimapped_mode,stranded,minaqual,quiet):


	gff = HTSeq.GFF_Reader( gtf_file, end_included=True )
	k = 0

	# t = time.process_time()

	# TODO: Mirar si con esto puedo evitar hacer el HACK (creo que si)
	# CIGAR match characters (including alignment match, sequence match, and
	# sequence mismatch
	# Ref https://drive5.com/usearch/manual/cigar.html
	com = ('M', '=', 'X')
	
	try:
		for f in gff:
			# feature type (3rd column in GFF file)
			
			# TODO Esto lo  puse solo para un test
			# Pero contig no va y tengo que ver porqué
			# No va con todos porque mapea en varios lados

			# TODO: Remove this
			# if f.type != 'contig':
			# if f.type == 'exon':
			
			if f.type == feature_type:
				try:
					feature_id = f.attr[id_attribute]
				except KeyError:
					raise ValueError(
							"Feature %s does not contain a '%s' attribute" %
							(f.name, id_attribute))
							
				
				features[ f.iv ] += feature_id
				features_properties[ f.name ] = f.iv.length,f.iv.start,f.iv.end
				# genes = Gene(f.name,f.iv.length,f.iv.start,f.iv.end)
		

				if stranded != "no" and f.iv.strand == ".":
					raise ValueError(
							"Feature %s at %s does not have strand information but you are "
							"running htseq-count in stranded mode. Use '--stranded=no'." %
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



	# elapsed_time = time.process_time() - t
	# print("ELAPSED1 : ",elapsed_time)
	
	# GET COUNTS
	bam = HTSeq.BAM_Reader(bam_file)
	counts = collections.Counter( )
	k = 0
	genes_set = set()


	# t = time.process_time()

	# TODO remove this line
	# and remove itertools
	# ~ for align in itertools.islice( bam, 100000):
	for align in bam:

		
		if k > 0 and k % 100000 == 0 and not quiet:
			sys.stderr.write("%d SAM alignment records processed.\n" % k)
			sys.stderr.flush()
		
		k += 1
		

		
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
			counts[ "__too_long" ] += 1
			# ~ for co in align.cigar:
				# ~ # cut long cigar strings, only the best of reads, maximum one variation
				# ~ if co.type not in com and co.size > 0:
					# ~ print("NO",co,co.type,co.size,co.ref_iv)
				# ~ else:
					# ~ print("SI",co,co.type,co.size,co.ref_iv)
			continue
		
		
		# ~ iv_seq = (co.ref_iv for co in align.cigar if co.type in com
                                  # ~ and co.size > 0)
 
		# ~ fs = set()

		# ~ for iv in iv_seq:
			# ~ if iv.chrom not in features.chrom_vectors:
				# ~ raise UnknownChrom
			# ~ for iv2, fs2 in features[iv].steps():
				# ~ fs = fs.union(fs2)
				
		if align.aQual < minaqual:
			counts[ "__too_low_aQual" ] += 1
			continue

		gene_ids = set()
		
		# interval, gen
		for iv_f, gene in features[ align.iv ].steps():
			# Add element to the set.
			# ~ print ("A",iv_f,gene)
			gene_ids |= gene

		# overlap_mode (same as htseq-count union)

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
	
	# elapsed_time = time.process_time() - t
	# print("ELAPSED2 : ",elapsed_time)
	
	# GET RESULTS
	# New dictionary sort on first element of the key and then on the second element
	# sorted_counts = collections.OrderedDict(sorted(counts.items(),key = lambda i: (i[0][0], i[0][1])))
	# ~ print ("Length: ",len(counts))
	
	# ~ for a in counts:
		# ~ print (a)
	
	# New list sort on gene_name and then on nt position on that gene
	sorted_counts = sorted(counts.items(),key = lambda i: (i[0][0], i[0][1]))

	last_gene_id = ''
	additional_data = []
	list_of_reads = []
	
	# t0 = time.process_time()
	# print (sorted_counts)
	my_iter = iter(sorted_counts)
	
	# TODO rename this variable to number of ??
	p = 0
	for (index, gene_pos_count) in enumerate(sorted_counts[:-1]):
		if index < len(sorted_counts):
			current, next_ = gene_pos_count, sorted_counts[index + 1]
			
			# print("CN",current, next_)			
			# Gene ID and read position
			gene_id_pos = current[0]
			# Read counts for current position
			read_counts = current[1]
			
	
		# t = time.process_time()
		
		# print("gene_id_pos", gene_id_pos)
		
		# _unmapped _no_feature or _ambiguous
		if isinstance(gene_id_pos, str):
			additional_data.append((gene_id_pos,read_counts))
			
		else:
			p += 1
			gene_id = gene_id_pos[0] 
			position = gene_id_pos[1]
			next_gene_id = next_[0][0]
			
			# TODO: Remove this
			# gene_id_length = features_properties.get(gene_id)[0]
			gene_id_start = features_properties.get(gene_id)[1]
			gene_id_end = features_properties.get(gene_id)[2]
		
			# print(gene_id,next_gene_id)

			# elapsed_time = time.process_time() - t
			# print("END0 : ",elapsed_time)
			# t = time.process_time()

		# Check if there is something mapped
		# TODO check this because it should continue to print _unmapped etc
		try: 
			gene_id
		except:
			sys.stderr.write("No genes mapped")
			raise SystemExit

		# Print only the Gene name
		if (last_gene_id != gene_id) :
			i = gene_id_start
			
			# print(gene_id)
			reads = []
			# reads.append(gene_id+"\t")
			reads.append(gene_id)
			last_gene_id = gene_id
		
		
			# elapsed_time = time.process_time() - t
			# print("END1 : ",elapsed_time)
			# t = time.process_time()

		# Still working with the same gene
		if (gene_id == next_gene_id):

			while (i < position):
				# print("pos 0: ", i, "pos ", position, "counts: ",0, end = '\n')
				# TODO llenar con 0 al principio
				# reads.extend(repeat(0, position - i))
				# i = position
				# ~ print(position-i,gene_id,i,position)
				reads.append(0)
				i+=1
			if i == position:
				# print("pos A: ", i, "pos ", position, "gene_id_start",gene_id_start, "counts: ",read_counts, end = '\n')
				# print(current)
				# reads.append(sorted_counts[gene_id_pos])
				reads.append(read_counts)
				# reads.insert(position, sorted_counts[gene_id_pos])
				i+=1
				
		# Last position(s) of the gene
		else:
			if i == position:
				# print("pos B: ", i, "pos ", position, "counts: ", sorted_counts[tuple_e], end = '\n')
				# reads.append(sorted_counts[gene_id_pos])
				reads.append(read_counts)
				# reads.insert(position, sorted_counts[gene_id_pos])

				i+=1
			while (i < gene_id_end):
				# print("pos X: ", i, "pos ", position, "counts: ",0, end = '\n')
				reads.append(0)
				i+=1
			# ~ else:
				# ~ # TODO
				# ~ # complete with 0s
				# ~ # reads.extend(repeat(0, gene_id_end - i))
				# ~ # i = gene_id_end
				# ~ while (i < gene_id_end):
					# ~ reads.append(0)
					# ~ i+=1
			
			# Add reads to the final
			#reads[0] = ''.join(str(reads[0:2]))

			list_of_reads.append(reads)
		
	# Last position of the file
	# TODO put this one
	gene_id_pos = next_[0]
	read_counts = next_[1]
	
	# ~ if isinstance(gene_id_pos, str):
		# ~ print("Entra?")
		# ~ end_table.append((gene_id_pos,read_counts))
		# ~ print(end_table)
		
		# ~ df = pd.DataFrame(data=end_table)
		# ~ # Append mode
		# ~ df.to_csv(FILE_TO_WRITE, mode='a', index=False,header=False)
		
	# ~ else:
		# ~ gene_id = gene_id_pos[0] 
		# ~ position = gene_id_pos[1]
		# ~ next_gene_id = next_[0][0]
		
		# ~ # TODO: Remove this
		# ~ # gene_id_length = features_properties.get(gene_id)[0]
		# ~ gene_id_start = features_properties.get(gene_id)[1]
		# ~ gene_id_end = features_properties.get(gene_id)[2]
	
		# ~ print(gene_id,next_gene_id)
		
		

			# elapsed_time = time.process_time() - t
			# print("END0 : ",elapsed_time)
			# t = time.process_time()

		# elapsed_time = time.process_time() - t
		# print("END2 : ",elapsed_time)

	# elapsed_time = time.process_time() - t0
	# print("END0 : ",elapsed_time)	
	
	# TODO: Remove pandas
	# Append mode
	# df = pd.DataFrame(data=list_of_reads)
	# df.to_csv(FILE_TO_WRITE, mode='a', index=False,header=False)
	
	# Append mode
	write_to_file(out_file,additional_data)
	write_to_file(out_file,list_of_reads)
	# To test I add header 
	# df.to_csv(FILE_TO_WRITE, index=False, mode='a')




# Write to file
def write_to_file (FILE_TO_WRITE,data_to_write):
	with open(FILE_TO_WRITE,"a") as f:
		wr = csv.writer(f,delimiter=';')
		wr.writerows(data_to_write)
		






# TODO add a method to write to a file
# write_to_samout

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

