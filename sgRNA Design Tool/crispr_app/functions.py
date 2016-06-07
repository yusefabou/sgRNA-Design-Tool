#Standard Library
import time
import sys
sys.path.insert(0, '/Users/Yusef/Documents/Kitematic/crisprapp/hello_django_docker/sgRNA Design Tool/crispr_app/Rule_Set_2_scoring_v1/analysis')

#Third-Party
from Bio.Seq import Seq
from Bio import Entrez, SeqIO, SeqFeature
from BCBio import GFF
from rs2_score_calculator import calculateScore
import cPickle as pickle
import dill

def chromosome_to_accession():
	return {'1':'NC_000001.11' ,'2':'NC_000002.12', '3':'NC_000003.12', '4':'NC_000004.12', 
	'5':'NC_000005.10', '6':'NC_000006.12', '7':'NC_000007.14', '8':'NC_000008.11', '9':'NC_000009.12',
	'10':'NC_000010.11', '11':'NC_000011.10', '12':'NC_000012.12', '13':'NC_000013.11', '14':'NC_000014.9',
	'15':'NC_000015.10', '16':'NC_000016.10', '17':'NC_000017.11', '18':'NC_000018.10', '19':'NC_000019.10',
	'20':'NC_000020.11', '21':' NC_000021.9', '22':'NC_000022.11', 'X':'NC_000023.11', 'Y':'NC_000024.10'}

def get_range(crispr_type, strand, transcription_start):
	start = ''
	stop = ''
	if crispr_type == 'activation' and strand == 1:
		start = transcription_start
		stop = transcription_start - 250
	elif crispr_type == 'activation' and strand == -1:
		start = transcription_start
		stop = transcription_start + 250
	elif crispr_type == 'repression' and strand == 1:
		start = transcription_start
		stop = transcription_start + 150
	elif crispr_type == 'repression' and strand == -1:
		start = transcription_start
		stop = transcription_start - 150
	return (start, stop)

def gene_to_TSS(gene_name):
	#Initialize variables
	transcription_start = ''
	strand = ''
	chromosome = ''
	found_gene = False

	#Open annotation file
	annotation_file = 'crispr_app/Homo_sapiens.GRCh38.84.gtf'
	limit_info = dict(
	         gff_type = ["transcript"])
	annotation_handle = open(annotation_file)

	#Parse through annotated data, searching for matching gene names
	for rec in GFF.parse(annotation_handle, limit_info=limit_info, target_lines=1):
		feature = rec.features[0]
		qualifiers = feature.qualifiers

		#Once matching gene is found, determine the transcription start site and chromosome
		if str(qualifiers['gene_name']).strip('[').strip(']').strip('\'') == gene_name:
			found_gene = True
			chromosome = rec.id
			strand = feature.strand
			if strand == 1:
				if not transcription_start:
					transcription_start = float('inf')
				transcription_start = min(int(feature.location.start), int(transcription_start))
			elif strand == -1:
				if not transcription_start:
					transcription_start = -1
				transcription_start = max(int(feature.location.end), int(transcription_start))
		elif found_gene == True:
			break
	annotation_handle.close()
	return (transcription_start, strand, chromosome)

def gene_to_early_exons(gene_name, num_exons):
	#Initialize variables
	exons = {}
	exonCount = 0
	maxExons = 0

	#Open annotation file
	annotation_file = 'crispr_app/Homo_sapiens.GRCh38.84.gtf'
	limit_info = dict(
	         gff_type = ["exon"])
	annotation_handle = open(annotation_file)

	#Parse through annotated data, searching for matching gene names & exons 
	strand = ''
	for rec in GFF.parse(annotation_handle, limit_info=limit_info, target_lines=1):
		feature = rec.features[0]
	 	qualifiers = feature.qualifiers

	 	#Once matching gene is found, determine the exon regions and chromosome
	 	if str(qualifiers['gene_name']).strip('[').strip(']').strip('\'') == gene_name:
	 		chromosome = rec.id
	 		strand = feature.strand

	 		#Get only first version of gene in annotated data
	 		exonNum = str(qualifiers['exon_number']).strip('[').strip(']').strip('\'')
	 		maxExons = max(maxExons, int(exonNum))
	 		exonCount +=1 
			if exonCount > maxExons:
				break
			if exonCount > num_exons:
				break
			exons[exonNum] = [int(feature.location.start), int(feature.location.end), strand]
	annotation_handle.close()
	return exons, chromosome

def sequence_from_NCBI(chromosome, start, stop):
	#Give Entrez email in order to access API
	Entrez.email = 'yabourem@ucsd.edu'

	#Request information and pause to avoid overusage of API
	chromosome_accessions = chromosome_to_accession()
	handle = Entrez.efetch(db='nucleotide',id=chromosome_accessions[chromosome],rettype='gb',
		retmode='text', seq_start=start, seq_stop=stop)
	time.sleep(1)

	#Determine sequence
	record = SeqIO.read(handle, 'gb')
	sequence = record.seq
	handle.close()
	sequence = str(sequence)
	if stop < start:
		mySeq = Seq(sequence)
		sequence = str(mySeq.reverse_complement())
	return sequence

def find_and_score(sequence):
	#Initialize arrays
	scores = []
	indices = []
	complements = []
	PAMs = []
	
	#Load scoring file
	model_file = open('crispr_app/V3_model_nopos.pickle', 'rb')
	model = pickle.load(model_file)

	#Score in 5-->3 direction
	for i in range(len(sequence) - 30):
		toScore = sequence[i:i+30]
		if len(toScore) == 30 and toScore[25:27] == 'GG':
			complements.append(toScore[4:24])
			PAMs.append(toScore[24:27])
			scores.append(calculateScore(toScore, model))
			indices.append(i+21)

	#Score in 3-->5 (Reverse complement) direction
	mySeq = Seq(sequence)
	reverseComp = str(mySeq.reverse_complement())
	for i in range(len(reverseComp) - 30):
		toScore = reverseComp[i:i+30]
		if len(toScore) == 30 and toScore[25:27] == 'GG':
			complements.append(toScore[4:24])
			PAMs.append(toScore[24:27])
			scores.append(calculateScore(toScore, model)) 
			indices.append(len(sequence)-(i+21))
	return scores, indices, complements, PAMs
