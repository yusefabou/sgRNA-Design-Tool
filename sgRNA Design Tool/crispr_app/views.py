#Standard Library
import sys
sys.path.insert(0, 'Users/Yusef/Documents/bowtie/')
import subprocess
import os 
from operator import itemgetter

#Third-Party
from Bio.Seq import Seq
from Bio import Entrez, SeqIO, SeqFeature
from functions import *

#Django
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from .models import RNA

def index(request):
	#Load HTML template
	template = loader.get_template('crispr_app/index.html')

	#Store accession numbers of every chromosome
	chromosome_accessions = chromosome_to_accession()
	
	#Initialize context variables
	sequence = ''
	sgRNAInfo = []
	species = ''
	gene_name = ''
	crispr_type = ''
	chromosome = ''

	#Populate variables when user submits form
	if request.method == 'POST':

		#Refresh page if fields empty
		if not request.POST.get('gene_name') and not request.POST.get('singleSeq'):
			return HttpResponse(template.render({}, request))

		#Get gene name if submitted
		if request.POST.get('gene_name'):
			gene_name = request.POST.get('gene_name')
			gene_name = gene_name.strip().upper()

		#Get sequence if submitted
		else:
			sequence = request.POST.get('singleSeq')
			sequence = sequence.strip().upper()
			#Refresh page if sequence not an appropriate length 
			if len(sequence) < 30:
				return HttpResponse(template.render({}, request))

		#Get species
		#NOTE: Only human hg38 supported; must download other genomes and annotation files 
		#to support other species
		species = request.POST.get('species')

		#Fetch sequence data if a gene name is used
		if gene_name:

			#Find transcription start site if CRISPR will be used for activation/repression
			#NOTE: Location data reported relative to reference genome (ex: hg38)
			crispr_type = request.POST.get('crispr_type')
			if crispr_type == 'activation' or crispr_type == 'repression':

				#Find genomic coordinates of gene
				transcription_start, strand, chromosome = gene_to_TSS(gene_name)

				#Set start and stop points
				#TODO: Coordinates working on reverse (-1) strand. Check on forward (+) strand.
				start, stop = get_range(crispr_type, strand, transcription_start)

				#Fetch sequence data from NCBI Nucleotide database
				if start and stop:
					sequence = sequence_from_NCBI(chromosome, start, stop)

				#Calculate scores of sequence, and collect other sgRNA information
				scores, indices, complements, PAMs = find_and_score(sequence)

				#Calculate alignments of complement+PAM sequence via bowtie
				mismatch_dict = {}
				indexname = '/Users/Yusef/Documents/bowtie/indexes/GCA_000001405.15_GRCh38_no_alt_analysis_set'
				for i in range(len(complements)):
					molecule = complements[i] + PAMs[i]
					if RNA.objects.filter(sequence=molecule, species=species).count() == 1:
						continue
					mismatch_dict[molecule] = [-1, 0, 0, 0]
					alignments = subprocess.Popen(['/Users/Yusef/Documents/bowtie/bowtie', '-c', '-a', '-v 3', indexname, molecule], stdout=subprocess.PIPE)
					for line in alignments.stdout:
						mismatches = line.count(':')
						mismatch_dict[molecule][mismatches] += 1
					#Add the complement+PAM and alignment information to database
					rna = RNA.objects.create(sequence=molecule, species=species, 
						zero_mismatch_count=mismatch_dict[molecule][0], 
						one_mismatch_count=mismatch_dict[molecule][1], 
						two_mismatch_count=mismatch_dict[molecule][2], 
						three_mismatch_count=mismatch_dict[molecule][3])

				#Visualize cut points in the sequence
				for cut in sorted(indices, reverse=True):
					sequence = sequence[:cut] + '['  + str(cut) + ']' + sequence[cut:]

				#Put all sgRNA info into one array
				#Note to self: careful with indices, as the user input starts at 1
				for i in range(len(scores)):
					rna = RNA.objects.get(sequence=complements[i]+PAMs[i], species=species)
					temp = [complements[i], PAMs[i], scores[i], indices[i], rna.zero_mismatch_count, 
					rna.one_mismatch_count, rna.two_mismatch_count, rna.three_mismatch_count]
					sgRNAInfo.append(temp)
			
			if crispr_type == 'knock_out':
				#Find genomic cooardinates of exons
				exons, chromosome = gene_to_early_exons(gene_name, 3)

				#Fetch sequences from NCBI Nucleotide database
				Entrez.email = 'yabourem@ucsd.edu'
				for exon, info in exons.iteritems():
					start = exons[exon][0]
					stop = exons[exon][1]
					sequence = sequence_from_NCBI(chromosome, start, stop)
					exons[exon].append(sequence)

					#Calculate scores of sequence, and collect other sgRNA information
					scores, indices, complements, PAMs = find_and_score(sequence)

					#Calculate alignments of complement+PAM sequence via bowtie
					mismatch_dict = {}
					indexname = '/Users/Yusef/Documents/bowtie/indexes/GCA_000001405.15_GRCh38_no_alt_analysis_set'
					for i in range(len(complements)):
						molecule = complements[i] + PAMs[i]
						if RNA.objects.filter(sequence=molecule, species=species).count() == 1:
							continue
						mismatch_dict[molecule] = [-1, 0, 0, 0]
						alignments = subprocess.Popen(['/Users/Yusef/Documents/bowtie/bowtie', '-c', '-a', '-v 3', indexname, molecule], stdout=subprocess.PIPE)
						for line in alignments.stdout:
							mismatches = line.count(':')
							mismatch_dict[molecule][mismatches] += 1
						#Add the complement+PAM and alignment information to database
						rna = RNA.objects.create(sequence=molecule, species=species, 
							zero_mismatch_count=mismatch_dict[molecule][0], 
							one_mismatch_count=mismatch_dict[molecule][1], 
							two_mismatch_count=mismatch_dict[molecule][2], 
							three_mismatch_count=mismatch_dict[molecule][3])

					#Put all sgRNA info into one array
					for i in range(len(scores)):
						rna = RNA.objects.get(sequence=complements[i]+PAMs[i], species=species)
						temp = [complements[i], PAMs[i], scores[i], indices[i], rna.zero_mismatch_count, 
						rna.one_mismatch_count, rna.two_mismatch_count, rna.three_mismatch_count]
						sgRNAInfo.append(temp)
					exons[exon].append(sorted(sgRNAInfo, key=itemgetter(2), reverse=True))

				context = {
				'gene_name' : gene_name,
				'exons' : exons,
				'species' : species
				}
				return HttpResponse(template.render(context, request))
		
		#If an individual sequence is given, skip fetching for sequence through NCBI
		elif sequence:
			gene_name = 'N/A'

			#Calculate scores of sequence, and collect other sgRNA information
			scores, indices, complements, PAMs = find_and_score(sequence)

			#Calculate alignments of complement+PAM sequence via bowtie
			mismatch_dict = {}
			indexname = '/Users/Yusef/Documents/bowtie/indexes/GCA_000001405.15_GRCh38_no_alt_analysis_set'
			for i in range(len(complements)):
				molecule = complements[i] + PAMs[i]
				if RNA.objects.filter(sequence=molecule, species=species).count() == 1:
					continue
				mismatch_dict[molecule] = [-1, 0, 0, 0]
				alignments = subprocess.Popen(['/Users/Yusef/Documents/bowtie/bowtie', '-c', '-a', '-v 3', indexname, molecule], stdout=subprocess.PIPE)
				for line in alignments.stdout:
					mismatches = line.count(':')
					mismatch_dict[molecule][mismatches] += 1
				#Add the complement+PAM and alignment information to database
				rna = RNA.objects.create(sequence=molecule, species=species, 
					zero_mismatch_count=mismatch_dict[molecule][0], 
					one_mismatch_count=mismatch_dict[molecule][1], 
					two_mismatch_count=mismatch_dict[molecule][2], 
					three_mismatch_count=mismatch_dict[molecule][3])

			#Visualize cut points in the sequence
			for cut in sorted(indices, reverse=True):
				sequence = sequence[:cut] + '['  + str(cut) + ']' + sequence[cut:]

			#Put all sgRNA info into one array
			#Note to self: careful with indices, as the user input starts at 1
			for i in range(len(scores)):
				rna = RNA.objects.get(sequence=complements[i]+PAMs[i], species=species)
				temp = [complements[i], PAMs[i], scores[i], indices[i], rna.zero_mismatch_count, 
				rna.one_mismatch_count, rna.two_mismatch_count, rna.three_mismatch_count]
				sgRNAInfo.append(temp)

	context = {
		'gene_name' : gene_name,
		'seq' : sequence,
		'sgRNAInfo' : sorted(sgRNAInfo, key=itemgetter(2), reverse=True),
		'species' : species
	}
	
	return HttpResponse(template.render(context,request))

def form(request):
	template = loader.get_template('crispr_app/form.html')
	context = {
	}
	return HttpResponse(template.render(context, request))


def guidelines(request):
	template = loader.get_template('crispr_app/guidelines.html')
	context = {
	}
	return HttpResponse(template.render(context, request))
