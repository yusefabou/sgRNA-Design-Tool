from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from .models import GeneID
from .forms import NumberForm
from Bio.Seq import Seq
from Bio import Entrez, SeqIO
import sys
sys.path.insert(0, '/Users/Yusef/Documents/Kitematic/crisprapp/hello_django_docker/sgRNA Design Tool/crispr_app/Rule_Set_2_scoring_v1/analysis')
from rs2_score_calculator import calculateScore
import cPickle as pickle
import dill
import time

Entrez.email='yabourem@ucsd.edu'

def index(request):
	template = loader.get_template('crispr_app/index.html')
	number = ''
	if request.method == 'POST':
		# locus = request.POST.get('locus')
		# start = request.POST.get('start')
		# end = request.POST.get('end')
		# print start, end, locus
		print request.POST['number']

	context = {
		'number' : number,	
	}
	output = number
	
	return HttpResponse(template.render(context,request))

def form(request):
	#Load Template
	template = loader.get_template('crispr_app/form.html')

	#Initialize context variables
	sequence = ''
	sgRNAInfo = []

	if request.method == 'POST':
		#Refresh page if fields empty
		# if not gi and not start and not stop:
		# 	return HttpResponse(template.render({}, request))

		#Get user input 
		# gi = request.POST.get('locus')
		# start = request.POST.get('start')
		# stop = request.POST.get('end')

		#Fetch sequence data
		Entrez.email = 'yabourem@ucsd.edu'
		handle = Entrez.efetch(db='nucleotide',id='KP641121',rettype='gb',retmode='text', seq_start=1, seq_stop=200)
		#handle = Entrez.efetch(db='nucleotide',id=gi,rettype='gb',retmode='text', seq_start=start, seq_stop=stop)
		time.sleep(1)
		record = SeqIO.read(handle, 'gb')
		sequence = record.seq
		handle.close()
		sequence = str(sequence)

		#Calculate scores of sequence, and collect other sgRNA information
		model_file = open('crispr_app/V3_model_nopos.pickle', 'rb')
		model = pickle.load(model_file)
		scores = []
		indices = []
		complements = []
		PAMs = []
		#5-->3
		for i in range(len(sequence) - 30):
			toScore = sequence[i:i+30]
			if len(toScore) == 30 and toScore[25:27] == 'GG':
				complements.append(toScore[4:24])
				PAMs.append(toScore[24:27])
				scores.append(calculateScore(toScore, model)) 
				indices.append(i+21)

		#3-->5 (Reverse complement)
		mySeq = Seq(sequence)
		reverseComp = str(mySeq.reverse_complement())
		for i in range(len(reverseComp) - 30):
			toScore = reverseComp[i:i+30]
			if len(toScore) == 30 and toScore[25:27] == 'GG':
				complements.append(toScore[4:24])
				PAMs.append(toScore[24:27])
				scores.append(calculateScore(toScore, model)) 
				indices.append(200-(i+21))
		
		#Visualize cut points in the sequence
		for cut in sorted(indices, reverse=True):
			sequence = sequence[:cut] + '['  + str(cut) + ']' + sequence[cut:]

		#Put all sgRNA info into one array
		#Note to self: careful with indices, as the user input starts at 1
		for i in range(len(scores)):
			temp = [scores[i], indices[i], complements[i], PAMs[i]]
			sgRNAInfo.append(temp)
	
	context = {
	'seq' : sequence,
	'sgRNAInfo' : sgRNAInfo

	}
	return HttpResponse(template.render(context, request))


def result(request):
	template = loader.get_template('crispr_app/result.html')
	context = {
	}
	return HttpResponse(template.render(context, request))
