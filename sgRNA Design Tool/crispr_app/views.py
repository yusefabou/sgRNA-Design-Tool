from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from .models import GeneID
from .forms import NumberForm
from Bio.Seq import Seq
from Bio import Entrez, SeqIO
import sys
sys.path.insert(0, '/Users/Yusef/Documents/Kitematic/crisprapp/hello_django_docker/Project/crispr_app/Rule_Set_2_scoring_v1/analysis')
from rs2_score_calculator import calculateScore
import cPickle as pickle
import dill

Entrez.email='yabourem@ucsd.edu'
# Create your views here.

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
	template = loader.get_template('crispr_app/form.html')
	sequence = ''
	scores = ''
	indices = ''
	if request.method == 'POST':
		#Get sequence
		# gi = request.POST.get('locus')
		# start = request.POST.get('start')
		# stop = request.POST.get('end')
		Entrez.email = 'yabourem@ucsd.edu'
		handle = Entrez.efetch(db='nucleotide',id='EU490707.1',rettype='gb',retmode='text', seq_start=1, seq_stop=200)
		record = SeqIO.read(handle, 'gb')
		sequence = record.seq
		handle.close()
		print type(sequence)
		sequence = str(sequence)
		#Get scores of sequence
		model_file = open('crispr_app/V3_model_nopos.pickle', 'rb')
		model = pickle.load(model_file)
		scores = []
		indices = []
		for i in range(len(sequence) - 30):
			toScore = sequence[i:i+30]
			if len(toScore) == 30 and toScore[25:27] == 'GG':
				print i
				scores.append(calculateScore(toScore, model)) 
				indices.append(i)
		
	context = {
	'seq' : sequence,
	'scores' : scores,
	'indices' : indices
	}
	return HttpResponse(template.render(context, request))


def result(request):
	template = loader.get_template('crispr_app/result.html')
	context = {
	'hi' : 1
	}
	return HttpResponse(template.render(context, request))
