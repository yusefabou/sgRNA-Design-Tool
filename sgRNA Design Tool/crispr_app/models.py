from django.db import models
from django.forms import ModelForm

class GeneID(models.Model):
	geneID_text = models.CharField(max_length=200)
	def __str__(self):
		return self.geneID_text

sgRNA_PURPOSES = (('Knockout', 'KNOCKOUT'), 
	('Regulation', 'REGULATION'),)

class Purpose(models.Model):
	geneID = models.ForeignKey(GeneID)
	purpose = models.CharField(max_length=15, choices=sgRNA_PURPOSES)
	def __str__(self):
		return self.purpose

# class sgRNAForm(ModelForm):
# 	class Meta:
# 		model = GeneID
# 		fields = ['geneID_text']
# 	