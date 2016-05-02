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

class RNA(models.Model):
	sequence = models.CharField(max_length=30)
	species = models.CharField(max_length=20)
	zero_mismatch_count = models.IntegerField(blank=True, null=True)
	one_mismatch_count = models.IntegerField(blank=True, null=True)
	two_mismatch_count = models.IntegerField(blank=True, null=True)
	three_mismatch_count = models.IntegerField(blank=True, null=True)
	four_mismatch_count = models.IntegerField(blank=True, null=True)
	five_mismatch_count = models.IntegerField(blank=True, null=True)
	six_mismatch_count = models.IntegerField(blank=True, null=True)
	off_target_score = models.DecimalField(blank=True, null=True, max_digits=5, decimal_places=2)

	def __str__(self):
		return self.sequence + ' ' + self.species


	
	