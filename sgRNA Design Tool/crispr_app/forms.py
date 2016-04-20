from django import forms

class NumberForm(forms.Form):
	number = forms.CharField(label='Number', max_length=10)