# sgRNA Design Tool 

## Step 1: Install Dependencies
Download and install: </br>
bowtie (not bowtie2) </br> 
django 1.9.5 </br>
python 2.7.9 </br>
biopython (python module)

## Step 2: Change user-specific directories 
In crispr_app/views.py, change all user specific references to bowtie to your directory of bowtie.
In crispr_app/functions.py, change all user specific references to crispr_app/Rule_Set_2_scoring_v1/analysis to your own directory with the analysis executable.

## Step 3: Run django server
Navigate to the parent directory in the terminal and run "python manage.py runserver".

## Step 4: Use the tool
Navigate to the default web page, given in the terminal. The default web page should be "http://127.0.0.1:8000". Once there, enter a gene symbol of interest and a CRISPR analysis category. Alternatively, enter a DNA sequence of appropriate length (>30bp) with at least one PAM sequence (NGG). Then, once the gene or sequence name has been entered, hit "Search" to process the data. Notes: The default gene loaded into the SQLite database is "APOB", and should return results quickly, and without processing any alignments. 