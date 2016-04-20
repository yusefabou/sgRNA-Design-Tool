#Calculates the Rule set 2 score for the given 30-mer
#Input: 1. 30mer sgRNA+context sequence, NNNN[sgRNA sequence]NGGNNN
#       2. Amino acid cut position, for full model prediction only
#       3. Percent peptide, for full model prediction only
#Output: Rule set 2 score

import pandas as pd
import csv, sys
import pickle
import model_comparison

def calculateScore(seq, model):
    if seq[25:27] == 'GG':
        score = model_comparison.predict(seq, -1, -1, model)
        #print 'Rule set 2 score: %.4f'% (score)
        return str(score)
    else:
        print >> sys.stderr, 'Calculates on-target scores for sgRNAs with NGG PAM only.'

if __name__ == '__main__':
    holder = 'holder'