# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 12:50:23 2021

@author: makow
"""

import random
random.seed(16)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sc
import seaborn as sns
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score, mean_absolute_error, mean_absolute_percentage_error
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_validate as cv
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import itertools
from sklearn.linear_model import Ridge

#colormap used in figures
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
cmap = plt.cm.get_cmap('plasma')
colormap9= np.array([cmap(0.25),cmap(0.77)])
cmap9 = LinearSegmentedColormap.from_list("mycmap", colormap9)

variant_colors = ['red', 'navy', 'blueviolet', 'green', 'aqua', 'magenta', 'yellow']
cmap_var = ListedColormap(variant_colors)

def lighten_color(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


wt_seq = 'NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST'
variant_seqs = {'WT':'NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST', 'B.1.1.7':'NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKST', 'B.1.351': 'NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVKGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKST', 'B.1.429': 'NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYRYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST', 'P.1': 'NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGTIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVKGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKST', 'B.1.617.2': 'NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYRYRLFRKSNLKPFERDISTEIYQAGSKPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST', 'B.1.529': 'NITNLCPFDEVFNATRFASVYAWNRKRISNCVADYSVLYNLAPFFTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKVSGNYNYLYRLFRKSNLKPFERDISTEIYQAGNKPCNGVAGFNCYFPLRSYSFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST'}
variant_seqs = pd.DataFrame.from_dict(variant_seqs, orient = 'index')

wt_dna = list('AATATTACAAACTTGTGCCCTTTTGGTGAAGTTTTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGATTATTCTGTCCTATATAATTCCGCATCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAAAAGTCTACT')
delta_dna = list('AATATTACAAACTTGTGCCCTTTTGGTGAAGTTTTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGATTATTCTGTCCTATATAATTCCGCATCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCGGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAAAAGTCTACT')

# indices of wt AA that when mutated to S or T introduce glysolylation sites
glycan_end_st = [2,5,14,25,31,41,59,65,88,93,108,110,111,119,121,131,148,152,158,172]
# indices of wt AA that when mutaed away from N, S, or T introduce glycosylation sites
glycan_begin = [0,2,12,14]
# indices of wt AA that when mutaed to N introduce glycosylation
glycan_end_n = [16,26,33,38,40,42,43,50,52,60,66,82,97,105,110,126,136,137,144,145,161,167,181,190,197,198]

dna_table = {
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'C': ['TGC', 'TGT'],
    'D': ['GAC', 'GAT'],
    'E': ['GAA', 'GAG'],
    'F': ['TTC', 'TTT'],
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'H': ['CAC', 'CAT'],
    'I': ['ATA', 'ATC', 'ATT'],
    'K': ['AAA', 'AAG'],
    'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'],
    'Q': ['CAA', 'CAG'],
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
    'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'W': ['TGG'],
    'Y': ['TAC', 'TAT'],
    '_': ['TAA', 'TAG', 'TGA'],}

def get_key(val):
    for key, value in dna_table.items():
         if val == value:
             return key
 
    return "key doesn't exist"

#integer encoding amino acids for one hot encoding RBD sequences
alph_letters = np.array(sorted('ACDEFGHIKLMNPQRSTVWY'))
le = LabelEncoder()
integer_encoded_letters = le.fit_transform(alph_letters)
integer_encoded_letters = integer_encoded_letters.reshape(len(integer_encoded_letters), 1)
one = OneHotEncoder(sparse = False)
ohe_letters = one.fit_transform(integer_encoded_letters)

#one hot encoding wt sequence
wt_enc = le.transform(list(wt_seq))
ohe_let = pd.DataFrame(one.transform(wt_enc.reshape(len(list(wt_seq)), 1)))
wt_ohe = ohe_let.values.flatten()
wt_ohe = pd.DataFrame(wt_ohe).T


def ace_binding_prepro(ace_binding_github):
    ace_binding_group = ace_binding_github.groupby(['aa_substitutions']).mean()
    ace_binding_group.dropna(axis = 0, inplace = True)
    ace_binding = []
    for index, row in ace_binding_group.iterrows():
        muts = index.split()
        seq = list(wt_seq)
        for i in muts:
            mut = list(i)
            remove = i[0]
            add = i[-1]
            mut.pop(0)
            mut.pop(-1)
            res = int(''.join(mut)) - 1
            seq[res] = add
        ace_binding.append([''.join(seq), row[2], row[8], row[1], muts])
    ace_binding = pd.DataFrame(ace_binding)
    ace_binding = ace_binding[ace_binding[1] > 6]
    ace_binding = ace_binding[ace_binding[1] < 13]
    return ace_binding


def pAb_escape_prepro(data):
    pAb_escape_group = data.groupby(['aa_substitutions']).mean()
    pAb_escape_group.dropna(axis = 0, inplace = True)
    pAb_escape = []
    for index, row in pAb_escape_group.iterrows():
        muts = index.split()
        seq = list(wt_seq)
        for i in muts:
            mut = list(i)
            remove = i[0]
            add = i[-1]
            mut.pop(0)
            mut.pop(-1)
            res = int(''.join(mut)) - 1
            seq[res] = add
        pAb_escape.append([''.join(seq), row[0], row[5], row[3], muts])
    pAb_escape = pd.DataFrame(pAb_escape)
    return pAb_escape


def ohe_encode(sequence_df):
    enc= []
    ohe = []
    for i in sequence_df.iloc[:,0]:
        chars = le.transform(list(i))
        enc.append(chars)
    enc = pd.DataFrame(enc)
    for index, row in enc.iterrows():
        enc_row = np.array(row)
        let = enc_row.reshape(201,1)
        ohe_let = pd.DataFrame(one.transform(let))
        ohe.append(ohe_let.values.flatten())
    ohe = np.stack(ohe)
    
    return ohe







