# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 15:16:27 2021

@author: makow
"""

from utils import *
import itertools

#ace data import and training
ace_github = pd.read_csv("C:\\Users\\makow\\Documents\\GitHub\\SARS-CoV-2-RBD_DMS\\results\\binding_Kds\\binding_Kds.csv", header = 0, index_col = 0)
ace_binding = ace_binding_prepro(ace_github)
ace_binding = ace_binding[ace_binding[3] > 53]
ace_binding.reset_index(drop = True, inplace = True)
ace_binding_ohe = ohe_encode(ace_binding)

#ace model training
ace_ridge = Ridge(alpha = 1.90)
ace_ridge.fit(ace_binding_ohe, ace_binding.iloc[:,1])
ace_binding_predict = pd.DataFrame(ace_ridge.predict(ace_binding_ohe))

#pAb data import and processing
pAb_github = pd.read_csv("C:\\Users\\makow\\Documents\\GitHub\\SARS-CoV-2-RBD_MAP_HAARVI_sera\\results\\escape_scores\\scores.csv", header = 0, index_col = 0)
pAb_escape = pAb_escape_prepro(pAb_github)
pAb_escape = pAb_escape[pAb_escape[3] > 12]
pAb_escape.reset_index(drop = True, inplace = True)
pAb_escape_ohe = ohe_encode(pAb_escape)

#pAb model trianing
pAb_ridge = Ridge(alpha = 12.07)
pAb_ridge.fit(pAb_escape_ohe, pAb_escape.iloc[:,1])
pAb_escape_predict = pd.DataFrame(pAb_ridge.predict(pAb_escape_ohe))

#variant one hot encoding and model prediction
variant_ohe = ohe_encode(variant_seqs)
ace_variant = ace_ridge.predict(variant_ohe)
pAb_variant = pAb_ridge.predict(variant_ohe)


#%%
mut_count_ace = pd.read_csv("C:\\Users\\makow\\Documents\\Research\\Data Analysis\\3.12.21_covid_comp\\mut_count_ace.csv", header = 0, index_col = 0)
mut_count_ace.columns = np.arange(0,201)
mut_count_pAb = pd.read_csv("C:\\Users\\makow\\Documents\\Research\\Data Analysis\\3.12.21_covid_comp\\mut_count_pAb.csv", header = 0, index_col = 0)
mut_count_pAb.columns = np.arange(0,201)

sites = np.arange(331,532,1)
sampling = pd.DataFrame(index = sites, columns = sites)
for i in sites:
    for ii in sites:
        if i !=ii:
            sum_res1 = sum(mut_count_ace[i-331] > 6)
            sum_res2 = sum(mut_count_ace[ii-331] > 6)
            sampling.loc[i,ii] = (sum_res1*sum_res2)


sites = np.arange(331,532,1)
sampling_voc = pd.DataFrame(index = [0,1,2,3,4], columns = sites)
for ii in sites:
        sum_res2 = sum(mut_count_ace[ii-331] > 6)
        sum_res1 = sum(mut_count_pAb[ii-331] > 6)
        res1 = min([sum_res1, sum_res2])
        sampling_voc.loc[0,ii] = res1
        sampling_voc.loc[1,ii] = res1
        sampling_voc.loc[2,ii] = res1
        sampling_voc.loc[3,ii] = res1
        sampling_voc.loc[4,ii] = res1
sampling_voc = sampling_voc.astype(int)

site_numbers = pd.DataFrame(np.arange(0,201,1), columns = ['index'])
                    
var_scan_seqs = []
for k in variant_seqs.iterrows():
    for i in np.arange(0,201,1):
        for j in alph_letters:
            var_mut = list(k[1][0]).copy()
            res = var_mut[i]
            if res != j:
                if mut_count_ace.loc[j,i] > 6:
                    if mut_count_pAb.loc[j,i] > 6:
                        var_mut[i] = j
                        var_scan_seqs.append([''.join(var_mut), k[0], res, i, j, k[1][0]])
var_scan_seqs = pd.DataFrame(var_scan_seqs)


single_nuc = []
for index, row in var_scan_seqs.iterrows():
    orig_dna = str(''.join(wt_dna[(row[3]*3):((row[3]*3)+3)]))
    sub_dna = dna_table[row[4]]
    hams = []
    for i in sub_dna:
        hams.append(sc.spatial.distance.hamming(list(i), list(orig_dna))*3)
    single_nuc.append(min(hams))
single_nuc = pd.DataFrame(single_nuc)
single_nuc.set_index(var_scan_seqs.index, inplace = True)
var_scan_seqs = pd.concat([var_scan_seqs, single_nuc], axis = 1, ignore_index = False)
var_scan_seqs.columns = np.arange(0,7)

var_scan_ohe = ohe_encode(var_scan_seqs)


#%%
### need to run multi-objective_opt script first if want double mutant background ###
ace_variant_scan = pd.DataFrame(ace_ridge.predict(var_scan_ohe))
pAb_variant_scan = pd.DataFrame(pAb_ridge.predict(var_scan_ohe))

variant_scan_predictions = pd.concat([ace_variant_scan, pAb_variant_scan], axis = 1)
variant_scan_predictions.columns = ['ACE', 'pAb']

for i in np.arange(0,6):
    seq_ind = np.arange(0,20223,2889)
    variant_scan_predictions_sub = variant_scan_predictions.iloc[seq_ind[i]:seq_ind[i+1]]
    variant_scan_predictions_both = variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i])),:]
    plt.figure()
    plt.scatter(co_op_samp.iloc[0:5000,0], co_op_samp.iloc[0:5000,1]*100, s = 75, c = 'white', edgecolor = 'k', linewidth = 0.25)
    plt.scatter(variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i])  & (var_scan_seqs.iloc[:,6] == 1)),'ACE'], variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i]) & (var_scan_seqs.iloc[:,6] == 1)),'pAb']*100, s = 125, c = lighten_color(variant_colors[i], 0.15), edgecolor = 'k', linewidth = 0.25)
    plt.scatter(ace_variant[i], pAb_variant[i]*100, c = variant_colors[i], s = 230, edgecolor = 'k', linewidth = 0.25)
    plt.scatter(ace_variant[0], pAb_variant[0]*100, c = variant_colors[0], s = 230, edgecolor = 'k', linewidth = 0.25)
    plt.scatter(max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'ACE']), variant_scan_predictions_both.loc[(variant_scan_predictions_both['ACE'] == max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'ACE'])),'pAb']*100, s = 240, c = lighten_color(variant_colors[i], 0.85), edgecolor = 'k', linewidth = 0.25, marker = '^')
    plt.scatter(variant_scan_predictions_both.loc[(variant_scan_predictions_both['pAb'] == max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'pAb'])),'ACE'], max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'pAb'])*100, s = 200, c = lighten_color(variant_colors[i], 0.85), edgecolor = 'k', linewidth = 0.25, marker = 's')
    plt.xlim(10.55, 12.25)
    plt.ylim(3.25, 14.25)
    plt.xticks([10.8, 11.2, 11.6, 12.0], fontsize = 26)
    plt.yticks(fontsize = 26)

#%%
### WT inset
for i in np.arange(0,1):
    seq_ind = np.arange(0,20223,2889)
    variant_scan_predictions_sub = variant_scan_predictions.iloc[seq_ind[i]:seq_ind[i+1]]
    variant_scan_predictions_both = variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i])),:]
    plt.figure()
    plt.scatter(co_op_samp.iloc[0:5000,0], co_op_samp.iloc[0:5000,1]*100, s = 50, c = 'white', edgecolor = 'k', linewidth = 0.25)
    plt.scatter(variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i])  & (var_scan_seqs.iloc[:,6] == 1)),'ACE'], variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i]) & (var_scan_seqs.iloc[:,6] == 1)),'pAb']*100, s = 125, c = lighten_color(variant_colors[i], 0.15), edgecolor = 'k', linewidth = 0.25)
    plt.scatter(ace_variant[i], pAb_variant[i]*100, c = variant_colors[i], s = 230, edgecolor = 'k', linewidth = 0.25)
    plt.scatter(max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'ACE']), variant_scan_predictions_both.loc[(variant_scan_predictions_both['ACE'] == max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'ACE'])),'pAb']*100, s = 240, c = lighten_color(variant_colors[i], 0.85), edgecolor = 'k', linewidth = 0.25, marker = '^')
    plt.scatter(variant_scan_predictions_both.loc[(variant_scan_predictions_both['pAb'] == max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'pAb'])),'ACE'], max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'pAb'])*100, s = 200, c = lighten_color(variant_colors[i], 0.85), edgecolor = 'k', linewidth = 0.25, marker = 's')
    plt.xlim(10.95, 11.45)
    plt.ylim(3.5, 9.25)
    plt.xticks([11.0, 11.1, 11.2, 11.3, 11.4], fontsize = 26)
    plt.yticks(fontsize = 26)

#%%
### b.1.351 inset
for i in np.arange(2,3):
    seq_ind = np.arange(0,20223,2889)
    variant_scan_predictions_sub = variant_scan_predictions.iloc[seq_ind[i]:seq_ind[i+1]]
    variant_scan_predictions_both = variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i])),:]
    plt.figure()
    plt.scatter(co_op_samp.iloc[0:5000,0], co_op_samp.iloc[0:5000,1]*100, s = 50, c = 'white', edgecolor = 'k', linewidth = 0.25)
    plt.scatter(variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i])  & (var_scan_seqs.iloc[:,6] == 1)),'ACE'], variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i]) & (var_scan_seqs.iloc[:,6] == 1)),'pAb']*100, s = 125, c = lighten_color(variant_colors[i], 0.3), edgecolor = 'k', linewidth = 0.25)
    plt.scatter(ace_variant[i], pAb_variant[i]*100, c = variant_colors[i], s = 230, edgecolor = 'k', linewidth = 0.25)
    plt.scatter(max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'ACE']), variant_scan_predictions_both.loc[(variant_scan_predictions_both['ACE'] == max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'ACE'])),'pAb']*100, s = 240, c = lighten_color(variant_colors[i], 0.90), edgecolor = 'k', linewidth = 0.25, marker = '^')
    plt.scatter(variant_scan_predictions_both.loc[(variant_scan_predictions_both['pAb'] == max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'pAb'])),'ACE'], max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'pAb'])*100, s = 200, c = lighten_color(variant_colors[i], 0.90), edgecolor = 'k', linewidth = 0.25, marker = 's')
    plt.xlim(10.6, 11.1)
    plt.ylim(8.5, 14)
    plt.xticks([10.6,10.7, 10.8, 10.9, 11.0, 11.1], fontsize = 26)
    plt.yticks(fontsize = 26)
    
#%%
### delta inset
for i in np.arange(5,6):
    seq_ind = np.arange(0,20223,2889)
    variant_scan_predictions_sub = variant_scan_predictions.iloc[seq_ind[i]:seq_ind[i+1]]
    variant_scan_predictions_both = variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i])),:]
    plt.figure()
    plt.scatter(co_op_samp.iloc[0:5000,0], co_op_samp.iloc[0:5000,1]*100, s = 75, c = 'white', edgecolor = 'k', linewidth = 0.25)
    plt.scatter(variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i])  & (var_scan_seqs.iloc[:,6] == 1)),'ACE'], variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i]) & (var_scan_seqs.iloc[:,6] == 1)),'pAb']*100, s = 125, c = lighten_color(variant_colors[i], 0.15), edgecolor = 'k', linewidth = 0.25)
    plt.scatter(ace_variant[i], pAb_variant[i]*100, c = variant_colors[i], s = 230, edgecolor = 'k', linewidth = 0.25)
    plt.scatter(max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'ACE']), variant_scan_predictions_both.loc[(variant_scan_predictions_both['ACE'] == max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'ACE'])),'pAb']*100, s = 240, c = lighten_color(variant_colors[i], 0.85), edgecolor = 'k', linewidth = 0.25, marker = '^')
    plt.scatter(variant_scan_predictions_both.loc[(variant_scan_predictions_both['pAb'] == max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'pAb'])),'ACE'], max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'pAb'])*100, s = 200, c = lighten_color(variant_colors[i], 0.85), edgecolor = 'k', linewidth = 0.25, marker = 's')
    plt.xlim(11.6, 12.05)
    plt.ylim(4.25, 10)
    plt.xticks([11.6, 11.7, 11.8, 11.9, 12.0], fontsize = 26)
    plt.yticks(fontsize = 26)


#%%
### delta inset
for i in np.arange(5,6):
    seq_ind = np.arange(0,20223,2889)
    variant_scan_predictions_sub = variant_scan_predictions.iloc[seq_ind[i]:seq_ind[i+1]]
    variant_scan_predictions_both = variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i])),:]
    plt.figure()
    plt.scatter(co_op_samp.iloc[0:5000,0], co_op_samp.iloc[0:5000,1]*100, s = 75, c = 'white', edgecolor = 'k', linewidth = 0.25)
    plt.scatter(variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i])  & (var_scan_seqs.iloc[:,6] == 1)),'ACE'], variant_scan_predictions_sub.loc[((variant_scan_predictions_sub.iloc[:,0] > ace_variant[i]) & (variant_scan_predictions_sub.iloc[:,1] > pAb_variant[i]) & (var_scan_seqs.iloc[:,6] == 1)),'pAb']*100, s = 125, c = lighten_color(variant_colors[i], 0.15), edgecolor = 'k', linewidth = 0.25)
    plt.scatter(ace_variant[i], pAb_variant[i]*100, c = variant_colors[i], s = 230, edgecolor = 'k', linewidth = 0.25)
    plt.scatter(max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'ACE']), variant_scan_predictions_both.loc[(variant_scan_predictions_both['ACE'] == max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'ACE'])),'pAb']*100, s = 240, c = lighten_color(variant_colors[i], 0.85), edgecolor = 'k', linewidth = 0.25, marker = '^')
    plt.scatter(variant_scan_predictions_both.loc[(variant_scan_predictions_both['pAb'] == max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'pAb'])),'ACE'], max(variant_scan_predictions_both.loc[(var_scan_seqs.iloc[:,6] == 1),'pAb'])*100, s = 200, c = lighten_color(variant_colors[i], 0.85), edgecolor = 'k', linewidth = 0.25, marker = 's')
    plt.xlim(11.6, 12.05)
    plt.ylim(4.25, 10)
    plt.xticks([11.6, 11.7, 11.8, 11.9, 12.0], fontsize = 26)
    plt.yticks(fontsize = 26)
    
    
    
    
    


