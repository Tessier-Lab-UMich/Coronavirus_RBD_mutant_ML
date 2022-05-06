# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 14:40:23 2021

@author: makow
"""

from utils import *

#ace data import and training
ace_github = pd.read_csv(".\\binding_Kds.csv", header = 0, index_col = 0)
ace_binding = ace_binding_prepro(ace_github)
ace_binding = ace_binding[ace_binding[3] > 53]
ace_binding.reset_index(drop = True, inplace = True)
ace_binding_ohe = ohe_encode(ace_binding)

#ace model training
ace_ridge = Ridge(alpha = 1.90)
ace_ridge.fit(ace_binding_ohe, ace_binding.iloc[:,1])
ace_binding_predict = pd.DataFrame(ace_ridge.predict(ace_binding_ohe))

#pAb data import and processing
pAb_escape = pd.read_csv(".\\scores.txt", header = 0, index_col = 0)
pAb_escape = pAb_escape[pAb_escape.iloc[:,3] > 15]
pAb_escape.reset_index(drop = True, inplace = True)
pAb_escape_ohe = ohe_encode(pAb_escape)

#pAb model trianing
pAb_ridge = Ridge(alpha = 6.2)
pAb_ridge.fit(pAb_escape_ohe, pAb_escape.iloc[:,1])
pAb_escape_predict = pd.DataFrame(pAb_ridge.predict(pAb_escape_ohe))

#variant one hot encoding and model prediction
variant_ohe = ohe_encode(variant_seqs)
ace_variant = ace_ridge.predict(variant_ohe)
pAb_variant = pAb_ridge.predict(variant_ohe)


#%%
mut_count_ace = pd.read_csv(".\\mut_count_ace.csv", header = 0, index_col = 0)
mut_count_ace.columns = np.arange(0,201)
mut_count_pAb = pd.read_csv(".\\mut_count_pAb.csv", header = 0, index_col = 0)
mut_count_pAb.columns = np.arange(0,201)

sites = np.arange(331,532,1)
sampling = pd.DataFrame(index = sites, columns = sites)
for i in sites:
    for ii in sites:
        if i !=ii:
            sum_res1 = sum(mut_count_ace[i-331] > 6)
            sum_res2 = sum(mut_count_ace[ii-331] > 6)
            sampling.loc[i,ii] = (sum_res1*sum_res2)


wt_scan_seqs = []
sites = np.arange(0,201,1)
for i in sites:
    for ii in sites:
        if i !=ii:
            for j in alph_letters:
                for jj in alph_letters:
                    wt_mut = list(wt_seq).copy()
                    res1 = wt_mut[i]
                    res2 = wt_mut[ii]
                    wt_mut[i] = j
                    wt_mut[ii] = jj
                    if res1 != j:
                        if res2 != jj:
                            ham = sc.spatial.distance.hamming(wt_mut, list(wt_seq))
                            if ham > 0.005:
                                if mut_count_ace.loc[j,i] > 6:
                                    if mut_count_ace.loc[jj,ii] > 6:
                                        if mut_count_pAb.loc[j,i] > 6:
                                            if mut_count_pAb.loc[jj,ii] > 6:
                                                wt_scan_seqs.append([''.join(wt_mut), i,res1,j, ii,res2,jj])
wt_scan_seqs = pd.DataFrame(wt_scan_seqs)

#%%
wt_scan = pd.DataFrame()
wt_scan_predictions = pd.DataFrame()
for i in np.arange(0,8300000,100000):
    wt_scan_seqs_used = wt_scan_seqs.iloc[i:i+100000,:]
    wt_scan_ohe = ohe_encode(wt_scan_seqs_used)
    wt_scan_ace_predict = pd.DataFrame(ace_ridge.predict(wt_scan_ohe))
    wt_scan_pAb_predict = pd.DataFrame(pAb_ridge.predict(wt_scan_ohe))
    wt_predictions = pd.concat([wt_scan_ace_predict, wt_scan_pAb_predict], axis = 1)
    wt_predictions.index = wt_scan_seqs_used.index
    wt_scan = pd.concat([wt_scan, wt_scan_seqs_used], axis = 0)
    wt_scan_predictions = pd.concat([wt_scan_predictions, wt_predictions], axis = 0)

wt_scan_predictions.index = wt_scan_seqs.index
wt_scan_predictions_samp = wt_scan_predictions.sample(5000)
wt_scan_predictions_samp.columns = ['ACE', 'pAb']

#%%
#isolating variants in the 'Region of concern'
concern_region_predictions = wt_scan_predictions.loc[(wt_scan_predictions.iloc[:,0] > ace_variant[2]) & (wt_scan_predictions.iloc[:,1] > 0),:]
concern_region = wt_scan_seqs.loc[(wt_scan_predictions.iloc[:,0] > ace_variant[2]) & (wt_scan_predictions.iloc[:,1] > 0),:]
concern_region_predictions.columns = [0,1]
concern_region.drop_duplicates(subset=0, keep='first', inplace = True)
concern_region_predictions.drop_duplicates(subset=0, keep='first', inplace = True)
concern_region_predictions_samp = concern_region_predictions.sample(1000)

single_nuc = []
for index, row in concern_region.iterrows():
    orig_dna = str(''.join(wt_dna[(row[1]*3):((row[1]*3)+3)]))
    sub_dna = dna_table[row[3]]
    hams = []
    for i in sub_dna:
        hams.append(sc.spatial.distance.hamming(list(i), list(orig_dna))*3)
    ham1 = min(hams)
    hams = []
    orig_dna = str(''.join(wt_dna[(row[4]*3):((row[4]*3)+3)]))
    sub_dna = dna_table[row[6]]
    for i in sub_dna:
        hams.append(sc.spatial.distance.hamming(list(i), list(orig_dna))*3)
    single_nuc.append([ham1, min(hams), ham1+min(hams)])
single_nuc = pd.DataFrame(single_nuc)
single_nuc.set_index(concern_region.index, inplace = True)

#%%
x = [ace_variant[2], 12.25]
y1 = [18,18]
y2 = [0,0]

plt.figure()
plt.fill_between(x, y1, y2, color='darkgray', alpha=0.25, edgecolor = 'k', linewidth = 0.25)
plt.scatter(wt_scan_predictions_samp.loc[(wt_scan_predictions_samp.iloc[:,0] < ace_variant[2]) | (wt_scan_predictions_samp.iloc[:,1] < 0),'ACE'], wt_scan_predictions_samp.loc[(wt_scan_predictions_samp.iloc[:,0] < ace_variant[2]) | (wt_scan_predictions_samp.iloc[:,1] < 0),'pAb']*100, s = 75, c = 'white', edgecolor = 'k', linewidth = 0.25)
plt.scatter(concern_region_predictions_samp.iloc[:,0], concern_region_predictions_samp.iloc[:,1]*100, s = 75, c = 'dimgray', edgecolor = 'k', linewidth = 0.25)

plt.scatter(ace_variant, pAb_variant*100, c = np.arange(0,7), cmap = cmap_var, s = 150, edgecolor = 'k', linewidth = 0.25)
plt.ylim(-2, 17)
plt.xlim(10.25, 12.25)
plt.xticks([10.5, 11.0, 11.5, 12.0], fontsize = 26)
plt.yticks([0, 5, 10, 15], fontsize = 26)


#%%
#isolating variants at the Pareto frontier (Figure 3A)
pareto_predictions = wt_scan_predictions.loc[(((wt_scan_predictions.iloc[:,1])+(0.1186*(wt_scan_predictions.iloc[:,0]))>1.43872) & ((wt_scan_predictions.iloc[:,0])>ace_variant[2])),:]
pareto = wt_scan_seqs.loc[(((wt_scan_predictions.iloc[:,1])+(0.1186*(wt_scan_predictions.iloc[:,0]))>1.43872) & ((wt_scan_predictions.iloc[:,0])>ace_variant[2])),:]

pareto_predictions.columns = [0,1]
pareto_predictions.drop_duplicates(keep = 'first', inplace = True)
pareto.drop_duplicates(subset = 0, keep = 'first', inplace = True)

single_nuc = []
for index, row in pareto.iterrows():
    orig_dna = str(''.join(wt_dna[(row[1]*3):((row[1]*3)+3)]))
    sub_dna = dna_table[row[3]]
    hams = []
    for i in sub_dna:
        hams.append(sc.spatial.distance.hamming(list(i), list(orig_dna))*3)
    ham1 = min(hams)
    hams = []
    orig_dna = str(''.join(wt_dna[(row[4]*3):((row[4]*3)+3)]))
    sub_dna = dna_table[row[6]]
    for i in sub_dna:
        hams.append(sc.spatial.distance.hamming(list(i), list(orig_dna))*3)
    single_nuc.append([ham1, min(hams), ham1+min(hams)])
single_nuc = pd.DataFrame(single_nuc)
single_nuc.set_index(pareto.index, inplace = True)

### need to add every site two before an S or a T and remove it ist mutated to N
glycan = []
for index, row in pareto.iterrows():
    if np.isin(row[1], glycan_begin):
        glycan.append(1)
    elif np.isin(row[4], glycan_begin):
        glycan.append(1)
    elif (np.isin(row[1], glycan_end_st)) & (row[3] == 'T'):
        glycan.append(1)
    elif (np.isin(row[1], glycan_end_st)) & (row[3] == 'S'):
        glycan.append(1)
    elif (np.isin(row[4], glycan_end_st)) & (row[6] == 'T'):
        glycan.append(1)
    elif (np.isin(row[4], glycan_end_st)) & (row[6] == 'S'):
        glycan.append(1)
    elif (np.isin(row[1], glycan_end_n)) & (row[3] == 'N'):
        glycan.append(1)
    elif (np.isin(row[4], glycan_end_n)) & (row[6] == 'N'):
        glycan.append(1)
    else:
        glycan.append(0)
glycan = pd.DataFrame(glycan)
glycan.index = pareto.index

pareto_predictions = pd.concat([pareto_predictions, single_nuc, glycan], axis = 1, ignore_index = False)
pareto_predictions.columns = np.arange(0,6)
pareto = pareto.loc[(pareto_predictions.iloc[:,4]==2) & (pareto_predictions.iloc[:,5]==0)]
pareto_predictions = pareto_predictions.loc[(pareto_predictions.iloc[:,4]==2)  & (pareto_predictions.iloc[:,5]==0)]

x2 = [ace_variant[2],12.1, 12.25]
y5 = [18,0,0]
y6 = [17.1,17.1,17.1]

plt.figure()
plt.fill_between(x, y1, y2, color='darkgray', alpha=0.25, edgecolor = 'k', linewidth = 0.25)
plt.fill_between(x2, y5, y6, color='darkgray', alpha=0.5, edgecolor = 'k', linewidth = 0.25)
plt.scatter(wt_scan_predictions_samp.iloc[:,0], wt_scan_predictions_samp.iloc[:,1]*100, s = 75, c = 'white', edgecolor = 'k', linewidth = 0.25)
plt.scatter(pareto_predictions.iloc[:,0], pareto_predictions.iloc[:,1]*100, c = 'dimgray', edgecolor = 'k', linewidth = 0.25, s = 75)
plt.scatter(ace_variant[0], pAb_variant[0]*100, c = variant_colors[0], s = 150, edgecolor = 'k', linewidth = 0.25)

plt.ylim(-2, 17)
plt.xlim(10.25, 12.25)
plt.xticks([10.5, 11.0, 11.5, 12.0], fontsize = 26)
plt.yticks([0, 5, 10, 15], fontsize = 26)

co_op_samp = wt_scan_predictions.sample(5000)


#%%
colormap12 = ['white', 'blue']
cmap12 = LinearSegmentedColormap.from_list("mycmap", colormap12)
extreme_mut_heatmap_ace = pd.read_csv(".\\8.20.21_heatmap_muts_ace2.csv", index_col = 0, header = 0)
plt.figure(figsize = (12,4))
sns.heatmap(extreme_mut_heatmap_ace, cmap = 'bwr', annot = True, linewidths = 0.1, linecolor = 'silver', annot_kws = {'fontsize': 13}, fmt = '.1f', vmin = 10.15)
plt.xticks(rotation = 90, fontsize = 22)
plt.yticks(fontsize = 22)


extreme_mut_heatmap_pAb = pd.read_csv(".\\8.20.21_heatmap_muts_pAb.csv", index_col = 0, header = 0)
plt.figure(figsize = (12,4))
sns.heatmap(extreme_mut_heatmap_pAb, cmap = 'bwr', annot = True, linewidths = 0.1, linecolor = 'silver', annot_kws = {'fontsize': 13}, fmt = '.1f', vmin = -14, vmax = 20)
plt.xticks(rotation = 90, fontsize = 22)
plt.yticks(fontsize = 22)


