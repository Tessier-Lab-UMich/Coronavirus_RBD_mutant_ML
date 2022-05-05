# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 12:54:17 2021

@author: makow
"""

from utils import *


"""
#data processing if scores are imported from Bloom Github
pAb_github = pd.read_csv("...\\GitHub\\SARS-CoV-2-RBD_MAP_HAARVI_sera\\results\\escape_scores\\scores.csv", header = 0, index_col = 0)
pat_name = ['12C', '13_', '1C_', '22C', '23C', '23_', '24C', '25C', '25_', '6C_', '7C_']

pAb_github['pat'] = pAb_github.index.str[:3]

pAb_pat = []
for i in pat_name:
    pAb_pat_i = pAb_github[pAb_github['pat'] == i]
    pAb_pat.append(pAb_escape_prepro(pAb_pat_i))
pAb_escape = pAb_escape_prepro(pAb_github[pAb_github['pat'].isin(pat_name)])
"""

pAb_escape = pd.read_csv(".\\scores.txt", header = 0, index_col = 0)
print(sc.stats.spearmanr(pAb_escape.iloc[:,2], 100*pAb_escape.iloc[:,1]))

#%%
#tuning sequence observation cutoff
count_test_spearman_coefs = []
count_train_spearman_coefs = []
for j in np.arange(1,30,0.5):
    pAb_escape = pAb_escape[pAb_escape[3] > j]
    pAb_escape.reset_index(drop = True, inplace = True)
    pAb_escape_ohe = ohe_encode(pAb_escape)
    train_spearman = []
    test_spearman = []
    for i in np.arange(0,5):
        pAb_escape_train, pAb_escape_test, target_train, target_test = train_test_split(pAb_escape_ohe, pAb_escape.iloc[:,1])
        covid_ridge = Ridge(alpha = 1)
        covid_ridge.fit(pAb_escape_train, target_train)
        covid_ridge_predict_train = pd.DataFrame(covid_ridge.predict(pAb_escape_train))
        covid_ridge_predict_test = pd.DataFrame(covid_ridge.predict(pAb_escape_test))
        train_spearman.append(sc.stats.spearmanr(covid_ridge_predict_train.iloc[:,0], target_train)[0])
        test_spearman.append(sc.stats.spearmanr(covid_ridge_predict_test.iloc[:,0], target_test)[0])
    count_train_spearman_coefs.append(np.mean(train_spearman))
    count_test_spearman_coefs.append(np.mean(test_spearman))

plt.scatter(np.arange(1,30,0.5), count_test_spearman_coefs, c = 'gray', s = 65)
plt.scatter(np.arange(1,30,0.5), count_train_spearman_coefs,c = 'blue', s = 65)
plt.xticks(fontsize = 16)
plt.yticks([0.68, 0.72, 0.76, 0.80, 0.84], fontsize = 16)
#plt.ylim(0.66, 0.86)

#%%
#tuning alpha
alpha_test_spearman_coefs = []
alpha_train_spearman_coefs = []
for j in np.arange(0.2,20,0.2):
    pAb_escape = pAb_escape[pAb_escape[3] > 15]
    pAb_escape.reset_index(drop = True, inplace = True)
    pAb_escape_ohe = ohe_encode(pAb_escape)
    train_spearman = []
    test_spearman = []
    for i in np.arange(0,5):
        pAb_escape_train, pAb_escape_test, target_train, target_test = train_test_split(pAb_escape_ohe, pAb_escape.iloc[:,1])
        covid_ridge = Ridge(alpha = j)
        covid_ridge.fit(pAb_escape_train, target_train)
        covid_ridge_predict_train = pd.DataFrame(covid_ridge.predict(pAb_escape_train))
        covid_ridge_predict_test = pd.DataFrame(covid_ridge.predict(pAb_escape_test))
        train_spearman.append(sc.stats.spearmanr(covid_ridge_predict_train.iloc[:,0], target_train)[0])
        test_spearman.append(sc.stats.spearmanr(covid_ridge_predict_test.iloc[:,0], target_test)[0])
    alpha_train_spearman_coefs.append(np.mean(train_spearman))
    alpha_test_spearman_coefs.append(np.mean(test_spearman))

plt.scatter(np.arange(0.2,20,0.2), alpha_test_spearman_coefs, c = 'gray', s = 65)
plt.scatter(np.arange(0.2,20,0.2), alpha_train_spearman_coefs, c = 'blue', s = 65)
plt.xticks(fontsize = 16)
plt.yticks([0.68, 0.72, 0.76, 0.80, 0.84], fontsize = 16)

#%%
#final model train/test
pAb_escape = pAb_escape[pAb_escape[3] > 15]
pAb_escape_ohe = ohe_encode(pAb_escape)
pAb_ridge = Ridge(alpha = 6.2)
pAb_escape_train, pAb_escape_test, target_train, target_test = train_test_split(pAb_escape_ohe, pAb_escape.iloc[:,1])
pAb_ridge.fit(pAb_escape_train, target_train)
pAb_ridge_train_predict = pd.DataFrame(pAb_ridge.predict(pAb_escape_train))
pAb_ridge_test_predict = pd.DataFrame(pAb_ridge.predict(pAb_escape_test))
print(sc.stats.spearmanr(pAb_ridge_train_predict.iloc[:,0], target_train))
print(mean_absolute_error(pAb_ridge_train_predict.iloc[:,0], target_train))
print(mean_absolute_percentage_error(pAb_ridge_train_predict.iloc[:,0], target_train))
print(sc.stats.spearmanr(pAb_ridge_test_predict.iloc[:,0], target_test))
print(mean_absolute_error(pAb_ridge_test_predict.iloc[:,0], target_test))
print(mean_absolute_percentage_error(pAb_ridge_test_predict.iloc[:,0], target_test))
pAb_ridge_wt_predict = pd.DataFrame(pAb_ridge.predict(wt_ohe))

plt.figure(0)
plt.scatter(pAb_ridge_train_predict.iloc[0:5000,0]*100, target_train[0:5000]*100, c = 'blue', s = 75, edgecolor = 'k', linewidth = 0.25)
plt.scatter(pAb_ridge_test_predict.iloc[0:2500,0]*100, target_test[0:2500]*100, c = 'darkgray', s = 75, edgecolor = 'k', linewidth = 0.25)
plt.scatter(pAb_ridge_wt_predict.iloc[0,0]*100, 0, c = 'red', s = 100, edgecolor = 'k', linewidth = 0.5)
plt.yticks(fontsize = 26)
plt.xticks(fontsize = 26)
plt.xlabel('')
plt.ylabel('')
plt.ylim(-2, 62)
plt.xlim(-2, 62)

#%%
"""
pAb_pat_filt = []
for i in pAb_pat:
    pAb_pat_filt.append(i[i.iloc[:,3]>0])

pAb_pat_ohe = []
for i in pAb_pat_filt:
    pAb_pat_ohe.append(ohe_encode(i))

pat_spearmans = []
for i in np.arange(11):
    pAb_pat_predict = pd.DataFrame(pAb_ridge.predict(pAb_pat_ohe[i]))
    pat_spearmans.append(sc.stats.spearmanr(pAb_pat_predict.iloc[:,0], pAb_pat_filt[i].iloc[:,1])[0])
    plt.figure()
    plt.scatter(pAb_pat_predict.iloc[0:2500,0]*100, pAb_pat_filt[i].iloc[0:2500,1]*100, c = 'lightskyblue', s = 75, edgecolor = 'k', linewidth = 0.25)
    plt.scatter(pAb_ridge_wt_predict.iloc[0,0]*100, 0, c = 'red', s = 100, edgecolor = 'k', linewidth = 0.5)
    plt.yticks(fontsize = 26)
    plt.xticks(fontsize = 26)
    plt.xlabel('')
    plt.ylabel('')
    plt.ylim(-2, 62)
    plt.xlim(-2, 62)

data = pd.read_csv("C:\\Users\\makow\\Documents\\Research\\Manuscripts\\2021 Covid Mutation Analysis\\PloS Comp Bio Submission\\pAb_subject.csv", header = 0, index_col = None)

plt.figure(figsize = (8,2))
sns.heatmap(data.iloc[:,0:1].T, cmap = 'bwr', square = True, vmin = 0.3, vmax = 0.6)

"""

