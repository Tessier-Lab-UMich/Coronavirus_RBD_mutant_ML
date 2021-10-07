# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 12:54:17 2021

@author: makow
"""

from utils import *

pAb_github = pd.read_csv(".\\scores.csv", header = 0, index_col = 0)
pAb_escape = pAb_escape_prepro(pAb_github)

#tuning sequence observation cutoff
count_test_spearman_coefs = []
count_train_spearman_coefs = []
for j in np.arange(5,50,1):
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


#tuning alpha
alpha_test_spearman_coefs = []
alpha_train_spearman_coefs = []
for j in np.arangenp.arange(0.1,10,0.1):
    pAb_escape = pAb_escape[pAb_escape[3] > 12]
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

#%%
#final model train/test
pAb_ridge = Ridge(alpha = 12.07)
pAb_ridge.fit(pAb_escape_train, target_train)
pAb_ridge_train_predict = pd.DataFrame(pAb_ridge.predict(pAb_escape_train))
pAb_ridge_test_predict = pd.DataFrame(pAb_ridge.predict(pAb_escape_test))
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

