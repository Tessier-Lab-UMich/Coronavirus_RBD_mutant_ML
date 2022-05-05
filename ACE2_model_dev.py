# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 12:53:48 2021

@author: makow
"""

from utils import *

ace_github = pd.read_csv(".\\binding_Kds.csv", header = 0, index_col = 0)
ace_lib1 = ace_github[ace_github.index == 'lib1']
ace_lib2 = ace_github[ace_github.index == 'lib2']
ace_binding = ace_binding_prepro(ace_github)
ace_binding_lib1 = ace_binding_prepro(ace_lib1)
ace_binding_lib2 = ace_binding_prepro(ace_lib2)
print(sc.stats.spearmanr(ace_binding.iloc[:,2], ace_binding.iloc[:,1]))

#%%
#tuning sequence observation cutoff
count_test_spearman_coefs = []
count_train_spearman_coefs = []
for j in np.arange(5,125,1):
    ace_binding = ace_binding[ace_binding[3] > j]
    ace_binding.reset_index(drop = True, inplace = True)
    ace_binding_ohe = ohe_encode(ace_binding)
    train_spearman = []
    test_spearman = []
    for i in np.arange(0,5):
        ace_binding_train, ace_binding_test, target_train, target_test = train_test_split(ace_binding_ohe, ace_binding.iloc[:,1])
        covid_ridge = Ridge(alpha = 1.9)
        covid_ridge.fit(ace_binding_train, target_train)
        covid_ridge_predict_train = pd.DataFrame(covid_ridge.predict(ace_binding_train))
        covid_ridge_predict_test = pd.DataFrame(covid_ridge.predict(ace_binding_test))
        train_spearman.append(sc.stats.spearmanr(covid_ridge_predict_train.iloc[:,0], target_train)[0])
        test_spearman.append(sc.stats.spearmanr(covid_ridge_predict_test.iloc[:,0], target_test)[0])
    count_train_spearman_coefs.append(np.mean(train_spearman))
    count_test_spearman_coefs.append(np.mean(test_spearman))


#tuning alpha
alpha_test_spearman_coefs = []
alpha_train_spearman_coefs = []
for j in np.arangenp.arange(0.1,10,0.1):
    ace_binding = ace_binding[ace_binding[3] > 53]
    ace_binding.reset_index(drop = True, inplace = True)
    ace_binding_ohe = ohe_encode(ace_binding)
    train_spearman = []
    test_spearman = []
    for i in np.arange(0,5):
        ace_binding_train, ace_binding_test, target_train, target_test = train_test_split(ace_binding_ohe, ace_binding.iloc[:,1])
        covid_ridge = Ridge(alpha = j)
        covid_ridge.fit(ace_binding_train, target_train)
        covid_ridge_predict_train = pd.DataFrame(covid_ridge.predict(ace_binding_train))
        covid_ridge_predict_test = pd.DataFrame(covid_ridge.predict(ace_binding_test))
        train_spearman.append(sc.stats.spearmanr(covid_ridge_predict_train.iloc[:,0], target_train)[0])
        test_spearman.append(sc.stats.spearmanr(covid_ridge_predict_test.iloc[:,0], target_test)[0])
    alpha_train_spearman_coefs.append(np.mean(train_spearman))
    alpha_test_spearman_coefs.append(np.mean(test_spearman))

#%%
#final model train/test
ace_binding = ace_binding[ace_binding[3] > 53]
ace_binding_ohe = ohe_encode(ace_binding)

ace_binding_lib1 = ace_binding_lib1[ace_binding_lib1[3] > 1]
ace_binding_lib1_ohe = ohe_encode(ace_binding_lib1)

ace_binding_lib2 = ace_binding_lib2[ace_binding_lib2[3] > 1]
ace_binding_lib2_ohe = ohe_encode(ace_binding_lib2)

#%%
ace_binding_train, ace_binding_test, target_train, target_test = train_test_split(ace_binding_ohe, ace_binding.iloc[:,1])
ace_ridge = Ridge(alpha = 1.90)
ace_ridge.fit(ace_binding_train, target_train)
ace_ridge_train_predict = pd.DataFrame(ace_ridge.predict(ace_binding_train))
ace_ridge_test_predict = pd.DataFrame(ace_ridge.predict(ace_binding_test))
print(sc.stats.spearmanr(ace_ridge_train_predict.iloc[:,0], target_train))
print(mean_absolute_error(ace_ridge_train_predict.iloc[:,0], target_train))
print(mean_absolute_percentage_error(ace_ridge_train_predict.iloc[:,0], target_train))
print(sc.stats.spearmanr(ace_ridge_test_predict.iloc[:,0], target_test))
print(mean_absolute_error(ace_ridge_test_predict.iloc[:,0], target_test))
print(mean_absolute_percentage_error(ace_ridge_test_predict.iloc[:,0], target_test))
ace_ridge_wt_predict = pd.DataFrame(ace_ridge.predict(wt_ohe))
ace_ridge_lib1_predict = pd.DataFrame(ace_ridge.predict(ace_binding_lib1_ohe))
ace_ridge_lib2_predict = pd.DataFrame(ace_ridge.predict(ace_binding_lib2_ohe))
print(sc.stats.spearmanr(ace_ridge_lib1_predict.iloc[:,0], ace_binding_lib1.iloc[:,1]))
print(sc.stats.spearmanr(ace_ridge_lib2_predict.iloc[:,0], ace_binding_lib2.iloc[:,1]))

plt.figure(0)
plt.scatter(ace_ridge_train_predict.iloc[0:5000,0], target_train[0:5000], c = 'blue', s = 75, edgecolor = 'k', linewidth = 0.25)
plt.scatter(ace_ridge_test_predict.iloc[0:2500,0], target_test[0:2500], c = 'darkgray', s = 75, edgecolor = 'k', linewidth = 0.25)
plt.scatter(ace_ridge_wt_predict.iloc[0,0], 10.79, c = 'red', s = 100, edgecolor = 'k', linewidth = 0.5)
plt.yticks(fontsize = 26)
plt.xlabel('')
plt.ylabel('')
plt.xticks(fontsize = 26)
plt.xlim(-0.25,12.25)

plt.figure(1)
plt.scatter(ace_ridge_lib1_predict.iloc[0:5000,0], ace_binding_lib1.iloc[0:5000,1], c = 'lightskyblue', s = 75, edgecolor = 'k', linewidth = 0.25)
plt.scatter(ace_ridge_wt_predict.iloc[0,0], 10.79, c = 'red', s = 100, edgecolor = 'k', linewidth = 0.5)
plt.yticks(fontsize = 26)
plt.xlabel('')
plt.ylabel('')
plt.xticks(fontsize = 26)
plt.xlim(-0.25,12.25)

plt.figure(2)
plt.scatter(ace_ridge_lib2_predict.iloc[0:5000,0], ace_binding_lib2.iloc[0:5000,1], c = 'lightskyblue', s = 75, edgecolor = 'k', linewidth = 0.25)
plt.scatter(ace_ridge_wt_predict.iloc[0,0], 10.79, c = 'red', s = 100, edgecolor = 'k', linewidth = 0.5)
plt.yticks(fontsize = 26)
plt.xlabel('')
plt.ylabel('')
plt.xticks(fontsize = 26)
plt.xlim(-0.25,12.25)

    
    
