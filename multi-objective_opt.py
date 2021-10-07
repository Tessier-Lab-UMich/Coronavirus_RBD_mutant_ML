# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 13:48:18 2021

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
pAb_github = pd.read_csv(".\\scores.txt", header = 0, index_col = 0)
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


#co-optimization dataframe
ace_binding_predict.index = ace_binding.iloc[:,0]
pAb_escape_predict['mut_num'] = pAb_escape.iloc[:,2]
pAb_escape_predict.index = pAb_escape.iloc[:,0]

co_op = pd.concat([ace_binding_predict, pAb_escape_predict], axis = 1, ignore_index = False)
co_op.dropna(inplace = True)
co_op.columns= ['ACE2', 'pAb Escape', 'Mut Num']

co_op_samp = co_op.sample(10000)

#%%
#co-optimization figure
plt.figure(0)
plt.scatter(co_op.loc[co_op['Mut Num'] > 2, 'ACE2'], co_op.loc[co_op['Mut Num'] > 2, 'pAb Escape']*100, s = 75, c = 'blue', edgecolor = 'k', linewidth = 0.25)
plt.scatter(co_op_samp.loc[co_op_samp['Mut Num'] == 2, 'ACE2'], co_op_samp.loc[co_op_samp['Mut Num'] == 2, 'pAb Escape']*100, s = 75, c = 'darkgray', edgecolor = 'k', linewidth = 0.25)

plt.scatter(ace_variant, pAb_variant*100, c = np.arange(0,6), cmap = cmap_var, s = 150, edgecolor = 'k', linewidth = 0.25)
plt.xlim(5.8, 13)
plt.ylim(-5.2, 45)
plt.xticks([6, 8, 10, 12], fontsize = 26)
plt.yticks([0, 10, 20, 30, 40], fontsize = 26)


#%%
#region of concern figure
x = [ace_variant[2],13]
y1 = [20,20]
y2 = [0,0]

plt.figure()
plt.fill_between(x, y1, y2, color='darkgray', alpha=0.25, edgecolor = 'k', linewidth = 0.25)
plt.scatter(co_op.loc[co_op['Mut Num'] > 2, 'ACE2'], co_op.loc[co_op['Mut Num'] > 2, 'pAb Escape']*100, s = 75, c = 'blue', edgecolor = 'k', linewidth = 0.25)
plt.scatter(co_op_samp.loc[co_op_samp['Mut Num'] == 2, 'ACE2'], co_op_samp.loc[co_op_samp['Mut Num'] == 2, 'pAb Escape']*100, s = 75, c = 'darkgray', edgecolor = 'k', linewidth = 0.25)
plt.scatter(ace_variant, pAb_variant*100, c = np.arange(0,6), cmap = cmap_var, s = 150, edgecolor = 'k', linewidth = 0.25)
plt.ylim(-2, 17)
plt.xlim(10.25, 12.25)
plt.xticks([10.5, 11.0, 11.5, 12.0], fontsize = 26)
plt.yticks([0, 5, 10, 15], fontsize = 26)




