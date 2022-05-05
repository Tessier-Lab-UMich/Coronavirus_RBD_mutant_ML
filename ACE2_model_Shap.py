# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 12:53:48 2021

@author: makow
"""

from utils import *
import shap
import sklearn.metrics as metrics
from sklearn.ensemble import RandomForestClassifier as RFC
import xgboost
import pycircos
import matplotlib.pyplot as plt
import collections

ace_github = pd.read_csv(".\\binding_Kds.csv", header = 0, index_col = 0)
ace_lib1 = ace_github[ace_github.index == 'lib1']
ace_lib2 = ace_github[ace_github.index == 'lib2']
ace_binding = ace_binding_prepro(ace_github)
ace_binding = ace_binding[ace_binding[3] > 100]
ace_binding_ohe = ohe_encode(ace_binding)
ace_binding['label'] = 0
for index, row in ace_binding.iterrows():
    if row[1] > 10.5:
        ace_binding.loc[index, 'label'] = 1
ace_ohe_train, ace_ohe_test, ace_ohe_target_train, ace_ohe_target_test = train_test_split(ace_binding_ohe, ace_binding['label'])

#%%
rfc_ace = xgboost.XGBClassifier(n_estimators = 50, max_depth = 50)
#cv_results = cv(rfc_ace, ace_binding_ohe, ace_binding.iloc[:,1], cv = 2, return_train_score = True)
#print(np.mean(cv_results['train_score']))
#print(np.mean(cv_results['test_score']))
rfc_ace.fit(ace_ohe_train, ace_ohe_target_train)
rfc_ace_predict_proba = rfc_ace.predict_proba(ace_ohe_test)
rfc_ace_predict = rfc_ace.predict(ace_ohe_test)

ace_fpr, ace_tpr, threshold = metrics.roc_curve(ace_ohe_target_test, rfc_ace_predict)
ace_roc_auc = metrics.auc(ace_fpr, ace_tpr)
plt.plot(ace_fpr, ace_tpr, 'b', label = 'ACE AUC = %0.3f' % ace_roc_auc, linewidth = 3)
plt.legend(loc = 'lower right', fontsize = 16)
plt.plot([0, 1], [0, 1],'k-')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.ylabel(' ', fontsize = 22)
plt.xlabel(' ', fontsize = 22)
plt.xticks(fontsize = 26)
plt.yticks(fontsize = 26)
plt.tight_layout()



