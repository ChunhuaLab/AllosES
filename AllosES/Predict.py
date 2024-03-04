# -*- coding: utf-8 -*-
"""
Created on Fri Mar 3 2023 15:44:20

@author: Fangrui Hu
"""
import pandas as pd
import numpy as np
import joblib
import warnings
warnings.filterwarnings('ignore')

def testdataPrepare(X,mean,std):
    nor_X = (X-mean)/std
    return nor_X

def predict_pro(model_path,data):
    final_proba = []
    for model in joblib.load(model_path):
        prepared_testset = testdataPrepare(data, model[1], model[2])
        testdata = np.array(prepared_testset)
        proba = model[0].predict_proba(np.array(testdata))
        final_proba.append(pd.DataFrame(proba)[1])
    final_proba = pd.concat(final_proba, axis=1)
    final_proba['pro_ave'] = final_proba.apply(lambda x: x.mean(), axis=1)  # Probability Mean
    return final_proba


