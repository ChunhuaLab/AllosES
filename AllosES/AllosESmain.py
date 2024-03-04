# -*- coding: utf-8 -*-
"""
Created on Fri Mar 3 2023 17:13:13

@author: Fangrui Hu
"""
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from Predict import predict_pro
from utils import AllosESutils,parser
if __name__ == "__main__":
    args = parser.parse_args()
    PDBID = args.PDBID
    chain = args.CHAIN
    path = r'./'
    original_features = AllosESutils(path, PDBID,chain)
    feature = original_features.original_features()
    model_path = path + 'models.m'
    data = feature.drop(columns=['residues', 'pocket_name']).astype(float)
    final_proba = predict_pro(model_path,data)
    pocket_res = pd.concat([final_proba, feature['residues']],axis=1)
    result = pocket_res.sort_values(by='pro_ave', ascending=False, axis=0)[['residues','pro_ave']]
    result.to_csv(path + PDBID.upper() + '_' + chain.upper() + '_result.csv', index=True)
    print('******Start prediction...******')
    print('******Prediction completed!******')
    print('→The Pocket Rankings:')
    indeval = 0
    for inde in result['residues']:
        indeval += 1
        print('  Prediction Pocket' + str(indeval) + ':\t' + inde)