#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 10:17:56 2019

@author: Sandro-Hurtado
"""
import os
from scipy import stats as sst
import pandas as pd
import pre_processing.find_stable as fs
import pre_processing.pre_analisis as pa
import pre_processing.mongo_connection as mongo

# path = '/home/antonio/Anaconda_Projects/Python_Projects/Gene_Regulation_Networks'
# os.chdir(path)

"""
This class preprocess the data and create a geneExpression dataframe in a proper form to work with 
"""


class GeneExpressions:

    def __init__(self, ids, labels, percentage):
        self.id_list = ids
        self.label_list = labels
        self.percentage = percentage
        self.dfz = pd.DataFrame()

    def load_gene_expressions(self):
        os.getcwd()
        res = pa.run_norm_mongodb(self.id_list)
        cont = 0
        df = pd.DataFrame(index=res[cont].tolist())

        # If no option label is selected the name of the sample will be displayed
        if len(self.label_list) < 2:
            self.label_list = []
            for id_sample in self.id_list:
                mongodb = mongo.MongoConnection()
                label = mongodb.get_sample_name(id_sample)
                self.label_list.append(label)

        for i in self.label_list:
            cont += 1
            df[i] = res[cont].tolist()

        # obtain unstable genes
        n_genes = int(df.index.size * self.percentage)

        unstable_genes = fs.findUnstable(df.values, n_genes)

        # select unstable genes from gene expression to plot
        M = pd.DataFrame()
        for g in unstable_genes:
            M = pd.concat([M, df.iloc[[g]]])
        # transform to zscope in axis 1 (rows)
        dfzs1 = sst.zscore(M.values, axis=1, ddof=1)
        self.dfz = pd.DataFrame(index=M.index, data=dfzs1, columns=self.label_list)
        return self.dfz
