import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns

from distributed import Client, LocalCluster

import matplotlib.pyplot as plt

### title: Run scenic for each sample in sequential cutaneous melanoma dataset
### author: Yiping Wang date: 11/08/2022

def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]

if __name__=="__main__":

    #prefixes = ["r310_pre","r310_on","r310_on_later"]
    prefixes = ["r310_on","r310_on_later"]
    for aprefix in prefixes:
        tf_names = load_tf_names("hs_hgnc_curated_tfs.txt")
        dbs = [RankingDatabase(fname="hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather", name=name("hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"))]
        print("read1")
        df_r310 = pd.read_csv(aprefix+"_counts.csv")
        print("read2")
        client = Client(LocalCluster())
        adjacencies = grnboost2(df_r310.T, tf_names=tf_names, verbose=True)
        pickle.dump(adjacencies,open(aprefix+"_adjacencies.pickle","wb"))

        modules = list(modules_from_adjacencies(adjacencies, df_r310.T))
        df = prune2df(dbs, modules, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
        df.to_csv(aprefix+"_motifs.csv")
        regulons = df2regulons(df)
        pickle.dump(regulons, open(aprefix+"_regulons.pickle","wb"))
        auc_mtx = aucell(df_r310.T, regulons, num_workers=4)

        auc_mtx.to_csv(aprefix+"_auc_mtx.csv")
        
        plt.figure(figsize=(20,20))
        sns.clustermap(auc_mtx)
        plt.savefig(aprefix+"_auc_mtx_clustermap.pdf")
