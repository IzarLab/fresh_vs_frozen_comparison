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

if __name__=="__main__":

    tf_names = load_tf_names("hs_hgnc_curated_tfs.txt")
    dbs = [RankingDatabase(fname="motifs-v9-nr.hgnc-m0.001.o0.0.tbl", name=name("motifs-v9-nr.hgnc-m0.001-o0.0.tbl"))]
    df_r310_pre = pd.read_csv("/data/r310_pre_counts.csv")
    client = Client(LocalCluster())
    adjacencies = grnboost2(df_r310_pre.T, tf_names=tf_names, verbose=True)
