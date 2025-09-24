import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.decomposition import PCA
import numpy as np
from pyensembl import EnsemblRelease
import matplotlib.pyplot as plt

D11 = sc.read_h5ad('/path/to/D11.h5') 
D30 = sc.read_h5ad('/path/to/D30.h5') 
D52 = sc.read_h5ad('/path/to/D52.h5') 

D11_FPP = D11[(D11.obs['celltype'].isin(['FPP', 'NB'])) & (D11.obs['treatment'] == 'NONE')].copy()
D30_SERT = D30[(D30.obs['celltype'].isin(['FPP', 'Sert'])) & (D30.obs['treatment'] == 'NONE')].copy()
D52_SERT = D52[(D52.obs['celltype'].isin(['P_Sert', 'Sert'])) & (D52.obs['treatment'] == 'NONE')].copy()

# read donor list - open sourced
donor_df = pd.read_csv(
    "/path/to/opensoruce_donorlist.txt",
    delim_whitespace=True,
    dtype=str
)
open_donors = set(donor_df["donor_id"].tolist())

# find donors that are present in all three AnnData objects
common_donors = (
    set(D11_FPP.obs["donor_id"]) &
    set(D30_SERT.obs["donor_id"]) &
    set(D52_SERT.obs["donor_id"])
)

# restrict to donors that are both in the donorlist AND in all 3 datasets
final_donors = open_donors & common_donors

# filter each dataset
D11 = D11_FPP[D11_FPP.obs["donor_id"].isin(final_donors)].copy()
D30 = D30_SERT[D30_SERT.obs["donor_id"].isin(final_donors)].copy()
D52 = D52_SERT[D52_SERT.obs["donor_id"].isin(final_donors)].copy()


adata_combined = ad.concat([D11, D30, D52], 
                           join="outer",       # keep all genes
                           label="batch",      # new column in .obs called "batch"
                           keys=["D11","D30T","D52"], 
                           index_unique="-")  # keep original cell barcodes
adata_combined.obs_names_make_unique()
