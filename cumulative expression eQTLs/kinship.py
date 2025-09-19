import numpy as np
import pandas as pd
from pandas_plink import read_plink

def load_genotypes_plink(plink_prefix):
    """
    Load PLINK bed/bim/fam and return numpy genotype matrix.

    Parameters
    ----------
    plink_prefix : str
        Prefix of PLINK binary files (.bed/.bim/.fam)

    Returns
    -------
    bim : DataFrame (variant info)
    fam : DataFrame (sample info)
    G   : dask.array (raw genotype matrix)
    genotypes : ndarray (samples × snps) with missing filled as 2.0
    """
    bim, fam, G = read_plink(plink_prefix, verbose=False)

    # G: dask array [variants × samples]; transpose → [samples × variants]
    bed = G.compute().T

    # Replace missing (NaN) with 2.0
    genotypes = np.where(np.isnan(bed), 2.0, bed)

    return bim, fam, G, genotypes

def generate_kinship(genotypes):
    """
    Compute kinship matrix from genotype matrix.

    Parameters
    ----------
    genotypes : ndarray (samples × snps)

    Returns
    -------
    kinship : ndarray (samples × samples)
    """
    # Standardize genotypes (per SNP)
    G_std = (genotypes - np.nanmean(genotypes, axis=0)) / np.nanstd(genotypes, axis=0)
    G_std = np.nan_to_num(G_std)

    # Kinship = G * G^T / M  (M = number of SNPs)
    M = G_std.shape[1]
    kinship = (G_std @ G_std.T) / M
    return kinship


if __name__ == "__main__":
    geno_prefix = "/path/to/merged-vcf"
    output_filename = geno_prefix + "_kinship.txt"

    bim, fam, bed, genotype_mat = load_genotypes_plink(geno_prefix)
    kinship_mat = generate_kinship(genotype_mat)

    kinship_df = pd.DataFrame(data=kinship_mat, index=fam["iid"], columns=fam["iid"])

