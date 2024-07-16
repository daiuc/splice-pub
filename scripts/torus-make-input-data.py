def compute_t_statistic(slope, p_value, df):
    # Calculate the two-tailed critical value
    critical_value = stats.t.ppf((1 + (1 - p_value)) / 2, df)

    # Calculate the t-statistic
    t_statistic = np.sign(slope) * abs(critical_value)

    return t_statistic


if __name__ == "__main__":

    import argparse
    import gzip

    import numpy as np
    import pandas as pd
    from scipy import stats

    opts = argparse.ArgumentParser(description="Generate input data for Torus")

    opts.add_argument("--tissue", type=str, required=True, help="Tissue name")
    opts.add_argument("--fdr", type=float, default=0.1, help="FDR threshold")
    opts.add_argument("--outprefix", type=str, default="torus", help="Output prefix")

    opt = opts.parse_args()
    tissue, fdr, outprefix = opt.tissue, opt.fdr, opt.outprefix

    print(f"processing {tissue}...")

    # Load the data
    sgene_df = []
    for i in range(1, 23):
        sgene_file = f"/project/yangili1/cdai/SpliFi/code/results/qtl/noisy/GTEx/{tissue}/separateNoise/cis_100000/perm/chr{str(i)}.addQval.txt.gz"
        df = pd.read_csv(sgene_file, sep=" ")
        df = df[(df["q"] < 0.1) & (df["phenotype_id"].str.contains("PR|UP"))]
        sgene_df.append(df)
    sgene_df = pd.concat(sgene_df)

    # Compute the t-statistic
    sgene_df["t_stat"] = [
        compute_t_statistic(s, p, d)
        for s, p, d in zip(
            sgene_df["slope"], sgene_df["pval_nom"], sgene_df["dof_true"]
        )
    ]

    # round up slope, t_stat, and pval_nom to 25 decimal places
    sgene_df["slope"] = sgene_df["slope"].round(5)
    sgene_df["t_stat"] = sgene_df["t_stat"].round(5)
    sgene_df["pval_nom"] = sgene_df["pval_nom"].round(5)

    pr_sgene_df = sgene_df[sgene_df["phenotype_id"].str.contains("PR")]
    up_sgene_df = sgene_df[sgene_df["phenotype_id"].str.contains("UP")]

    torus_cols = ["SNP", "gene", "beta", "t-stat", "p-value"]
    out_pr_sgene = pr_sgene_df[
        ["best_genotype_id", "phenotype_id", "slope", "t_stat", "pval_nom"]
    ]
    out_up_sgene = up_sgene_df[
        ["best_genotype_id", "phenotype_id", "slope", "t_stat", "pval_nom"]
    ]

    out_pr_sgene.columns = torus_cols
    out_up_sgene.columns = torus_cols
    with gzip.open(f"{outprefix}_{tissue}_p-sqtl.txt.gz", "wt") as f:
        out_pr_sgene.to_csv(f, sep="\t", index=False, header=True)
    with gzip.open(f"{outprefix}_{tissue}_u-sqtl.txt.gz", "wt") as f:
        out_up_sgene.to_csv(f, sep="\t", index=False, header=True)

    print(f"{tissue} done!")
