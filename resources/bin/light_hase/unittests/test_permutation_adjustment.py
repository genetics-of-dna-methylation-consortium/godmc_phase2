from __future__ import print_function

import argparse
import os
import sys

#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#import seaborn as sns
from scipy.stats import genpareto

#sns.set_style(style="ticks")


def anderson_darling_gpd(sorted_values, genpareto_dist, genpareto_shape, table):
    N = len(sorted_values)
    # Get the cumulative density function, log transformed
    logcdf = genpareto_dist.logcdf(sorted_values)
    # Get the survival function, log transformed
    logsf = genpareto_dist.logsf(sorted_values)
    i = np.arange(1, N + 1)
    anderson_darling = (
            -N - np.sum((2*i - 1.0) / N * (logcdf + logsf[::-1]), axis=0))
    significance_values_character = np.array(table.columns[1:table.shape[1]])
    critical_values = np.array([np.interp(genpareto_shape, table["k"], table[col])
                                for col in significance_values_character])
    significance_values = significance_values_character.astype(float)
    return anderson_darling, critical_values, significance_values


def plot_genpareto(standard_scores, genpareto_fit, n_exceedances, label):
    abs_standard_scores_ordered = np.sort(np.abs(standard_scores))[::-1]
    # Calculate the threshold value
    exceedances_threshold = (
        abs_standard_scores_ordered[n_exceedances]
        + abs_standard_scores_ordered[n_exceedances+1]) / 2
    # Calculate the z argument for the generalized pareto distribution
    exceedence_values_sorted_ascending = np.sort(
        abs_standard_scores_ordered[(abs_standard_scores_ordered > exceedances_threshold)])
    z_test_statistics_for_fit = (
            exceedence_values_sorted_ascending
            - exceedances_threshold)
    # Calculate empirical CDF
    empirical_cdf = (
            np.arange(1, len(exceedence_values_sorted_ascending)+1)
            / float(len(exceedence_values_sorted_ascending)))
    # X values
    x_values = np.arange(z_test_statistics_for_fit.min(),
                         z_test_statistics_for_fit.max(), 0.1)
    theoretical_cdf = genpareto_fit.cdf(x_values)
    sns.lineplot(x=exceedence_values_sorted_ascending, y=empirical_cdf, label="empirical CDF")
    threshold = x_values + exceedances_threshold
    sns.lineplot(x=threshold, y=theoretical_cdf, label="theoretical CDF")
    sns.rugplot(a=exceedence_values_sorted_ascending)
    plt.text(exceedances_threshold, 1, label, horizontalalignment='left', size='medium', color='black',
             weight='semibold')
    plt.legend()
    plt.show()

def plot_pdf_over_maf(genpareto_distributions, statistics, maf_bins):
    # Every genpareto distribution starts at a different threshold
    x_values = np.arange(0,
                         np.concatenate(statistics.values()).max(), 0.1)
    dataframe_dict = dict()
    for maf_bin, (genpareto_fit, n_exceedances) in genpareto_distributions.items():
        # First calculate the exceedances threshold.
        abs_standard_scores_ordered = np.sort(np.abs(statistics[maf_bin]))[::-1]
        print(maf_bin)
        print(abs_standard_scores_ordered[n_exceedances - 1])
        print(abs_standard_scores_ordered[n_exceedances])
        # Calculate the threshold value
        exceedances_threshold = (
            abs_standard_scores_ordered[n_exceedances - 1]
            + abs_standard_scores_ordered[n_exceedances]) / 2
        density_function = genpareto_fit.pdf(x_values - exceedances_threshold)
        dataframe_dict[maf_bin] = pd.DataFrame({
            "x": x_values,
            "y": density_function})
    concatenated = pd.concat(
        dataframe_dict.values(),
        keys=dataframe_dict.keys()).rename_axis(['maf', 'i']).reset_index()
    concatenated["maf"] = pd.Categorical(concatenated["maf"],
                                         categories=maf_bins, ordered=True)
    ax = sns.lineplot(x='x', y='y', hue='maf', data=concatenated,
                      palette="viridis", hue_order=maf_bins)
    ax.set(ylim=(0, 0.1))
    plt.show()

class PermutationPvalueEstimator:
    def __init__(self, test_stat, anderson_darling_table):
        self.perm_stat = "perm_test_stat"
        self.gpd_fit = None
        self.n_exceedances = None
        self.anderson_darling_table = anderson_darling_table
        self.test_stat = test_stat
    def fit_permutations(self, permutation_results):
        self.gpd_fit, self.n_exceedances = (
            self.fit_distribution(permutation_results[self.test_stat],
                                  self.anderson_darling_table))
    def fit_distribution(self, test_statistics, anderson_darling_table,
                         n_exceedances_init=250):
        # Compute the permutation p-value using the
        # P_gpd method
        abs_standard_scores_ordered = np.sort(np.abs(test_statistics))[::-1]
        return self._fit_distribution(
            abs_standard_scores_ordered,
            n_exceedances_init,
            anderson_darling_table)
    def _fit_distribution(self, abs_standard_scores_ordered,
                          n_exceedances, anderson_darling_table):
        # Calculate the threshold value
        exceedances_threshold = (
            abs_standard_scores_ordered[n_exceedances - 1]
            + abs_standard_scores_ordered[n_exceedances]) / 2
        # Calculate the z argument for the generalized pareto distribution
        exceedence_values_sorted_ascending = np.sort(
            abs_standard_scores_ordered[(abs_standard_scores_ordered > exceedances_threshold)])
        z_test_statistics_for_fit = (
                exceedence_values_sorted_ascending
                - exceedances_threshold)
        # Fit the generalized pareto distribution
        fitted_parameters = genpareto.fit(z_test_statistics_for_fit)
        # Generate the probability density function
        genpareto_dist = genpareto(c=fitted_parameters[0],
                                   loc=fitted_parameters[1],
                                   scale=fitted_parameters[2])
        ad_statistic, crit_values, pvalues = anderson_darling_gpd(
            z_test_statistics_for_fit, genpareto_dist,
            fitted_parameters[0], anderson_darling_table)
        goodness_of_fit_pval = (pvalues[ad_statistic >= crit_values]).min(initial=1)
        print(n_exceedances, goodness_of_fit_pval)
        if goodness_of_fit_pval <= 0.05 and n_exceedances >= 10:
            return self._fit_distribution(
                abs_standard_scores_ordered, n_exceedances - 10,
                anderson_darling_table)
        elif goodness_of_fit_pval > 0.05:
            return genpareto_dist, n_exceedances
        else:
            return None, 0
    def estimate_pvalues(self,
                         permutation_results,
                         empirical_results):
        # Get the permuted dataframe for this maf bin
        permutation_data = (
            permutation_results.loc[:, np.array([self.test_stat])]
            .rename(columns={self.test_stat: self.perm_stat})
            .sort_values(by=[self.perm_stat]))
        empirical_data = empirical_results.sort_values(by=[self.test_stat])
        # Exceedances threshold
        exceedances_threshold = (
                (permutation_data[self.perm_stat].iloc[-self.n_exceedances]
                 + permutation_data[self.perm_stat].iloc[-(self.n_exceedances + 1)])
                / 2)
        # Get the number of permuted values
        n_permuted = float(permutation_data.shape[0])
        # Compute M
        permutation_data['M'] = (
            permutation_data[self.perm_stat].rank(ascending=False))
        pvalue_dataframe_bin = pd.merge_asof(
            empirical_data,
            permutation_data,
            left_on=self.test_stat, right_on=self.perm_stat,
            direction="forward"
        )
        pvalue_dataframe_bin.loc[pvalue_dataframe_bin['M'].isna(), 'M'] = 0
        # Compute P'ecdf for all empirical results
        pvalue_dataframe_bin['P_ecdf'] = (
                pvalue_dataframe_bin['M'] / n_permuted)
        pvalue_dataframe_bin['P_ecdf_pseudocount'] = (
                (pvalue_dataframe_bin['M'] + 1) / n_permuted)
        # For those empirical results that have an M < 10
        pvalue_dataframe_bin['P_gpd'] = -1
        # Select rows for which we can use the GPD
        gpd_rows = (pvalue_dataframe_bin['M'] < 10)
        # Get the generalized pareto distribution pdf for the rows
        # that are applicable for the GPD
        gpd_cdf = self.gpd_fit.cdf(
            pvalue_dataframe_bin.loc[gpd_rows, self.test_stat]
            - exceedances_threshold)
        # calculate the accurate permuted p-value
        pvalue_dataframe_bin.loc[gpd_rows, 'P_gpd'] = (
                self.n_exceedances / float(n_permuted) * (1 - gpd_cdf))
        pvalue_dataframe_bin['P_accurate'] = (
            pvalue_dataframe_bin['P_ecdf']
        )
        pvalue_dataframe_bin.loc[gpd_rows, 'P_accurate'] = (
            pvalue_dataframe_bin.loc[gpd_rows, 'P_gpd']
        )
        print(pvalue_dataframe_bin)
        print(pvalue_dataframe_bin.columns)
        return pvalue_dataframe_bin

def estimate_pvalues_over_maf_bins(
        empirical_results, permutation_results, maf_bins, pval_estimator):
    pval_estimators = dict()
    accurate_pvalue_estimations_list = list()
    for maf_bin in maf_bins:
        # Get variants of interest
        if np.sum(permutation_results["MAF_BIN"] == maf_bin) < 1:
            continue
        permutation_results_bin = permutation_results.loc[
            permutation_results["MAF_BIN"] == maf_bin]
        empirical_results_bin = empirical_results.loc[
            empirical_results["MAF_BIN"] == maf_bin]
        pval_estimator.fit_permutations(
            permutation_results_bin)
        accurate_pvalue_estimations_list.append(pval_estimator.estimate_pvalues(
            permutation_results_bin, empirical_results_bin))
        pval_estimators[maf_bin] = pval_estimator
    accurate_pvalue_estimations = pd.concat(accurate_pvalue_estimations_list)
    return(accurate_pvalue_estimations)

def plot_accurate_pvalues(accurate_pvalue_estimations, output_prefix):
    # ax = sns.countplot(accurate_pvalue_estimations["MAF_BIN"])
    # ax.tick_params(axis='x', rotation=90)
    # plt.show()
    #
    # # ax = sns.countplot(permutation_dataframe["MAF_BIN"])
    # # ax.tick_params(axis='x', rotation=90)
    # # plt.show()
    #
    accurate_pvalue_estimations["P_accurate_log10"] = -np.log10(
        accurate_pvalue_estimations["P_accurate"])
    accurate_pvalue_estimations["P_ecdf_psdoc_log10"] = -np.log10(
        accurate_pvalue_estimations["P_ecdf_pseudocount"])
    accurate_pvalue_estimations["P_gpd_log10"] = -np.log10(
        accurate_pvalue_estimations["P_gpd"])
    accurate_pvalue_estimations["P_emp_log10"] = -np.log10(
        accurate_pvalue_estimations["Pemp"])
    accurate_pvalue_estimations["MAF_BIN_I"] = \
        accurate_pvalue_estimations["MAF_BIN"].cat.codes
    # In theory, the smaller the MAF, the further the tail should reach.
    grid = sns.FacetGrid(accurate_pvalue_estimations,
                         col="MAF_BIN_I", col_wrap=4)
    grid.map(sns.scatterplot, "P_emp_log10", "P_accurate_log10",
             marker=".", alpha=0.5, linewidth=0)
    grid.fig.tight_layout(w_pad=1)
    grid.set(
        xlabel="empirical p-values (-log10 transformed)",
        ylabel="accurate perm p-values (-log10 transformed)")
    plt.savefig("{}.empirical_pval_panels.png".format(output_prefix))
    comb = sns.scatterplot(
        data=accurate_pvalue_estimations,
        x="P_emp_log10", y="P_accurate_log10",
        hue="MAF_BIN_I", palette="viridis",
        marker=".", alpha=0.5, linewidth=0)
    color_mapping = plt.cm.ScalarMappable(
        cmap="viridis", norm=plt.Normalize(vmin=0, vmax=19))
    color_mapping.set_array(accurate_pvalue_estimations["MAF_BIN_I"])
    comb.set(
        xlabel="empirical p-values (-log10 transformed)",
        ylabel="accurate perm p-values (-log10 transformed)")
    comb.get_legend().remove()
    comb.figure.colorbar(color_mapping)
    plt.savefig("{}.empirical_pval_overlaid.png".format(output_prefix))
    grid = sns.FacetGrid(accurate_pvalue_estimations, col="MAF_BIN", col_wrap=5)
    grid.map(sns.scatterplot, "P_ecdf_psdoc_log10", "P_accurate_log10",
             marker=".", alpha=0.5, linewidth=0)
    grid.fig.tight_layout(w_pad=1)
    grid.set(
        xlabel="ecdf p-values (with psdoc; -log10 transformed)",
        ylabel="accurate perm p-values (-log10 transformed)")
    plt.savefig("{}.ecdf_comparison.png".format(output_prefix))


def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Define argument parser
    parser = argparse.ArgumentParser(description="Correct p-values")
    parser.add_argument("--AD-test-table", dest="anderson_darling_path",
                        help="Anderson-Darling table for determining "
                             "goodness-of-fit p-values of generalized "
                             "pareto distributions")
    parser.add_argument("--emp",
                        help="File with empirical association statistics")
    parser.add_argument("--perm",
                        help="File with permutation association statistics")
    parser.add_argument("--maf", type=str,
                        help="File with minor allele frequencies")
    parser.add_argument("--breakpoints",
                        help="File with breakpoints for maf bins")
    parser.add_argument("--variants",
                        help="File with uncorrelated variants per maf bin")
    parser.add_argument("--output-prefix",
                        help="output prefix to write files to")

    # Parse args
    args = parser.parse_args(argv[1:])

    print(args)

    # Load the anderson-darling table
    anderson_darling_table = (
        pd.read_csv(args.anderson_darling_path).iloc[1::2, :])

    # Load unrelated variants (labeled by maf bin)
    uncorrelated_variants = pd.read_csv(
        args.variants, sep=" ", header=None, names=["SNP", "MAF_BIN_IDX"])

    # Breakpoints
    breakpoints = np.insert(
        pd.read_csv(args.breakpoints, header=None)[0].values, 0,0)

    # Get all MAF bin labels as R would also produce
    maf_bins = ["{0:.3g}-{1:.3g}".format(
        breakpoints[i], breakpoints[i+1])
        for i in range(len(breakpoints) - 1)]

    test_stat = "test_stat"

    maf_table = pd.read_csv(args.maf, sep="\t")
    maf_table["MAF"] = np.minimum(maf_table["AF"], 1 - maf_table["AF"])
    maf_table["MAF_BIN"] = pd.cut(
        maf_table["MAF"], breakpoints,
        right=True, include_lowest=True,
        labels=maf_bins)
    maf_table["SNP"] = maf_table["ID"]
    maf_table = maf_table[["SNP", "MAF", "MAF_BIN"]]
    maf_table_pruned = maf_table.merge(
        uncorrelated_variants[["SNP"]], how="inner", on="SNP")

    # Check if output path is available
    output_file_path = "{}.pval_estimations.csv.gz".format(args.output_prefix)
    if os.path.exists(output_file_path):
        raise argparse.ArgumentError(
            args.output_prefix,
            "There is already output present associated with the supplied prefix.")

    # Load data
    print("Loading Empirical results...")
    empirical_results = pd.read_csv(args.emp, sep="\t")
    print("Empirical results loaded!")

    print("Loading permutation results...")
    permutation_results = pd.read_csv(args.perm, sep="\t")
    print("Permutation results loaded!")

    print("Processing empirical results")
    empirical_results = process_assocation_results(
        maf_table, empirical_results, test_stat)
    print("Processing permutation results")
    permutation_results = process_assocation_results(
        maf_table_pruned, permutation_results, test_stat)

    print(permutation_results.columns)

    pval_estimator = PermutationPvalueEstimator(
        test_stat, anderson_darling_table)

    # Start iterating over genes

    write_header = True
    print("Estimating accurate p-values for:")
    phenotypes_to_loop_through = empirical_results["phenotype"].unique()

    for i, phenotype in enumerate(phenotypes_to_loop_through):
        print(
            "   {}, {}/{}\r"
            .format(phenotype, i, len(phenotypes_to_loop_through)))

        # Estimate pvalues for this phenotype
        accurate_pvalue_estimations = estimate_pvalues_over_maf_bins(
            empirical_results.loc[empirical_results["phenotype"] == phenotype],
            permutation_results.loc[permutation_results["phenotype"] == phenotype],
            maf_bins, pval_estimator)

        # plot_accurate_pvalues(
        #     accurate_pvalue_estimations,
        #     "{}.debug_{}".format(args.output_prefix, phenotype))

        accurate_pvalue_estimations.to_csv(
            "{}.pval_estimations.csv.gz".format(args.output_prefix),
            sep='\t', index=False,
            mode="a", header=write_header)

        write_header = False

    return 0


def process_assocation_results(maf_table, permutation_results, test_stat):
    permutation_results[test_stat] = (
            permutation_results["beta"] / permutation_results["se"]).abs().pow(3)
    permutation_results[test_stat] = (
            permutation_results["beta"] / permutation_results["se"]).abs().pow(3)
    print("Merging minor allele frequency table...")
    permutation_results = (permutation_results.merge(
        maf_table, how="inner", left_on="SNP", right_on="SNP")
    [["phenotype", "SNP", "P", "N", "MAF_BIN", "MAF", test_stat]])
    print("Merging done!")
    return permutation_results


if __name__ == "__main__":
    sys.exit(main())
