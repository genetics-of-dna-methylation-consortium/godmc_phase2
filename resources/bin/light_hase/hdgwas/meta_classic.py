from __future__ import print_function

import os
import gc

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.feather as ft

from hdgwas.hdregression import B4, \
    get_a_inverse_extended, hase_supporting_interactions, HASE
from hdgwas.tools import Timer, HaseAnalyser, select_identifiers, merge_genotype

class HaseException(Exception):
    pass


class ClassicMetaAnalyser:
    PYARROW_SCHEMA_META = pa.schema(
        [("phenotype", pa.string()),
         ("variant_index", pa.int64()),
         ("beta", pa.float64()),
         ("standard_error", pa.float64()),
         ("i_squared", pa.float64()),
         ("sample_size", pa.float64())])

    def __init__(self, meta_phen, meta_pard, sample_intersection, sample_indices, study_names, out,
                 selected_covariates=None, maf_threshold=0.0, t_statistic_threshold=None,
                 variants_full_log=None,
                 pheno_full_log=None, variant_phenotype_pairs_full_log=None,
                 t_threshold_full_log=None):

        self.meta_analysis_passed_t_threshold = None
        self.t_threshold_full_log = t_threshold_full_log
        self.variant_phenotype_pairs_full_log = variant_phenotype_pairs_full_log
        self.pheno_full_log = (
            pheno_full_log if pheno_full_log is not None else np.array([]))
        self.variants_full_log = (
            variants_full_log if variants_full_log is not None else np.array([]))
        self.t_statistic_threshold = t_statistic_threshold
        self.meta_pard = meta_pard
        self.meta_phen = meta_phen
        self.out = out
        self.cohort_list = self.construct_cohort_list(
            meta_pard, sample_intersection, sample_indices,
            study_names, selected_covariates=selected_covariates,
            maf_threshold=maf_threshold)
        self.output_type = "feather"
        self.results = None
        self.results_index = 0
        self._debug = True

    def construct_cohort_list(self, meta_pard, sample_intersection, sample_indices,
                              study_names, selected_covariates,
                              maf_threshold=0.0):
        # Build a list of studies / cohorts
        cohort_list = list()
        # Loop through studies / cohort names
        for study_index, study_name in enumerate(study_names):
            # Determine sample size
            sample_size = len(meta_pard.pd[study_name].folder._data.metadata['id'])

            # Now get the row indices matching with this studies samples,
            # so we can extract genotype and phenotype samples for this
            # study only.
            study_row_indices = select_identifiers(
                study_index, sample_indices, sample_intersection)
            
            # Get if covariate indices for this study
            # If the selected covariates file is non, we will set covariates_this_study to false,
            # This will make the covariate selector (CohortAnalyser.get_covariate_indices()) default to using
            # all covariates for which PDs were calculated
            covariates_this_study = False
            if selected_covariates is not None:
                # If covariate indices are not goven for this study
                # we can skip these. We make sure a keyerror is not thrown
                # by using None as a default. This will throw an error later on
                covariates_this_study = selected_covariates.get(
                    study_name, None)

            # Initialize the cohort analysis
            cohort = CohortAnalyser(
                study_index, study_name, study_row_indices, sample_size,
                selected_covariates=covariates_this_study, maf_threshold=maf_threshold)
            cohort.analyser.out = self.out

            # Add a dictionary as an index matching phenotypes to
            # indices in the partial derivatives C and b_cov
            cohort.set_partial_derivatives_phenotype_index(meta_pard)
            cohort_list.append(cohort)
        return cohort_list

    def analyse_genotype_chunk(self, genotype, variant_names, variant_indices,
                               chunk=None, node=None, chromosomes=None):

        for phenotype, phenotype_names, phenotype_indices in self.meta_phen:

            for cohort in self.cohort_list:
                # Now we need to do the analysis with the selected
                # phenotypes and genotypes
                cohort.set_phenotype_names(phenotype_names)
                cohort.set_phenotype_indices(phenotype_indices[cohort.study_index])

            self._analyse_genotype_chunk(
                variant_indices, variant_names,
                genotype, phenotype, phenotype_names,
                chunk, node, chromosomes)

    def _analyse_genotype_chunk(self, variant_indices, variant_names, genotype, phenotype, phenotype_names,
                                chunk=None, node=None, chromosomes=None):

        # Define matrices to do incremental meta-analysis
        weights_total = np.zeros((len(variant_names), len(phenotype_names)))
        weighted_effect_sizes_total = np.zeros(weights_total.shape)
        weighted_effect_sizes_squared_total = np.zeros(weights_total.shape)
        number_of_cohorts = np.zeros(weights_total.shape)

        # For every cohort, perform analysis, and add to the incremental meta-analysis
        for cohort in self.cohort_list:
            # Run the analysis for this cohort
            cohort.analyse(variant_indices, variant_names,
                           genotype, phenotype, self.meta_pard)

            self.save_results_per_cohort(cohort, variant_names, phenotype_names)

            # Perform incremental meta-analysis
            self.incremental_meta_analysis(
                cohort, number_of_cohorts,
                weighted_effect_sizes_squared_total, weighted_effect_sizes_total, weights_total)

            # Clear the per-cohort results
            cohort.analyser.t_stat = None
            cohort.analyser.standard_error = None

        # Finalize meta-analysis
        self.finalize_meta_analysis(
            number_of_cohorts, phenotype_names, variant_names,
            weighted_effect_sizes_squared_total, weighted_effect_sizes_total, weights_total)

        # Now do the classic meta-analysis
        self.save_results(chunk=chunk, node=node)
        self.results = None
        self.results_index += 1

    def finalize_meta_analysis(self, number_of_cohorts, phenotype_names, variant_names,
                               weighted_effect_sizes_squared_total, weighted_effect_sizes_total, weights_total):

        # The total weighted effect size should be divided by the
        # total weight of each of the studies.
        meta_analysed_effect_sizes = weighted_effect_sizes_total / weights_total
        meta_analysed_standard_error = np.sqrt(1 / weights_total)

        # Get the Q-statistic as indicated by Borenstein (p. 109, equation 16.3)
        q_statistic = weighted_effect_sizes_squared_total - (weighted_effect_sizes_total ** 2 / weights_total)

        # Get the correct degrees of freedom
        degrees_of_freedom = number_of_cohorts - 1

        # Calculate the I2 statistic
        with np.errstate(divide='ignore', invalid='ignore'):
            i_squared = np.maximum(0, (q_statistic - degrees_of_freedom) / q_statistic) * 100
            i_squared[np.isnan(i_squared)] = 101

        # Get the meta-analysed effect sizes
        meta_analysis_sample_sizes = np.nansum(
            self.get_sample_sizes(), axis=0)

        # Now create the cartesian product of all variant, phenotype names
        test_combinations = np.vstack(
            (np.repeat(variant_names, phenotype_names.shape[0]),
             np.tile(phenotype_names, variant_names.shape[0]))).T

        # Set the results
        self.set_results(meta_analysed_effect_sizes.flatten(),
                         meta_analysed_standard_error.flatten(),
                         meta_analysis_sample_sizes.flatten(),
                         i_squared.flatten(),
                         test_combinations)

    def incremental_meta_analysis(self, cohort, number_of_cohorts,
                                  weighted_effect_sizes_squared_total, weighted_effect_sizes_total, weights_total):

        # We rely on the methods below to return a 2d matrix wherein the
        # 0nd dimension represents the variants, and the
        # 1st dimension represents the phenotypes.
        effect_sizes = cohort.analyser.get_betas()
        inverse_variance_weights = 1 / (cohort.analyser.standard_error ** 2)

        # Update the totals
        weights_total += np.nan_to_num(inverse_variance_weights)
        weighted_effect_sizes_total += np.nan_to_num(effect_sizes * inverse_variance_weights)
        weighted_effect_sizes_squared_total += np.nan_to_num(
            inverse_variance_weights * effect_sizes ** 2)

        # Missing values are expected to be NaN
        # We should check which phenotype, variant combinations are NaN in all
        # cohorts. Then we can remove these from the variant and phenotype names
        number_of_cohorts += ~np.isnan(inverse_variance_weights)

        return

    def prepare_cohorts(self, variant_indices, variant_names):
        # First, prepare the genotype data for this chunk, per cohort
        for cohort in self.cohort_list:
            print("Preparing cohort named '{}'".format(cohort.study_name))
            cohort.set_variant_names(variant_names)
            # Get the partial derivatives for each study
            with Timer() as t_pd:
                # Extract the partial derivatives from the data sources
                a_test, b_cov, C, a_cov = self.meta_pard.get_single_study(
                    study_name=cohort.study_name,
                    study_index=cohort.study_index,
                    variant_indices=variant_indices)
                covariate_names = self.meta_pard.get_covariate_names(cohort.study_name)

            print("Time to get PD {}s".format(t_pd.secs))

            # Add the partial derivatives for this cohort
            cohort.set_partial_derivatives(a_test, b_cov, C, a_cov, covariate_names)
            # Add the maf filter for this cohort
            # cohort.set_observed_maf(self.meta_pard.minor_allele_frequencies_study(
            #     SNPs_index=variant_indices,
            #     study_name=cohort.study_name,
            #     study_index=cohort.study_index))
            cohort.set_variant_indices(variant_indices[cohort.study_index])
            # cohort.set_variant_filter(cohort.maf_filter())
            # Apply filters on partial derivatives and calculate processed
            # data structures.
            cohort.finalize_partial_derivatives()

    def meta_analyse(self, variant_names, phenotype_names):
        """
        Meta analysis for all associations.

        For each association, we apply an inverse variance meta analysis.

        For combining the results between cohorts, we complete the matrices
        with all missing values for ease of use.
        :param chromosomes:
        :param variant_names:
        :param phenotype_names:

        """
        # We rely on the methods below to return a 3d matrix wherein the
        # 0th dimension represents the studies / cohorts, the
        # 1st dimension represents the variants, and the
        # 2nd dimension represents the phenotypes.
        effect_sizes = self.get_effect_sizes()
        standard_error = self.get_standard_error()

        # Missing values are expected to be NaN
        # We should check which phenotype, variant combinations are NaN in all
        # cohorts. Then we can remove these from the variant and phenotype names
        not_all_nan = ~np.all(np.isnan(effect_sizes), axis=0)
        effect_sizes = effect_sizes[:, not_all_nan]
        standard_error = standard_error[:, not_all_nan]

        # Now create the cartesian product of all variant, phenotype names
        test_combinations = np.vstack(
            (np.repeat(variant_names, phenotype_names.shape[0]),
             np.tile(phenotype_names, variant_names.shape[0]))).T

        # From these combinations we can remove all that are nan in all cohorts
        test_combinations = test_combinations[not_all_nan.flatten()]

        # The inverse variance weighted meta-analysis requires calculating
        # the variance of the effect estimate. We do this by calculating the
        # 2nd power of the standard error.
        variance = standard_error ** 2

        # Weight the effect sizes by dividing the effect size by
        # the variance of the effect estimate.
        # This is the same thing as multiplying the effect size by
        # (1 / variance).
        effect_size_divided_by_variance = effect_sizes / variance
        effect_size_divided_by_variance_total = np.nansum(
            effect_size_divided_by_variance, axis=0)

        # Calculate the total weight as the inverse of the variance of the effect
        # estimate.
        one_divided_by_variance = 1 / variance
        one_divided_by_variance_total = np.nansum(one_divided_by_variance, axis=0)

        # The total weighted effect size should be divided still by the
        # total weight of each of the studies.
        meta_analysed_effect_sizes = effect_size_divided_by_variance_total / one_divided_by_variance_total
        meta_analysed_standard_error = np.sqrt(1 / one_divided_by_variance_total)

        # We need to quantify heterogeneity. We do this by calculating I2 values.
        # The first step is calculating a Q-value:
        # The weighted sum of squared differences between individual study
        # effects and the pooled effects across studies.
        effect_size_deviations_from_mean = np.power(
            effect_sizes - meta_analysed_effect_sizes, 2)

        # Weight and sum all effect size deviations from the mean.
        effect_size_deviations = np.nansum(
            one_divided_by_variance * effect_size_deviations_from_mean,
            axis=0)

        # Expected effect size deviations from the mean
        degrees_of_freedom = (np.sum(~np.isnan(effect_size_deviations_from_mean), axis=0) - 1)
        with np.errstate(divide='ignore'):
            i_squared = (
                ((effect_size_deviations - degrees_of_freedom)
                 / effect_size_deviations
                 ) * 100)

            # If the i_squared happens to be NaN, set it to a hundred.
            # Lack of precision causing the effect_size_deviations to be slightly different from 0
            # make the I_squared also one hundred.
            i_squared[np.isnan(i_squared)] = 100

        meta_analysis_sample_sizes = np.nansum(
            self.get_sample_sizes(), axis=0)[not_all_nan]

        self.set_results(meta_analysed_effect_sizes,
                         meta_analysed_standard_error,
                         meta_analysis_sample_sizes,
                         i_squared,
                         test_combinations)

    def set_results(self, meta_analysed_effect_sizes, meta_analysed_standard_error,
                    meta_analysis_sample_sizes, i_squared,
                    test_combinations):

        results_dataframe = pd.DataFrame({
            "variant_index": test_combinations[:, 0],
            "phenotype": test_combinations[:, 1],
            "beta": meta_analysed_effect_sizes,
            "standard_error": meta_analysed_standard_error,
            "i_squared": i_squared,
            "sample_size": meta_analysis_sample_sizes})

        t_statistic_mask = np.full((results_dataframe.shape[0]), True)

        if self.t_statistic_threshold is not None and self.t_statistic_threshold != 0:
            t_statistic_mask = (
                    np.abs(results_dataframe['beta'] / results_dataframe['standard_error'])
                    > self.t_statistic_threshold)

        sufficient_sample_size_mask = results_dataframe.sample_size != 0
        selection = np.logical_and(t_statistic_mask, sufficient_sample_size_mask)

        self.results = results_dataframe[selection]

    def write_debug(self):
        for cohort in self.cohort_list:
            file_path = os.path.join(
                self.out, "failed_", cohort.study_name, "_a_singular_{}.txt")
            for variant_name, a_singular_matrix in cohort.a_singular.items():
                np.savetxt(
                    file_path.format(variant_name),
                    a_singular_matrix)

    def save_results(self, chunk=None, node=None):

        print('Saving results to {}'.format(self.out))

        # Possibly set a threshold on p-value or Z score or t-statistic
        results = self.results.loc[self.results.sample_size != 0]
        output_path = self.get_output_path(chunk, node)

        results.columns = [col.decode('utf-8') for col in results.columns]

        results['phenotype'] = results['phenotype'].str.decode('utf-8')

        for name, group in results.groupby("phenotype".decode('utf-8'), sort=False):
            output_path_group = os.path.join(self.out, "meta", "phenotype={}".format(name))
            try:
                os.makedirs(output_path_group)
            except OSError:
                pass
            filename = os.path.join(output_path_group, "file_{}.feather".format(self.results_index))
            ft.write_feather(
                group,
                filename)

    def save_results_per_cohort(
            self, cohort, variant_names, phenotype_names):

        print('Saving results to {}'.format(self.out))

        # We rely on the methods below to return a 2d matrix wherein the
        # 0nd dimension represents the variants, and the
        # 1st dimension represents the phenotypes.
        effect_sizes = cohort.analyser.get_betas()
        standard_errors = cohort.analyser.standard_error

        effects_to_report = np.full_like(effect_sizes, False, dtype=bool)

        # Set all variant, phenotype combinations to true that where passed by a custom file
        if self.variant_phenotype_pairs_full_log is not None:
            for index, (variant_index, phenotype) in self.variant_phenotype_pairs_full_log:
                effects_to_report[
                    np.where(variant_names == variant_index),
                    np.where(phenotype_names == phenotype)] = True

        # Set all variant, phenotype combinations to true passing a certain
        # threshold in the per-cohort analysis
        if self.t_threshold_full_log is not None and self.t_threshold_full_log != 0:
            t_statistic_mask_per_cohort = (
                    np.abs(effects_to_report / standard_errors)
                    > self.t_threshold_full_log)

            effects_to_report |= t_statistic_mask_per_cohort

        elif self.t_threshold_full_log is not None:
            effects_to_report[...] = True

        variant_selection = np.isin(variant_names, self.variants_full_log)
        phenotype_selection = np.isin(phenotype_names, self.pheno_full_log)

        # We rely on the methods below to return a 3d matrix wherein the
        # 0th dimension represents the studies / cohorts, the
        # 1st dimension represents the variants, and the
        # 2nd dimension represents the phenotypes.

        effects_to_report[np.dot(
            variant_selection[:, None],
            phenotype_selection[None, :])] = True

        per_cohort_threshold = 1 * 10 ** 7
        if np.any(effects_to_report):
            pass
        elif len(self.pheno_full_log) == 0:
            print("No phenotypes provided for per cohort results. "
                  + "Trying to write per cohort results for all phenotypes.")
            phenotype_selection = np.array([True] * len(phenotype_names))

            effects_to_report[np.dot(
                variant_selection[:, None],
                phenotype_selection[None, :])] = True

            if np.sum(effects_to_report) > per_cohort_threshold:
                Warning("Refusing to write over {} effects for each cohort.".format(per_cohort_threshold))
                return

        elif len(self.variants_full_log) == 0:
            variant_selection = np.array([True] * len(variant_names))

            effects_to_report[np.dot(
                variant_selection[:, None],
                phenotype_selection[None, :])] = True

            if np.sum(effects_to_report) > per_cohort_threshold:
                Warning("Refusing to write over {} effects for each cohort.".format(per_cohort_threshold))
                return

        print("Per cohort effects to report")
        print(effects_to_report)

        if not np.any(effects_to_report):
            return

        effect_sizes_to_report = (
            effect_sizes[effects_to_report])
        standard_error_to_report = (
            standard_errors[effects_to_report])

        indices = np.where(effects_to_report)

        results_per_cohort = pd.concat([
            pd.DataFrame(
                {name: statistic})
            for name, statistic in zip(["beta", "standard_error"], [effect_sizes_to_report, standard_error_to_report])],
            axis=1)

        results_per_cohort.index = pd.MultiIndex.from_arrays(
            [np.array([cohort.study_name] * len(indices[1])),
             variant_names[indices[0]],
             phenotype_names[indices[1]]],
            names=['cohort', 'variant_index', 'phenotype'])

        results_per_cohort = results_per_cohort.reset_index()

        results_per_cohort.columns = [col.decode('utf-8') for col in results_per_cohort.columns]

        results_per_cohort['phenotype'] = results_per_cohort['phenotype'].str.decode('utf-8')
        results_per_cohort['cohort'] = results_per_cohort['cohort'].str.decode('utf-8')

        for key, group in results_per_cohort.groupby(["phenotype".decode('utf-8'), "cohort".decode('utf-8')]):
            phenotype = key[0]
            cohort = key[1]

            output_path_group = os.path.join(
                self.out, "cohort",
                "cohort={}".format(cohort),
                "phenotype={}".format(phenotype))
            try:
                os.makedirs(output_path_group)
            except OSError:
                pass

            filename = os.path.join(output_path_group, "file_{}.feather".format(self.results_index))
            ft.write_feather(
                group,
                filename)

    def get_output_label(self, chunk, node, per_cohort=False):
        if chunk is None:
            return str(self.results_index) + ('_per_cohort.' if per_cohort else '_result.')
        else:
            return "_".join(['node', str(node), str(chunk[0]), str(chunk[1]), str(self.results_index),
                             ('per_cohort.' if per_cohort else 'result.')])

    def get_output_path(self, chunk, node, per_cohort=False):
        return os.path.join(self.out,
                            self.get_output_label(chunk, node, per_cohort)
                            + self.output_type)

    def get_effect_sizes(self):
        """
        Method that returns a matrix of effect sizes.
        The first dimension represents variants,
        the second dimension represents the studies / cohorts.

        Missing associations are NA
        :param association_indices: A 3d matrix with for each study an array of variant indices, indicating at what position
        variants are located.
        :return: matrix of effect sizes [assoc, studies]
        """
        betas = np.stack([cohort.analyser.get_betas() for cohort in self.cohort_list])
        return betas

    def get_standard_error(self):
        """
        Method that returns a matrix of standard errors of effect sizes.
        The first dimension represents associations,
        the second dimension represents the studies / cohorts.
        :param association_indices: A list with for each study an array of variant indices, indicating at what position
        variants are located
        """
        standard_errors = np.stack([cohort.analyser.standard_error for cohort in self.cohort_list])
        return standard_errors

    def get_sample_sizes(self):
        return np.stack([cohort.get_sample_sizes() for cohort in self.cohort_list])

    def get_cohort_names(self):
        return [cohort.study_name for cohort in self.cohort_list]

    def run(self, mapper, gen, allowed_missingness_rate=1,
            is_no_b4_present_in_partial_derivatives=False):

        for phenotype, phenotype_names, phenotype_indices in self.meta_phen:

            for cohort in self.cohort_list:
                # Now we need to do the analysis with the selected
                # phenotypes and genotypes
                cohort.set_phenotype_names(phenotype_names)
                cohort.set_phenotype_indices(phenotype_indices[cohort.study_index])

            variant_indices = None
            # Start looping over all genotype chunks
            # We use while true, since implementing an iterator is apparently
            # too much work...
            while True:
                # We start with an empty chunk...
                ch = None
                node = None

                # Then we check whether or not we are working in an environment
                # wherein analysis is split over multiple nodes.

                # This determines how the current chunk should be obtained.
                if mapper.cluster == 'n':
                    variant_indices, variant_names = mapper.get(
                        allowed_missingness_rate=allowed_missingness_rate, chromosomes=False)
                elif mapper.cluster == 'y':
                    ch = mapper.chunk_pop()
                    node = mapper.node[1]
                    if ch is None:
                        break
                    variant_indices, variant_names = mapper.get(
                        chunk_number=ch, allowed_missingness_rate=allowed_missingness_rate, chromosomes=False)

                # Now we can check if we have looped through all genotype chunks...
                if isinstance(variant_indices, type(None)):
                    break

                # And we can get the genotypes for this chunk...
                if is_no_b4_present_in_partial_derivatives:
                    with Timer() as t_g:
                        genotype = merge_genotype(gen, variant_indices, mapper)
                    print("Time to get G {}s".format(t_g.secs))

                # Now analyse the genotype chunk
                self._analyse_genotype_chunk(
                    variant_indices, variant_names, genotype,
                    phenotype, phenotype_names,
                    ch, node)


class CohortAnalyser:
    """
    Class that handles HASE analysis for a single cohort.

    Attempts to limit memory usage as much as possible. It is written
    with the following workflow in mind:

    Initialize an instance using basic metadata for the specific
    cohort.

    For each genotype chunk, we can do the following with this class

    1. Process and store partial derivative structures on a genotype level that
    require much time to process. By doing this, we limit having
    to recalculate these for every set of phenotypes.

    2. Do only store filters for variants and phenotypes. Do not store
    actual filtered data. This will accumulate to quite a bit of data
    when looping over phenotypes without clearing.

    3. Do the analysis, filtering data on the fly.


    """

    def __init__(self, study_index, study_name, row_index,
                 sample_size, selected_covariates=None,
                 maf_threshold=0.0,
                 a_test=None, b_cov=None, C=None, a_cov=None):
        self._encoded_covariates = None
        self._selected_covariates = selected_covariates
        self.partial_derivatives_phenotype_index = None
        self.unknown_dimension_size = 1
        self.sample_size = sample_size
        self._a_cov = a_cov
        self._a_test = a_test
        self.a_inv = None
        self.a_singular = dict()
        self.phenotype_names = None
        self.variant_indices = None
        self.phenotype_indices = None
        self.C = C
        self._b_cov = b_cov
        self.study_name = study_name
        self.study_index = study_index
        self.sample_indices = row_index
        self.analyser = HaseAnalyser()
        self.number_of_variable_terms = None
        self.number_of_constant_terms = None
        self._finalized = False
        self.maf_threshold = maf_threshold

    def set_partial_derivatives(self, a_test, b_cov, C, a_cov, covariates):
        """
        Set partial derivatives.
        :param covariates:
        """
        self._b_cov = b_cov
        self.C = C
        self._a_test = a_test
        self._a_cov = a_cov
        self._encoded_covariates = ["_".join([self.study_name, covariate]) for covariate in covariates]

    def set_partial_derivatives_phenotype_index(self, meta_pard):
        """
        Set partial derivatives phenotype index. :param index: the index should be a dictionary with the names of the
        phenotype as keys, and the accompanying index in the partial derivatives as values.
        """
        phenotype_names = np.array(meta_pard.phen_names)
        phen_ind_series = pd.Series(
            np.arange(0, len(phenotype_names)),
            index=phenotype_names)
        self.partial_derivatives_phenotype_index = phen_ind_series

    def finalize_partial_derivatives(self):
        a_test = self.get_a_test()
        a_cov = self.get_a_cov()
        variant_names = self.get_variant_names()

        # Use new A_inverse that supports variable number of non-constant A parts
        a_inv, a_invertible = get_a_inverse_extended(a_cov, a_test)

        # Update variant indices with all variants for which the matrix
        # cannot be inverted
        self.a_singular = dict(zip(
            variant_names[~a_invertible],
            a_inv[~a_invertible]))

        self.add_variant_filter(a_invertible)
        self.a_inv = a_inv[a_invertible]

        self._set_number_of_variable_terms()
        self._set_number_of_constant_terms()
        self._finalized = True

    def add_variant_filter(self, filter):
        indices_ = self.variant_indices[self.variant_indices != -1]
        indices_[~filter] = -1
        self.variant_indices[self.variant_indices != -1] = indices_

    def get_a_cov(self):
        return self._a_cov[np.ix_(self.get_covariate_indices(),
                                  self.get_covariate_indices())]

    def get_a_test(self):
        if self.variant_indices is None:
            raise HaseException("variant indices must be set when working with partial derivatives")

        # When slicing a_test, we need to account for the variable terms, which
        # we need to all maintain.

        # We do this by first calculating which indices these should correspond
        # with:
        indices_of_variable_terms = np.arange(
            self._a_test.shape[1] - self._a_test.shape[2],  # The second dimension (1) of a_test corresponds to the
            # number of total independent variables, (constant terms, including the intercept, and variable terms)
            # The third dimension (2) corresponds to the number of variable terms. This is most commonly 1. but when
            # interaction terms are included (not fully supported at the time of writing this) this can exceed 1.
            # Subtracting both numbers should result in just the number of constant terms (including intercept)
            self._a_test.shape[1])  # Stop at the total number of independent variables (exclusive).

        # Now we
        covariate_indices_with_variable_terms = np.append(
            self.get_covariate_indices(), indices_of_variable_terms)

        return self._a_test[np.ix_(self.variant_indices != -1, covariate_indices_with_variable_terms)]

    def _set_number_of_constant_terms(self):
        self.number_of_constant_terms = self.get_a_cov().shape[0]

    def _set_number_of_variable_terms(self):
        self.number_of_variable_terms = self.get_a_test().shape[2]

    def set_observed_maf(self, maf):
        self.analyser.maf = maf

    def maf_filter(self):
        threshold = self.get_maf_threshold()
        variant_filter = np.array([True] * self.analyser.maf.shape[0])
        if threshold != 0:
            variant_filter = (self.analyser.maf >= threshold) \
                             & (self.analyser.maf <= (1 - threshold)) \
                             & (self.analyser.maf != 0.5)
            # This probably handles the large a_test matrices as well.
        return variant_filter

    def set_variant_indices(self, variant_indices):
        self.variant_indices = variant_indices

    def set_variant_filter(self, variant_filter):
        self.variant_indices[~variant_filter] = -1

    def set_phenotype_names(self, phenotype_names):
        self.phenotype_names = phenotype_names

    def set_phenotype_indices(self, phenotype_indices):
        self.phenotype_indices = phenotype_indices

    def apply_variant_filters(self, genotype):
        variant_indices = self.variant_indices
        if variant_indices is None:
            variant_filter = [True] * genotype.shape[0]
        else:
            variant_filter = self.variant_indices != -1

        filtered_genotype = genotype[np.ix_(variant_filter, self.sample_indices['genotype'])]

        if genotype.shape[0] == 0:
            raise HaseException("No SNPs pass variant filters")

        return filtered_genotype

    def apply_phenotype_filters_on_phenotype_matrix(self, phenotype):
        phenotype_indices = self.phenotype_indices
        if phenotype_indices is None:
            phenotype_filter = [True] * phenotype.shape[1]
        else:
            phenotype_filter = self.phenotype_indices != -1

        filtered_phenotype = phenotype[np.ix_(self.sample_indices['phenotype'], phenotype_filter)]

        print("Selected phenotype shape {}".format(filtered_phenotype.shape))

        return filtered_phenotype

    def get_a_inverse(self):
        return self.a_inv

    def get_b_cov(self):
        """
        This method returns the covariable part of b.
        :return: matrix b_cov
        """
        phenotype_indexes_for_partial_derivatives = self.get_phenotype_slicer_for_partial_derivatives()
        return self._b_cov[np.ix_(self.get_covariate_indices(), phenotype_indexes_for_partial_derivatives)]

    def get_C(self):
        """
        This method returns the C matrix
        :return: matrix C as a 2d numpy array.
        """
        phenotype_indexes_for_partial_derivatives = self.get_phenotype_slicer_for_partial_derivatives()
        return self.C[phenotype_indexes_for_partial_derivatives]

    def get_phenotype_slicer_for_partial_derivatives(self):
        """
        Method that returns an array of indexes with which phenotypes can
        be selected from the partial derivatives.

        We first need to apply the phenotype filters to only select those
        phenotypes that we wish to use in this cohort / analysis.

        Then we use this filtering stuff and the partial derivatives phenotype
        index to get the indices that match the phenotypes.

        :return: Array of indices.
        """
        phenotype_names = self.get_sliced_phenotype_names()
        # Now we match the filter phenotype names to individual indexes.
        if len(phenotype_names) == 0:
            raise HaseException('There is no common ids in phenotype files and PD data!')
        else:
            print('There are {} common ids in phenotype files and PD data!'.format(
                len(phenotype_names)))
        return self.partial_derivatives_phenotype_index[phenotype_names].values

    def get_sliced_phenotype_names(self):
        # First make a local copy of the phenotype filter
        phenotype_indices = self.phenotype_indices
        # Then, check if the phenotype filter has been set, if this is not the case,
        # overwrite the filter with True. This selects every phenotype
        if phenotype_indices is None:
            phenotype_filter = [True] * self.phenotype_names.shape[0]
        else:
            phenotype_filter = self.phenotype_indices != -1
        # Now, using the filter, select all names that we want to select
        phenotype_names = self.phenotype_names[phenotype_filter]
        return phenotype_names

    def analyse(self, variant_indices, variant_names, genotype, phenotype, meta_pard):
        print("Preparing cohort named '{}'".format(self.study_name))
        self.set_variant_names(variant_names)
        # Get the partial derivatives for each study
        with Timer() as t_pd:
            # Extract the partial derivatives from the data sources
            a_test, b_cov, C, a_cov = meta_pard.get_single_study(
                study_name=self.study_name,
                study_index=self.study_index,
                variant_indices=variant_indices)
            covariate_names = meta_pard.get_covariate_names(self.study_name)

        print("Time to get PD {}s".format(t_pd.secs))

        # Add the partial derivatives for this cohort
        self.set_partial_derivatives(a_test, b_cov, C, a_cov, covariate_names)
        # Add the maf filter for this cohort
        self.set_observed_maf(meta_pard.minor_allele_frequencies_study(
            SNPs_index=variant_indices,
            study_name=self.study_name,
            study_index=self.study_index))
        self.set_variant_indices(variant_indices[self.study_index])
        self.set_variant_filter(self.maf_filter())
        # Apply filters on partial derivatives and calculate processed
        # data structures.
        self.finalize_partial_derivatives()

        # Check if there are phenotypes in this chunk
        # If not, the _analyse function will fail.
        if len(self.get_sliced_phenotype_names()) == 0:

            # Therefore, should the number of phenotype names be 0,
            # we set the t_stats and standard error to be missing for all
            # phenotypes and genotypes.
            print("No phenotypes in cohort {} for this chunk. Skipping...".format(self.study_name))

            self.analyser.t_stat = \
                self.complete_output_with_missingness(np.array([]))
            self.analyser.standard_error = \
                self.complete_output_with_missingness(np.array([]))
        else:

            # If we have some phenotypes, we can continue.
            self._analyse(genotype, phenotype)

        self._cleanup()

    def _analyse(self, genotype, phenotype):

        # Get partial derivatives
        C_test = self.get_C()
        b_cov_test = self.get_b_cov()
        a_inverse = self.get_a_inverse()

        # Determine the number of constant terms
        number_of_constant_terms = self.get_number_of_constant_terms()

        # Calculate the B4 matrix
        b4 = self.calculate_b4_matrix(genotype, phenotype)
        b_variable = b4[np.newaxis, ...]

        # Calculate the degrees of freedom
        degrees_of_freedom = (self.sample_size - a_inverse.shape[1])

        print("B4 shape is {}".format(b4.shape))

        # t_stat, standard_error = HASE(
        #     b_variable[0, ...], a_inverse, b_cov_test, C_test,
        #     number_of_constant_terms, degrees_of_freedom)

        t_stat, standard_error = hase_supporting_interactions(
            b_variable, a_inverse, b_cov_test, C_test,
            number_of_constant_terms, degrees_of_freedom)

        self.analyser.t_stat = self.complete_output_with_missingness(t_stat[:, 0, :])
        self.analyser.standard_error = self.complete_output_with_missingness(standard_error[:, 0, :])
        # Now we have to store the t-stat, and standard_error in a data structure from
        # which we can do the meta analysis.

    def calculate_b4_matrix(self, genotype, phenotype):
        b4 = B4(self.apply_phenotype_filters_on_phenotype_matrix(phenotype),
                self.apply_variant_filters(genotype))
        return b4

    def complete_output_with_missingness(self, output_without_missingness):
        output_with_missingness = np.full((self.variant_indices.shape[0], self.phenotype_indices.shape[0]), np.NaN)
        output_with_missingness[
            np.ix_(self.variant_indices != -1, self.phenotype_indices != -1)] = output_without_missingness
        return output_with_missingness

    def get_number_of_constant_terms(self):
        return self.number_of_constant_terms

    def get_maf_threshold(self):
        return self.maf_threshold

    def set_variant_names(self, variant_names):
        self.analyser.rsid = variant_names

    def get_sample_sizes(self):
        sample_sizes = np.zeros((self.variant_indices.shape[0], self.phenotype_indices.shape[0]))
        sample_sizes[np.ix_(self.variant_indices != -1, self.phenotype_indices != -1)] = self.sample_size
        return sample_sizes

    def get_variant_names(self):
        variant_names = self.analyser.rsid
        variant_indices = self.variant_indices
        if variant_indices is None:
            variant_filter = [True] * variant_names
        else:
            variant_filter = self.variant_indices != -1

        return variant_names[variant_filter]

    def get_phenotype_names(self, phenotype_names):
        phenotype_indices = self.phenotype_indices
        if phenotype_indices is None:
            phenotype_filter = [True] * phenotype_names.shape[1]
        else:
            phenotype_filter = self.phenotype_indices != -1

        return phenotype_names[phenotype_filter]

    def get_betas(self, include_missing_values=False):
        betas = self.analyser.get_betas()
        if not include_missing_values:
            return betas[np.ix_(self.variant_indices != -1, self.phenotype_indices != -1)]
        return betas

    def get_covariate_indices(self):
        """
        Method that returns the covariate indices with a 0 included. This way
        it will also maintain the intercept thingy.
        :return:
        """
        if self._selected_covariates is None:
            raise ValueError(
                'Covariate selection not present for cohort: {}'.format(
                    self.study_name
                ))
            #return np.arange(self._a_cov.shape[0])
        elif self._selected_covariates is False:
            self._selected_covariates = self._encoded_covariates

        covariate_indices = np.where(np.in1d(self._encoded_covariates, self._selected_covariates))[0]

        if len(covariate_indices) < len(self._selected_covariates):
            missing_covariates = self._selected_covariates[~np.isin(self._selected_covariates,
                                                                    self._encoded_covariates)]
            raise ValueError(
                'Cannot find all covariates for study: {}'
                'Cannot find: {}'.format(
                    self.study_name, ", ".join(missing_covariates)
                ))

        if np.isin(covariate_indices[0], 0):
            return covariate_indices

        return np.insert(covariate_indices, obj=0, values=0, axis=0)

    def _cleanup(self):
        self._b_cov = None
        self.C = None
        self._a_test = None
        self._a_cov = None
        self._encoded_covariates = None
        self.a_singular = dict()
        self.a_inv = None

        gc.collect()