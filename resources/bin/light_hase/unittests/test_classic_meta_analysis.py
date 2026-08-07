#!/usr/bin/env python3

"""
Test script for the classic meta analysis
"""
import glob
import os
import unittest

import numpy as np
# TODO: implement test for missingness in phenotype thingy
# TODO: implement test for missingness in genotype thingy (new getter)
import pandas as pd
import pyarrow.parquet as pq
import pyarrow.feather as ft

from hase import load_mapper
from hdgwas import hdregression, meta_classic
from hdgwas.data import Reader, MetaPhenotype, MetaParData
from hdgwas.tools import merge_genotype, get_intersecting_individual_indices
from unittests.utils.hase_executor import HaseExecutor


def load_expected_variant_indices(path):
    return pd.read_csv(path, sep="\t", header=0, index_col="cohort")


class TestClassicMetaAnalyser(unittest.TestCase):

    def setUp(self):
        """
        Method that sets up testing data for this test case.

        Initial steps:
        Use generated data. Load this as if a hase run is started

        For testing the classic meta analyser
        1. Get a mapper ready,
        2. Get a meta phen ready
        3. Get a meta pard ready

        Make sure to set the hase chunking for both variants and genotypes
        to 10 to 20 or so to make sure chunking works as expected

        """
        # self.test_dir = tempfile.mkdtemp()
        self.resources_dir = os.path.join(os.path.dirname(__file__), "resources")
        self.project_root = os.path.dirname(os.path.dirname(__file__))
        self.dataset = "largensmallm"
        self.data_dir = os.path.join(self.resources_dir, "%s" % self.dataset)

        self.test_dir = os.path.join(self.project_root, "unittests", self.dataset, "out")
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir)

        self.study_names = ["test_0", "test_1"]

        self.hase_executor = HaseExecutor(
            self.data_dir,
            self.test_dir, self.project_root,
            self.study_names, new_mapper_options=False)

        # self.hase_executor.node_i = 2
        # self.hase_executor.nodes_n = 10
        #
        # self.hase_executor = HaseExecutor(
        #     self.test_dir, self.project_root,
        #     self.study_names, "1000G-30x_ref",
        #     "/Users/cawarmerdam/Documents/hdaa/hase_repository/unittests/resources/alleleflips/hase_reference",
        #     False)

        self.hase_executor.meta_analyse()

    def tearDown(self):
        pass

    def test_mapperGetNonMissing(self):
        """
        """

        expected_genotypes = self.expected_genotypes()

        mapper = load_mapper(
            mapper_chunk_size=30, study_name=self.hase_executor.study_names,
            ref_name=self.hase_executor.ref_name,
            mapper_folder=self.hase_executor.get_meta_mapper_folder(),
            encoded=[0] * len(self.hase_executor.study_names),
            cluster=None, node=None,
            snp_id_inc=None, snp_id_exc=None)

        gen = []
        for i, study_name in enumerate(self.hase_executor.study_names):
            gen.append(Reader('genotype'))
            gen[i].start(os.path.join(self.data_dir, study_name, "genotypes_hdf5"),
                         hdf5=True, study_name=study_name, ID=False)

        keys_list = list()
        actual_genotypes_list = list()

        while True:
            variant_indices, keys = mapper.get(
                allowed_missingness_rate=0)

            if isinstance(variant_indices, type(None)):
                break

            keys_list.append(keys)

            actual_genotypes_list.append(merge_genotype(gen, variant_indices, mapper))

        actual_genotypes = np.concatenate(actual_genotypes_list)
        actual_keys = np.concatenate(keys_list)

        self.assertTrue(np.array_equal(
            actual_genotypes,
            expected_genotypes.loc[actual_keys, :].values))

    def test_mapperGetLimitedMissingness(self):
        """
        Compare expected variant_indices,
        as created when generating testing data,
        with the actual variant_indices as returned by the
        mapper.get(allow_missingness=T) method.

        Make sure the
        """
        # TODO: test that the number of variants is according to what is expected

        expected_genotypes = self.expected_genotypes()

        mapper = load_mapper(
            mapper_chunk_size=30, study_name=self.hase_executor.study_names,
            ref_name=self.hase_executor.ref_name,
            mapper_folder=self.hase_executor.get_meta_mapper_folder(),
            encoded=[0] * len(self.hase_executor.study_names),
            cluster=None, node=None,
            snp_id_inc=None, snp_id_exc=None)

        gen = []
        for i, study_name in enumerate(self.hase_executor.study_names):
            gen.append(Reader('genotype'))
            gen[i].start(os.path.join(self.data_dir, study_name, "genotypes_hdf5"),
                         hdf5=True, study_name=study_name, ID=False)

        keys_list = list()
        actual_genotypes_list = list()

        while True:
            variant_indices, keys = mapper.get(
                allowed_missingness_rate=0.5)

            if isinstance(variant_indices, type(None)):
                break

            keys_list.append(keys)

            actual_genotypes_list.append(merge_genotype(gen, variant_indices, mapper))

        actual_genotypes = np.concatenate(actual_genotypes_list)
        actual_keys = np.concatenate(keys_list)

        self.assertTrue(np.array_equal(
            actual_genotypes,
            expected_genotypes.loc[actual_keys, :].values))


    def test_mapperGetMissingness(self):
        """
        Compare expected variant_indices,
        as created when generating testing data,
        with the actual variant_indices as returned by the
        mapper.get(allow_missingness=T) method.

        Make sure the
        """

        expected_genotypes = self.expected_genotypes()

        mapper = load_mapper(
            mapper_chunk_size=30, study_name=self.hase_executor.study_names,
            ref_name=self.hase_executor.ref_name,
            mapper_folder=self.hase_executor.get_meta_mapper_folder(),
            encoded=[0] * len(self.hase_executor.study_names),
            cluster=None, node=None,
            snp_id_inc=None, snp_id_exc=None)

        gen = []
        for i, study_name in enumerate(self.hase_executor.study_names):
            gen.append(Reader('genotype'))
            gen[i].start(os.path.join(self.data_dir, study_name, "genotypes_hdf5"),
                         hdf5=True, study_name=study_name, ID=False)

        keys_list = list()
        actual_genotypes_list = list()

        while True:
            variant_indices, keys = mapper.get(
                allowed_missingness_rate=1)

            if isinstance(variant_indices, type(None)):
                break

            keys_list.append(keys)

            actual_genotypes_list.append(merge_genotype(gen, variant_indices, mapper))

        actual_genotypes = np.concatenate(actual_genotypes_list)
        actual_keys = np.concatenate(keys_list)

        self.assertTrue(np.array_equal(
            actual_genotypes,
            expected_genotypes.loc[actual_keys, :].values))


    def test_CohortAnalyser(self):
        study_name = "test_2"
        study_index = 2
        max_missingness_rate = 1

        genotypes_path = os.path.join(self.data_dir, study_name, "genotypes.csv")
        genotypes = pd.read_csv(genotypes_path, delimiter="\t", index_col='ID')

        phenotypes_path = os.path.join(self.data_dir, study_name, "phenotypes", study_name + ".txt")
        phenotypes = pd.read_csv(phenotypes_path, delimiter="\t", index_col='id')

        number_of_participants = float(genotypes.values.shape[1])
        freq_aa = np.count_nonzero(genotypes.values == 0., axis=1) / number_of_participants
        freq_ab = np.count_nonzero(genotypes.values == 1., axis=1) / number_of_participants
        freq_bb = np.count_nonzero(genotypes.values == 2., axis=1) / number_of_participants
        freq_a = freq_aa + freq_ab * 0.5
        freq_b = freq_bb + freq_ab * 0.5
        maf_filter = np.logical_and(freq_a >= 0.05, freq_b >= 0.05)

        expected_b4_data_frame = pd.DataFrame(
            hdregression.B4(phenotypes.values, genotypes.values[maf_filter, :]),
            columns=phenotypes.columns, index=genotypes.index[maf_filter])

        expected_results = pd.read_csv(
            os.path.join(self.data_dir, study_name, "results.txt"),
            delimiter="\t")

        expected_betas = expected_results.pivot(index='variant', columns='phenotype', values='Coefficients')

        mapper = load_mapper(
            mapper_chunk_size=200, study_name=self.hase_executor.study_names,
            ref_name=self.hase_executor.ref_name,
            mapper_folder=self.hase_executor.get_meta_mapper_folder(),
            encoded=[1] * len(self.hase_executor.study_names),
            cluster=None, node=None,
            snp_id_inc=None, snp_id_exc=None)

        partial_derivatives_folders = []

        for i, j in enumerate(self.hase_executor.run_func_all_studies(
                self.hase_executor.get_partial_derivatives_meta_folder)):
            partial_derivatives_folders.append(Reader('partial'))
            partial_derivatives_folders[i].start(j, study_name=self.hase_executor.study_names[i])
            partial_derivatives_folders[i].folder.load()

        meta_pard = MetaParData(partial_derivatives_folders, self.hase_executor.study_names,
                                protocol=None,
                                allowed_missingness_rate=max_missingness_rate)

        phen = []
        for i, j in enumerate(self.hase_executor.run_func_all_studies(
                self.hase_executor.get_meta_phenotype_folder)):
            phen.append(Reader('phenotype'))
            phen[i].start(j)

        meta_phen = MetaPhenotype(phen,
                                  include=None, exclude=None, allowed_missingness_rate=max_missingness_rate)
        meta_phen.chunk_size = 5

        gen = []
        for i, j in enumerate(self.hase_executor.run_func_all_studies(
                self.hase_executor.get_meta_genotype_folder)):
            gen.append(Reader('genotype'))
            gen[i].start(j, hdf5=True, study_name=self.hase_executor.study_names[i], ID=False)

        # Create a dictionary with the datasets for obtaining
        # the indices of shared identifiers
        datasets = {"phenotype": tuple(i.folder._data for i in phen),
                    "genotype": tuple(i.folder._data for i in gen),
                    "partial_derivatives": tuple(i.folder._data.metadata for i in partial_derivatives_folders)}

        # Get common ids
        row_index, intersecting_identifiers = get_intersecting_individual_indices(datasets)

        classic_meta_analyser = meta_classic.ClassicMetaAnalyser(meta_phen, meta_pard, intersecting_identifiers,
                                                                 row_index, self.hase_executor.study_names,
                                                                 self.hase_executor.results_directory,
                                                                 maf_threshold=0.05, t_statistic_threshold=0)

        variant_indices, variant_names = mapper.get(
            allowed_missingness_rate=max_missingness_rate)
        genotype = merge_genotype(gen, variant_indices, mapper)

        classic_meta_analyser.prepare_cohorts(
            variant_indices, variant_names)

        for phenotype, phenotype_names, phenotype_indices in classic_meta_analyser.meta_phen:

            cohort = classic_meta_analyser.cohort_list[study_index]
            cohort.set_phenotype_names(phenotype_names)
            cohort.set_phenotype_indices(phenotype_indices[cohort.study_index])
            actual_b4_data_frame = pd.DataFrame(
                cohort.calculate_b4_matrix(genotype, phenotype),
                columns=cohort.get_phenotype_names(phenotype_names),
                index=cohort.get_variant_names())

            required_phenotype_names = phenotype_names[
                np.isin(phenotype_names, expected_b4_data_frame.columns.values)]

            self.compare_actual_matrix_to_expected(
                actual_b4_data_frame, expected_b4_data_frame[required_phenotype_names])

            cohort.analyse(variant_indices, variant_names,
                           genotype, phenotype,
                           classic_meta_analyser.meta_pard)
            actual_betas = pd.DataFrame(
                cohort.get_betas(),
                columns=cohort.get_phenotype_names(phenotype_names),
                index=cohort.get_variant_names())

            self.compare_actual_matrix_to_expected(
                actual_betas, expected_betas[required_phenotype_names])

    def test_per_cohort_results(self):
        max_missingness_rate = 1

        mapper = load_mapper(
            mapper_chunk_size=200, study_name=self.hase_executor.study_names,
            ref_name=self.hase_executor.ref_name,
            mapper_folder=self.hase_executor.get_meta_mapper_folder(),
            encoded=[1] * len(self.hase_executor.study_names),
            cluster=None, node=None,
            snp_id_inc=None, snp_id_exc=None)

        gen = []
        for i, j in enumerate(self.hase_executor.run_func_all_studies(
                self.hase_executor.get_meta_genotype_folder)):
            gen.append(Reader('genotype'))
            gen[i].start(j, hdf5=True, study_name=self.hase_executor.study_names[i], ID=False)

        variants_of_interest = np.array(["rs79067516", "21-24257984"])
        phenotypes_of_interest = np.array(["pheno_13", "pheno_09"])

        classic_meta_analyser = self.construct_classic_meta_analyser(
            max_missingness_rate, gen, variants_of_interest, phenotypes_of_interest)

        variant_indices, variant_names, chromosomes = mapper.get(
            allowed_missingness_rate=max_missingness_rate, chromosomes=True)
        genotype = merge_genotype(gen, variant_indices, mapper)

        for phenotype, phenotype_names, phenotype_indices in classic_meta_analyser.meta_phen:

            for cohort in classic_meta_analyser.cohort_list:
                # Now we need to do the analysis with the selected
                # phenotypes and genotypes
                cohort.set_phenotype_names(phenotype_names)
                cohort.set_phenotype_indices(phenotype_indices[cohort.study_index])
                cohort.analyse(variant_indices, variant_names,
                               genotype, phenotype,
                               classic_meta_analyser.meta_pard)

            # Now do the classic meta-analysis
            classic_meta_analyser.save_results_per_cohort(
                variant_names, phenotype_names)

    def construct_classic_meta_analyser(self, max_missingness_rate, gen,
                                        variants_of_interest = None,
                                        phenotypes_of_interest = None):
        partial_derivatives_folders = []
        for i, j in enumerate(self.hase_executor.run_func_all_studies(
                self.hase_executor.get_partial_derivatives_meta_folder)):
            partial_derivatives_folders.append(Reader('partial'))
            partial_derivatives_folders[i].start(j, study_name=self.hase_executor.study_names[i])
            partial_derivatives_folders[i].folder.load()
        meta_pard = MetaParData(partial_derivatives_folders, self.hase_executor.study_names,
                                protocol=None,
                                allowed_missingness_rate=max_missingness_rate)
        phen = []
        for i, j in enumerate(self.hase_executor.run_func_all_studies(
                self.hase_executor.get_meta_phenotype_folder)):
            phen.append(Reader('phenotype'))
            phen[i].start(j)
        meta_phen = MetaPhenotype(phen,
                                  include=None, exclude=None, allowed_missingness_rate=max_missingness_rate)
        meta_phen.chunk_size = 5

        # Create a dictionary with the datasets for obtaining
        # the indices of shared identifiers
        datasets = {"phenotype": tuple(i.folder._data for i in phen),
                    "genotype": tuple(i.folder._data for i in gen),
                    "partial_derivatives": tuple(i.folder._data.metadata for i in partial_derivatives_folders)}
        # Get common ids
        row_index, intersecting_identifiers = get_intersecting_individual_indices(datasets)
        classic_meta_analyser = meta_classic.ClassicMetaAnalyser(meta_phen, meta_pard, intersecting_identifiers,
                                                                 row_index, self.hase_executor.study_names,
                                                                 self.hase_executor.results_directory,
                                                                 variants_full_log=variants_of_interest,
                                                                 pheno_full_log=phenotypes_of_interest,
                                                                 maf_threshold=0.05, t_statistic_threshold=0)
        return classic_meta_analyser

    def compare_actual_matrix_to_expected(self, actual_data_frame, expected_data_frame):
        self.assertTrue(np.all(np.isin(
            expected_data_frame.index.values,
            actual_data_frame.index.values)))
        self.assertTrue(np.all(np.isin(
            expected_data_frame.columns.values,
            actual_data_frame.columns.values)))
        self.assertTrue(np.allclose(
            expected_data_frame.loc[
                actual_data_frame.index,
                actual_data_frame.columns].values,
            actual_data_frame.values))

    def test_MetaPhenotypeGetMissingness(self):
        expected_phenotypes = self.expected_phenotypes()

        meta_phen = self.setup_meta_phen()

        phenotype, phenotype_names = meta_phen.get()

        self.assertTrue(np.array_equal(
            expected_phenotypes.loc[:, phenotype_names].values,
            phenotype))

    def setup_meta_phen(self):
        phen = []
        for i, study_name in enumerate(self.hase_executor.study_names):
            phen.append(Reader('phenotype'))
            phen[i].start(os.path.join(self.data_dir, study_name, "phenotypes"))
        meta_phen = MetaPhenotype(
            phen,
            include=None, exclude=None, allowed_missingness_rate=1)
        return meta_phen

    def expected_genotypes(self):
        expected_genotypes_list = list()
        for i, study_name in enumerate(self.hase_executor.study_names):
            expected_genotypes_path = os.path.join(self.data_dir, study_name, "genotypes.csv")
            expected_genotypes_list.append(pd.read_csv(
                expected_genotypes_path, delimiter="\t", index_col='ID'))
        expected_genotypes = pd.concat(expected_genotypes_list, axis=1, sort=True)
        expected_genotypes = expected_genotypes.fillna(0)

        return expected_genotypes

    def test_metaAnalysisResults(self):

        # Load expected results
        output_list = list()

        for study_name in self.hase_executor.study_names:
            output_list.append(pd.read_csv(
                os.path.join(self.data_dir, study_name, "results.txt"),
                delimiter="\t"))

        # Load actual results
        file_names_pkl = glob.glob(os.path.join(self.hase_executor.results_directory, "*result.pkl"))
        file_names_parquet = glob.glob(os.path.join(self.hase_executor.results_directory, "*result.parquet"))
        file_names_feather = glob.glob(
            os.path.join(self.hase_executor.results_directory, "meta", "phenotype=*", "*.feather"))

        if len(file_names_parquet) > 0:
            actual_results = pd.concat([pq.read_table(f).to_pandas() for f in file_names_parquet])
        elif len(file_names_pkl) > 0:
            actual_results = pd.concat([pd.read_pickle(f) for f in file_names_pkl])
        elif len(file_names_feather) > 0:
            actual_results = pd.concat([ft.read_feather(f) for f in file_names_feather])
        else:
            self.fail("No output present")
        expected_results = pd.read_csv(
            os.path.join(self.data_dir, "meta.txt"), delimiter='\t')
        print(actual_results)

        if "variant" in actual_results.columns:
            actual_results_proc = actual_results
        elif "variant_index" in actual_results.columns:
            variant_reference = pd.read_csv(self.hase_executor.ref_file, sep=" ", header=0, usecols=[0, 1, 2, 3, 4], names=["variant", "bp", "ref_allele", "eff_allele", "chromosome"])
            actual_results["variant_index"] = actual_results["variant_index"].astype(np.int64)
            actual_results_proc = pd.merge(actual_results, variant_reference.loc[actual_results.variant_index.unique(), :],
                                           how="inner", left_on="variant_index", right_index=True)

        merged = pd.merge(
            actual_results_proc, expected_results,
            how="outer", on=["variant", "phenotype"],
            suffixes=('_actual', '_expected'))

        merged["beta_isclose"] = np.isclose(
            merged["beta_actual"], merged["beta_expected"])

        self.assertTrue(np.allclose(
            merged.standard_error_actual.values,
            merged.standard_error_expected.values))

    def expected_phenotypes(self):
        expected_phenotype_list = list()
        for i, study_name in enumerate(self.hase_executor.study_names):
            expected_phenotypes_path = os.path.join(self.data_dir, study_name, "phenotypes", study_name + ".txt")
            expected_phenotype_list.append(pd.read_csv(
                expected_phenotypes_path, delimiter="\t", index_col='id'))
        expected_phenotypes = pd.concat(expected_phenotype_list, axis=0, sort=True)
        expected_phenotypes = expected_phenotypes.fillna(0)

        return expected_phenotypes


if __name__ == '__main__':
    unittest.main()
