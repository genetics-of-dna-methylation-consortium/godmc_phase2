#!/usr/bin/env python3

import glob
import os
import shutil
import tempfile
import unittest
import numpy as np
import h5py

import pandas as pd

from hase import load_mapper
from hdgwas import hdregression, meta_classic
from hdgwas.data import Reader, MetaPhenotype, MetaParData, Hdf5Data
from hdgwas.tools import merge_genotype, get_intersecting_individual_indices
from tools import mapper
from unittests.test_classic_meta_analysis import load_expected_variant_indices
from unittests.test_hdregression import get_single_study_data
from unittests.utils.hase_executor import HaseExecutor


class TestInvalidEncodedSampleSlicing(unittest.TestCase):

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
        #self.test_dir = tempfile.mkdtemp()
        self.resources_dir = os.path.join(os.path.dirname(__file__), "resources")
        self.project_root = os.path.dirname(os.path.dirname(__file__))

        self.hase_executor_map = dict()

        self.study_name = "test_0"
        self.sample_to_exclude = "S700"

        self.dataset = os.path.join("complexdataset")
        self.data_dir = os.path.join(self.resources_dir, "%s" % self.dataset)

        self.base_test_dir = os.path.join(self.project_root, "unittests", "invalidsampleslicing", "out", self.dataset, self.study_name)

        # Get results without slicing out the sample
        test_dir_0 = os.path.join(self.base_test_dir, "actual")

        if not os.path.exists(test_dir_0):
            os.makedirs(test_dir_0)

        self.phenotype_matrix_alt = self.mutate_phenotype_sample()

        hase_executor_0 = HaseExecutor(
            test_dir_0, self.project_root, [self.study_name, ])
        hase_executor_0.meta_analyse(self.data_dir)
        self.hase_executor_map["default"] = hase_executor_0

        # Get results while slicing out the sample after encoding
        test_dir_1 = os.path.join(self.base_test_dir, "test")

        hase_executor_1 = HaseExecutor(
            test_dir_1, self.project_root, [self.study_name, ])
        hase_executor_0.copy_per_cohort_mapper_files(hase_executor_1)
        hase_executor_0.copy_per_cohort_analysis(hase_executor_1)

        self.remove_sample_from_encoded_data(hase_executor_1)

        hase_executor_1.meta_analyse(self.data_dir)
        self.hase_executor_map["test"] = hase_executor_1

        # Get results while slicing out the sample after encoding
        test_dir_2 = os.path.join(self.base_test_dir, "reference")

        hase_executor_2 = HaseExecutor(
            test_dir_2, self.project_root, [self.study_name, ])
        hase_executor_0.copy_per_cohort_mapper_files(hase_executor_2)

        self.per_cohort_analysis_with_sliced_sample(hase_executor_2)

        hase_executor_2.meta_analyse(self.data_dir)
        self.hase_executor_map["reference"] = hase_executor_2

        # Get results while slicing out the sample after encoding
        test_dir_3 = os.path.join(self.base_test_dir, "correct")

        hase_executor_3 = HaseExecutor(
            test_dir_3, self.project_root, [self.study_name, ])
        hase_executor_0.copy_per_cohort_mapper_files(hase_executor_3)

        if not hase_executor_3.per_cohort_results_exist(self.study_name):
            hase_executor_3.per_cohort_analysis(self.data_dir, self.study_name,
                                                phenotype_matrix=self.phenotype_matrix_alt)

        hase_executor_3.meta_analyse(self.data_dir)
        self.hase_executor_map["correct"] = hase_executor_3

    def mutate_phenotype_sample(self):
        phenotype_matrix = os.path.join(self.data_dir, self.study_name, "phenotypes")
        phenotype_matrix_alt = os.path.join(self.project_root, "unittests", "invalidsampleslicing", "phenotypes")
        if not os.path.exists(phenotype_matrix_alt):
            os.makedirs(phenotype_matrix_alt)
        phenotypes = pd.read_csv(os.path.join(phenotype_matrix, self.study_name + ".txt"), sep="\t")
        phenotypes.loc[self.sample_to_exclude == phenotypes.id, "id"] = "dummy"
        phenotypes.to_csv(os.path.join(phenotype_matrix_alt, self.study_name + ".txt"), sep="\t", index=False)
        return phenotype_matrix_alt

    def per_cohort_analysis_with_sliced_sample(self, hase_executor_2):
        # Start by running mapper?
        hdf5_genotype_directory = os.path.join(self.data_dir, self.study_name, "genotypes_hdf5")
        hase_executor_2.per_cohort_mapper_step(hdf5_genotype_directory,
                                               self.study_name)

        # Encode the genotype and phenotype files
        encoding_output = os.path.join(hase_executor_2.test_dir, "encoded")
        phenotype_matrix = os.path.join(self.data_dir, self.study_name, "phenotypes")
        cov_interaction_matrix = os.path.join(self.data_dir, self.study_name, "covariates")

        hase_executor_2.per_cohort_encoding_step(
            encoding_output, hdf5_genotype_directory,
            self.phenotype_matrix_alt, self.study_name)

        # Calculate partial derivatives
        partial_derivatives = os.path.join(hase_executor_2.test_dir, "pd")

        hase_executor_2.per_cohort_partial_derivatives(
            cov_interaction_matrix, hdf5_genotype_directory,
            partial_derivatives, phenotype_matrix, self.study_name)

        hase_executor_2.finalize_per_cohort(
            encoding_output, hdf5_genotype_directory,
            partial_derivatives, self.study_name)

    def tearDown(self):
        pass

    def test_InvalidEncodedSampleSlicing(self):

        # Load expected results
        results_dict = dict()

        # Now get the results for both executors. With the covariate pre-sliced
        # and with the covariate sliced within the mate-analysis
        for study_name, hase_executor in self.hase_executor_map.items():

            # Load actual results
            file_names = glob.glob(os.path.join(hase_executor.results_directory, "*.pkl"))
            results_dict[study_name] = pd.concat([pd.read_pickle(f) for f in file_names])

        concatenated = pd.concat(results_dict, ignore_index=False)
        concatenated = concatenated.droplevel(1)

        concatenated.to_csv(
            os.path.join(self.base_test_dir, "concatenated.csv"),
            sep="\t", index_label="run")

    def remove_sample_from_encoded_data(self, hase_executor):
        encoded_genotype_folder = hase_executor.get_meta_genotype_folder(self.study_name)
        h5_geno_files = os.listdir(os.path.join(encoded_genotype_folder, "genotype"))

        individuals_df = pd.read_hdf(
            os.path.join(encoded_genotype_folder, "individuals", self.study_name + '.h5'),
            'individuals')

        index_to_remove = individuals_df.individual == self.sample_to_exclude
        if np.any(index_to_remove):
            individuals_df.loc[index_to_remove] = 'dummy'

            individuals_df.to_hdf(
                os.path.join(encoded_genotype_folder, "individuals", self.study_name + '.h5'),
                key='individuals', format='table', min_itemsize=25, complib='zlib', complevel=9)

        # h5_gen_file = h5py.File(os.path.join(encoded_genotype_folder, "genotype", h5_geno_files[0]))
        #
        # atom = tables.Float64Atom()
        # self.genotype = self.h5_gen_file.create_carray(self.h5_gen_file.root, 'genotype', atom,
        #                                                (h5_gen_file.shape),
        #                                                title='Genotype', filters=self.pytable_filters)
        # self.genotype[:] = h5_gen_file
        # self.h5_gen_file.close()

