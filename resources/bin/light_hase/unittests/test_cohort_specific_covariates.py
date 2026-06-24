#!/usr/bin/env python3

import glob
import os
import unittest

import numpy as np
import pandas as pd

from unittests.utils.hase_executor import HaseExecutor


class TestCohortSpecificCovariates(unittest.TestCase):

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
        self.covariate_indices_path = os.path.join(
            self.resources_dir, "covariatesubsetting", "covariate_indices.txt")

        self.datasets = ["all_0", "sliced_0"]
        for study_name in self.datasets:

            self.dataset = os.path.join("covariatesubsetting")
            self.data_dir = os.path.join(self.resources_dir, "%s" % self.dataset)

            self.test_dir = os.path.join(self.project_root, "unittests", self.dataset, study_name, "out")
            if not os.path.exists(self.test_dir):
                os.makedirs(self.test_dir)

            hase_executor = HaseExecutor(
                self.test_dir, self.project_root, [study_name, ])
            hase_executor.set_covariate_indices_path(
                self.covariate_indices_path)
            hase_executor.meta_analyse(self.data_dir)
            self.hase_executor_map[study_name] = hase_executor

    def tearDown(self):
        pass

    def test_covariateSubsetting(self):

        # Load expected results
        results_dict = dict()

        # Now get the results for both executors. With the covariate pre-sliced
        # and with the covariate sliced within the mate-analysis
        for study_name, hase_executor in self.hase_executor_map.items():

            # Load actual results
            file_names = glob.glob(os.path.join(hase_executor.results_directory, "*.pkl"))
            results_dict[study_name] = pd.concat([pd.read_pickle(f) for f in file_names])

        merged = pd.merge(
            results_dict["all_0"], results_dict["sliced_0"],
            how="inner", on=["variant", "phenotype"],
            suffixes=('_tested', '_expected'))

        merged["beta_isclose"] = np.isclose(
            merged["beta_tested"], merged["beta_expected"])

        self.assertTrue(np.all(merged["beta_isclose"]))

        self.assertTrue(np.allclose(merged.standard_error_tested.values,
                                    merged.standard_error_expected.values))