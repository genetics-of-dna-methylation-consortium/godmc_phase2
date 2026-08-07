#!/usr/bin/env python3
import glob
import os
import shutil
import tempfile
import unittest
import distutils.core
import distutils.dir_util

from scipy import stats
from sklearn import linear_model

import hase
from tools import mapper

import numpy as np
import pandas as pd


def get_hase_results(results_directory):
    results_list = list()
    for file in os.listdir(results_directory):
        if file.endswith(".npy"):
            dict_with_nd_arrays = np.load(os.path.join(results_directory, file),
                                          allow_pickle=True).item()
            results_list.append(pd.DataFrame.from_dict(dict_with_nd_arrays))
    return pd.concat(results_list)


def get_sklearn_regression_difference(resources_directory, hase_results):
    genotypes = pd.read_csv(
        os.path.join(resources_directory, "exampledataset", "genotypes_csv", "dosage.csv"),
        header=0, index_col="ID")
    covariates = pd.read_csv(
        os.path.join(resources_directory, "exampledataset", "covariates", "example_study.csv"),
        sep="\t", header=0, index_col="id")
    phenotypes = pd.read_csv(
        os.path.join(resources_directory, "exampledataset", "phenotype", "example_study.csv"),
        sep="\t", header=0, index_col="id")

    for index, row in hase_results.iterrows():
        variant_id = row["index"]
        phenotype_id = row["phenotype"]
        data_frames = [genotypes.loc[variant_id], covariates]
        df_merged = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True), data_frames)
        sd, t_stat = fit_model(df_merged.values, phenotypes[phenotype_id][df_merged.index])
        print(row)
        print(sd, t_stat)
        # print(df_merged, phenotypes[phenotype_id][df_merged.index])


class HaseInteractionTestCase(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory
        #self.test_dir = tempfile.mkdtemp()
        self.test_dir = "/Users/cawarmerdam/Documents/hdaa/hase_repository/unittests/out"
        #self.test_dir = "/var/folders/ch/xw_1fk6j2v95btsrlchx63h40000gn/T/tmpvDZtAJ"
        self.resources_dir = os.path.join(os.path.dirname(__file__), "resources")
        self.project_root = os.path.dirname(os.path.dirname(__file__))

    def tearDown(self):
        pass
        #shutil.rmtree(self.test_dir)

    def test_get_sklearn_regression_results(self):
        get_sklearn_regression_difference(self.resources_dir,
                                          hase_results=get_hase_results(os.path.join(
                                           self.resources_dir, "exampledataset", "hase_results")))

    def test_meta_analysis(self):
        # Reference for mapper has to be downloaded
        assert os.path.isfile(
            os.path.join(self.project_root, "data", "1000Gp1v3.ref.gz"))
        assert os.path.isfile(
            os.path.join(self.project_root, "data", "1000Gp1v3.ref_info.h5"))

        # Start by running mapper?
        mapper_directory = os.path.join(self.test_dir, "mapper", "")
        hdf5_genotype_directory = os.path.join(self.resources_dir, "exampledataset", "genotypes_hdf5")
        # mapper.main(["-g", hdf5_genotype_directory,
        #              "-study_name", "dosage",
        #              "-o", mapper_directory])

        # Encode the genotype and phenotype files
        encoding_output = os.path.join(self.test_dir, "encoded")
        phenotype_matrix = os.path.join(self.resources_dir, "exampledataset", "phenotype")
        # hase.main(["-g", hdf5_genotype_directory,
        #            "-study_name", "dosage",
        #            "-o", encoding_output,
        #            "-mapper", mapper_directory,
        #            "-ph", phenotype_matrix,
        #            "-mode", "encoding"])

        # Calculate partial derivatives
        partial_derivatives = os.path.join(self.test_dir, "pd")
        # hase.main(["-g", hdf5_genotype_directory,
        #            "-study_name", "dosage",
        #            "-o", partial_derivatives,
        #            "-mapper", mapper_directory,
        #            "-ph", phenotype_matrix,
        #            "-cov", os.path.join(self.resources_dir, "exampledataset", "covariates"),
        #            "-mode", "single-meta"])

        # Create directory structure required for meta-analysis
        # Define directory names.
        meta = os.path.join(self.test_dir, "meta")
        genotype_meta = os.path.join(meta, "genotype_encoded")
        individuals_meta = os.path.join(genotype_meta, "individuals")
        actual_genotype_meta = os.path.join(genotype_meta, "genotype")
        phenotype_meta = os.path.join(meta, "phenotype")
        partial_derivatives_meta = os.path.join(meta, "pd_shared")

        # Copy the probes to genotype_meta. This calls os.makedirs() always,
        # and throws an exception if it already exists.
        # First we remove the genotype_meta if it does already exist
        if os.path.exists(meta):
            shutil.rmtree(meta)
        shutil.copytree(os.path.join(hdf5_genotype_directory, "probes"), os.path.join(genotype_meta, "probes"))

        # Make the required directories
        # os.makedirs(actual_genotype_meta)
        os.makedirs(individuals_meta)
        os.makedirs(partial_derivatives_meta)
        os.makedirs(phenotype_meta)

        # Copy the required data
        distutils.dir_util.copy_tree(os.path.join(encoding_output, "encode_individuals"), individuals_meta)
        distutils.dir_util.copy_tree(os.path.join(encoding_output, "encode_genotype"), actual_genotype_meta)
        distutils.dir_util.copy_tree(os.path.join(encoding_output, "encode_phenotype"), phenotype_meta)
        for file in glob.glob(os.path.join(partial_derivatives, "*.npy")):
            shutil.copy(file, partial_derivatives_meta)

        # Perform meta-analysis
        results_directory = os.path.join(self.test_dir, "results")
        hase.main(["-g", genotype_meta, "-study_name", "dosage",
                   "-ph", phenotype_meta, "-derivatives", partial_derivatives_meta,
                   "-mapper", mapper_directory, "-o", results_directory,
                   "-mode", "meta-classic", "-allow_missingness"])

        # (Generate the output file)

        # Read in the actual results
        hase_results = get_hase_results(results_directory)
        # Read in the expected results
        difference = get_sklearn_regression_difference(self.resources_dir, hase_results)


class HaseInteractionActualTestCase(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory
        #self.test_dir = tempfile.mkdtemp()
        self.test_dir = "/var/folders/ch/xw_1fk6j2v95btsrlchx63h40000gn/T/tmpvDZtAJ"
        self.resources_dir = os.path.join(os.path.dirname(__file__), "resources")
        self.project_root = os.path.dirname(os.path.dirname(__file__))

    def tearDown(self):
        pass
        #shutil.rmtree(self.test_dir)

    def test_get_sklearn_regression_results(self):
        get_sklearn_regression_difference(self.resources_dir,
                                          hase_results=get_hase_results(os.path.join(
                                           self.resources_dir, "exampledataset", "hase_results")))

    def test_meta_analysis(self):
        # Reference for mapper has to be downloaded
        assert os.path.isfile(
            os.path.join(self.project_root, "data", "1000Gp1v3.ref.gz"))
        assert os.path.isfile(
            os.path.join(self.project_root, "data", "1000Gp1v3.ref_info.h5"))

        # Start by running mapper?
        mapper_directory = os.path.join(self.test_dir, "mapper", "")
        hdf5_genotype_directory = os.path.join(self.resources_dir, "exampledataset", "genotypes_hdf5")
        # mapper.main(["-g", hdf5_genotype_directory,
        #              "-study_name", "dosage",
        #              "-o", mapper_directory])

        # Encode the genotype and phenotype files
        encoding_output = os.path.join(self.test_dir, "encoded")
        phenotype_matrix = os.path.join(self.resources_dir, "exampledataset", "phenotype")
        cov_interaction_matrix = os.path.join(self.resources_dir, "exampledataset", "covariates")
        # hase.main(["-g", hdf5_genotype_directory,
        #            "-study_name", "dosage",
        #            "-o", encoding_output,
        #            "-mapper", mapper_directory,
        #            "-ph", phenotype_matrix,
        #            "-interaction", cov_interaction_matrix,
        #            "-mode", "encoding"])

        # Calculate partial derivatives
        partial_derivatives = os.path.join(self.test_dir, "pd")
        # hase.main(["-g", hdf5_genotype_directory,
        #            "-study_name", "dosage",
        #            "-o", partial_derivatives,
        #            "-mapper", mapper_directory,
        #            "-ph", phenotype_matrix,
        #            "-cov", cov_interaction_matrix,
        #            "-interaction", cov_interaction_matrix,
        #            "-mode", "single-meta"])

        # Create directory structure required for meta-analysis
        # Define directory names.
        meta = os.path.join(self.test_dir, "meta")
        genotype_meta = os.path.join(meta, "genotype_encoded")
        individuals_meta = os.path.join(genotype_meta, "individuals")
        actual_genotype_meta = os.path.join(genotype_meta, "genotype")
        phenotype_meta = os.path.join(meta, "phenotype")
        interaction_meta = os.path.join(meta, "interaction")
        partial_derivatives_meta = os.path.join(meta, "pd_shared")

        # Copy the probes to genotype_meta. This calls os.makedirs() always,
        # and throws an exception if it already exists.
        # First we remove the genotype_meta if it does already exist
        if os.path.exists(meta):
            shutil.rmtree(meta)
        shutil.copytree(os.path.join(hdf5_genotype_directory, "probes"), os.path.join(genotype_meta, "probes"))

        # Make the required directories
        # os.makedirs(actual_genotype_meta)
        os.makedirs(individuals_meta)
        os.makedirs(partial_derivatives_meta)
        os.makedirs(phenotype_meta)
        os.makedirs(interaction_meta)

        # Copy the required data
        distutils.dir_util.copy_tree(os.path.join(encoding_output, "encode_individuals"), individuals_meta)
        distutils.dir_util.copy_tree(os.path.join(encoding_output, "encode_genotype"), actual_genotype_meta)
        distutils.dir_util.copy_tree(os.path.join(encoding_output, "encode_phenotype"), phenotype_meta)
        distutils.dir_util.copy_tree(os.path.join(encoding_output, "encode_interaction"), interaction_meta)
        for file in glob.glob(os.path.join(partial_derivatives, "*.npy")):
            shutil.copy(file, partial_derivatives_meta)

        # Perform meta-analysis
        results_directory = os.path.join(self.test_dir, "results")
        hase.main(["-g", genotype_meta,
                   "-study_name", "dosage",
                   "-ph", phenotype_meta,
                   "-interaction_encoded",
                   os.path.join(interaction_meta, "cov_0"),
                   os.path.join(interaction_meta, "cov_1"),
                   os.path.join(interaction_meta, "cov_2"),
                   "-derivatives", partial_derivatives_meta,
                   "-mapper", mapper_directory,
                   "-o", results_directory,
                   "-mode", "meta-stage"])

        # (Generate the output file)

        # Read in the actual results
        hase_results = get_hase_results(results_directory)
        # Read in the expected results
        difference = get_sklearn_regression_difference(self.resources_dir, hase_results)


def fit_model(X2, y):
    lm = linear_model.LinearRegression()
    lm.fit(X2, y)

    params = np.append(lm.intercept_, lm.coef_)
    predictions = lm.predict(X2)

    newX = pd.DataFrame({"Constant": np.ones(len(X2))}).join(pd.DataFrame(X2))
    MSE = (sum((y - predictions) ** 2)) / (len(newX) - len(newX.columns))

    # Note if you don't want to use a DataFrame replace the two lines above with
    # newX = np.append(np.ones((len(X),1)), X, axis=1)
    # MSE = (sum((y-predictions)**2))/(len(newX)-len(newX[0]))

    var_b = MSE * (np.linalg.inv(np.dot(newX.T, newX)).diagonal())
    sd_b = np.sqrt(var_b)
    ts_b = params / sd_b

    return [sd_b[-1], ts_b[-1]]


if __name__ == '__main__':
    unittest.main()
