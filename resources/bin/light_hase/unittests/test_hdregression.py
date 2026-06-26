#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from __future__ import with_statement
from __future__ import absolute_import
import sys

import unittest

import numpy as np
import pandas as pd
from sklearn import linear_model, preprocessing
from sklearn.preprocessing import PolynomialFeatures
from scipy import stats

from hdgwas.hdregression import A_covariates, A_tests, B_covariates, C_matrix, B4, A_inverse, HASE, \
    calculate_variant_dependent_a, get_a_inverse_extended, hase_supporting_interactions, expand_sample_size_matrix, \
    a_inverse_extended_allow_missingness


class EncoderCopy:
    def __init__(self, number_of_individuals):
        self.F = np.random.randint(
            1, 10, number_of_individuals * number_of_individuals).reshape(
            number_of_individuals, number_of_individuals)
        self.F_inv = np.linalg.inv(self.F)

    def encode_genotype_matrix(self, genotype_matrix):
        return np.dot(genotype_matrix, self.F)

    def encode_with_inverse(self, phenotype_matrix):
        return np.dot(self.F_inv, phenotype_matrix)


class TestIndividualShuffling(unittest.TestCase):
    number_of_individuals = 40
    number_of_phenotypes = 200

    def test_individual_shuffling(self):

        phenotype_matrix = np.random.randint(
            1, 10, self.number_of_individuals * self.number_of_phenotypes).reshape(
            self.number_of_individuals, self.number_of_phenotypes)
        phenotype_matrix_shuffled = np.copy(phenotype_matrix)
        np.random.shuffle(phenotype_matrix_shuffled)

        # Encode the data
        encoder = EncoderCopy(self.number_of_individuals)
        phenotype_matrix_encoded = encoder.encode_with_inverse(phenotype_matrix_shuffled)

        resolved = np.transpose(
            np.linalg.lstsq(np.transpose(phenotype_matrix), np.transpose(phenotype_matrix_encoded))[0])

        print(resolved)
        print(encoder.F_inv)


class TestInteractionSupportingHase(unittest.TestCase):
    number_of_individuals = 81
    number_of_variants = 8

    def test_without_interactions(self):
        genotype_matrix = get_genotype_matrix(self.number_of_individuals,
                                              self.number_of_variants)
        covariates_matrix = get_covariates_matrix(self.number_of_individuals, 1)  # One covariate
        covariates_matrix_extended = covariates_matrix.repeat(self.number_of_variants, axis=1)
        base_x = np.stack((covariates_matrix_extended, genotype_matrix.T), axis=1)
        matrix_with_determinants_x = get_x_no_interaction(base_x)

        phenotype_matrix = get_phenotype_matrix_no_interaction(matrix_with_determinants_x)[:, 0:2]

        # Encode the data
        encoder = EncoderCopy(self.number_of_individuals)
        genotype_matrix_encoded = encoder.encode_genotype_matrix(genotype_matrix)
        phenotype_matrix_encoded = encoder.encode_with_inverse(phenotype_matrix)

        # Calculate the A inverse the new way
        a_inverse_alternative = calculate_a_alternative(covariates_matrix, genotype_matrix)

        # Calculate the A inverse the old way
        a_inverse = calculate_a(covariates_matrix, genotype_matrix)

        b_cov = B_covariates(covariates_matrix, phenotype_matrix)

        C = C_matrix(phenotype_matrix)

        b4 = B4(phenotype_matrix_encoded, genotype_matrix_encoded)

        b_variable = b4[np.newaxis, ...]

        number_of_variable_terms = b_variable.shape[0]

        N_con = a_inverse_alternative.shape[1] - number_of_variable_terms

        DF = (self.number_of_individuals - a_inverse_alternative.shape[1])

        t_stat, SE = HASE(b4, a_inverse, b_cov, C, N_con, DF)

        t_stat2, SE2 = hase_supporting_interactions(b_variable, a_inverse_alternative, b_cov, C, N_con, DF)

        self.assertTrue(np.allclose(SE, SE2))
        self.assertTrue(np.allclose(t_stat, t_stat2))
        print("Standard error and t-statistics equal between old model and new model")

        model_table = fit_model(matrix_with_determinants_x[..., 0], phenotype_matrix[..., 0])

        self.assertTrue(np.isclose(model_table.loc[a_inverse_alternative.shape[1] - 1, u"Standard Errors"], SE[0, 0, 0]))
        self.assertTrue(np.isclose(model_table.loc[a_inverse_alternative.shape[1] - 1, u"t values"], t_stat[0, 0, 0]))
        print("Standard error and t-statistics equal between model and regular regression analysis")

        # TODO: check if this does test the hase_supporting_interactions function?

    def test_without_interaction_with_variant_subsetting(self):

        a_cov_1, variant_dependent_a_1, b_cov_1, C_1, b4_1, matrix_with_determinants_x_1, phenotype_matrix_1 = \
            get_single_study_data(self.number_of_individuals, self.number_of_variants)

        # Set 3rd variant to be missing in second study
        missing_variant = 3

        a_cov_2, variant_dependent_a_2, b_cov_2, C_2, b4_2, matrix_with_determinants_x_2, phenotype_matrix_2 = \
            get_single_study_data(self.number_of_individuals, self.number_of_variants)

        # Set data for missing variant to 0
        variant_dependent_a_2[missing_variant] = 0
        b4_2[missing_variant] = 0

        matrix_with_determinants_x = np.concatenate((matrix_with_determinants_x_1, matrix_with_determinants_x_2), axis=0)

        # print(np.isclose(matrix_with_determinants_x[100:200, ...],matrix_with_determinants_x_2))
        phenotype_matrix = np.concatenate((phenotype_matrix_1, phenotype_matrix_2), axis = 0)

        a_inverse_alternative_1 = a_inverse_extended_allow_missingness(a_cov_1, variant_dependent_a_1)

        b_variable_1 = b4_1[np.newaxis, ...]

        number_of_variable_terms = b_variable_1.shape[0]

        N_con = a_inverse_alternative_1.shape[1] - number_of_variable_terms

        DF = (self.number_of_individuals - a_inverse_alternative_1.shape[1])

        t_stat_1, SE_1 = hase_supporting_interactions(b_variable_1, a_inverse_alternative_1, b_cov_1, C_1, N_con, DF)

        a_cov = a_cov_1 + a_cov_2
        variant_dependent_a = variant_dependent_a_1 + variant_dependent_a_2
        b_cov = b_cov_1 + b_cov_2
        C = C_1 + C_2
        b4 = b4_1 + b4_2

        a_inverse_alternative, _ = get_a_inverse_extended(a_cov, variant_dependent_a)
        a_inverse_alternative[missing_variant,::], _ = get_a_inverse_extended(a_cov_1,
                                                                           variant_dependent_a_1)[missing_variant,::]

        b_variable = b4[np.newaxis, ...]

        number_of_variable_terms = b_variable.shape[0]

        N_con = a_inverse_alternative.shape[1] - number_of_variable_terms

        # DF = (self.number_of_individuals * 2 - a_inverse_alternative.shape[1])

        sampleCount = np.array([200 for i in range(b4.size)]).reshape(b4.shape)[:, np.newaxis, :]
        sampleCount[missing_variant,...] = 100

        print(sampleCount)

        DF = (sampleCount - a_inverse_alternative.shape[1])

        t_stat2, SE2 = hase_supporting_interactions(b_variable, a_inverse_alternative, b_cov_1, C_1, N_con, DF)

        model_table = fit_model(matrix_with_determinants_x[..., 0], phenotype_matrix[..., 0])
        #
        # self.assertTrue(np.isclose(
        #     model_table.loc[a_inverse_alternative.shape[1] - 1, u"Standard Errors"],
        #     SE2[0, 0, 0]))
        # self.assertTrue(np.isclose(
        #     model_table.loc[a_inverse_alternative.shape[1] - 1, u"t values"],
        #     t_stat2[0, 0, 0]))

        model_table_missing = fit_model(matrix_with_determinants_x_1[..., missing_variant], phenotype_matrix_1[..., 0])

        self.assertTrue(np.isclose(
            model_table_missing.loc[a_inverse_alternative.shape[1] - 1, u"Standard Errors"],
            SE_1[missing_variant, 0, 0]))

        self.assertTrue(np.isclose(
            model_table_missing.loc[a_inverse_alternative.shape[1] - 1, u"t values"],
            t_stat_1[missing_variant, 0, 0]))

        self.assertTrue(np.isclose(
            model_table_missing.loc[a_inverse_alternative.shape[1] - 1, u"Standard Errors"],
            SE2[missing_variant, 0, 0]))
        #
        # self.assertTrue(np.isclose(
        #     model_table_missing.loc[a_inverse_alternative.shape[1] - 1, u"t values"],
        #     t_stat2[missing_variant, 0, 0]))

        print("Standard error and t-statistics equal between model and regular regression analysis")

        # model_table = fit_model(matrix_with_determinants_x[..., 1], phenotype_matrix[..., 1])

        # self.assertTrue(np.isclose(model_table.loc[a_inverse_alternative.shape[1] - 1, u"Standard Errors"], SE2[0, 0, 0]))
        # self.assertTrue(np.isclose(model_table.loc[a_inverse_alternative.shape[1] - 1, u"t values"], t_stat2[0, 0, 0]))
        # print("Standard error and t-statistics equal between model and regular regression analysis")

    def test_with_interactions(self):
        # Get a genotype matrix
        genotype_matrix = get_genotype_matrix(self.number_of_individuals,
                                              self.number_of_variants)
        # Get the matrix of covariates.single
        covariates_matrix = get_covariates_matrix(self.number_of_individuals, 2)  # two covariates.single
        # Get a matrix of covariates.single that is extended to be the same shape as the
        # genotype matrix
        covariates_matrix_extended = np.tile(covariates_matrix[:, np.newaxis, :], (1, self.number_of_variants, 1))
        # Create the X matrix with base determinants genotype and the covariate
        base_X = np.concatenate((covariates_matrix_extended, genotype_matrix.T[:, :, np.newaxis]), axis=2).transpose(
            (0, 2, 1))
        # Get the X matrix with the
        X = get_x_with_interaction(base_X)
        # Get the phenotype matrix
        phenotype_matrix = get_phenotype_matrix_with_interaction(X)[:, 0:2]

        # Encode the stuff
        encoder = EncoderCopy(self.number_of_individuals)
        genotype_matrix_encoded = encoder.encode_genotype_matrix(genotype_matrix)
        phenotype_matrix_encoded = encoder.encode_with_inverse(phenotype_matrix)

        encoded_interaction_phenotype_matrices = list()
        # Calculate the phenotypes * interaction covariate
        for covariate_array in covariates_matrix.T:
            interaction_phenotype_matrix = np.einsum('ij,i->ij', phenotype_matrix, covariate_array)
            print(interaction_phenotype_matrix)
            encoded_interaction_phenotype_matrices.append(encoder.encode_with_inverse(
                interaction_phenotype_matrix))

        # Get the constant part of A
        a_cov = A_covariates(covariates_matrix)
        variant_dependent_A = calculate_variant_dependent_a(genotype_matrix,
                                                            covariates_matrix,
                                                            covariates_matrix)

        # Get the b part that corresponds to the interaction values
        b_interaction_list = list()
        for encoded_interaction_phenotype_matrix in encoded_interaction_phenotype_matrices:
            b_interaction_list.append(np.dot(genotype_matrix_encoded, encoded_interaction_phenotype_matrix))
        b_interaction = np.array(b_interaction_list)

        # Get the constant part of b
        b_cov = B_covariates(covariates_matrix, phenotype_matrix)

        C = C_matrix(phenotype_matrix)

        b4 = B4(phenotype_matrix_encoded, genotype_matrix_encoded)

        a_inv, _ = get_a_inverse_extended(a_cov, variant_dependent_A)

        b_variable = np.append(b_interaction, b4[np.newaxis, ...], axis=0)

        number_of_variable_terms = b_variable.shape[0]

        N_con = a_inv.shape[1] - number_of_variable_terms

        DF = (self.number_of_individuals - a_inv.shape[1])

        t_stat, SE = hase_supporting_interactions(b_variable, a_inv, b_cov, C, N_con, DF)

        print(X[..., 0])
        print(phenotype_matrix[..., 0])

        model_table = fit_model(X[..., 0], phenotype_matrix[..., 0])

        print(model_table)

        # Assert whether the HASE2 method returns equal results compared to the regular
        # regression analyses.
        self.assertTrue(np.isclose(model_table.loc[a_inv.shape[1] - 1, u"Standard Errors"], SE[0, 0, 0]))
        self.assertTrue(np.isclose(model_table.loc[a_inv.shape[1] - 1, u"t values"], t_stat[0, 0, 0]))
        print("Standard error and t-statistics equal between model and regular regression analysis")

#
# class TestMissingnessFunctions(unittest.TestCase):
#
#     def test_expand_sample_size_matrix(self):
#         # Define three cohorts, with 52, 23, and 81 samples respectively
#         sample_size = [52, 23, 81]
#
#         # We assume to have 8 variants available
#         # Every cohort has a different set of variants available (and missing)
#         variant_indices = [np.array([1, 3, 4, 6, -1, 9, 10, 12]),
#                            np.array([-1, -1, 3, 5, 6, 8, -1, 11]),
#                            np.array([1, 2, 3, 5, 7, -1, -1, 4])]
#
#         # Test if the returned sample sizes are equal to the expected sample sizes
#         observed_sample_sizes = expand_sample_size_matrix(sample_size, variant_indices)
#         expected_sample_sizes = np.array([133, 133, 156, 156, 104, 75, 52, 156])
#         self.assertTrue(np.array_equal(observed_sample_sizes, expected_sample_sizes))
#
#         # We assume to have 4 phenotypes now
#         # Every cohort has a different set of phenotypes available (and missing)
#         phenotype_indices = [np.array([1, 2, 3, -1]),
#                              np.array([1, 3, 4, 5]),
#                              np.array([-1, 2, 3, 4])]
#
#         # Test if the returned sample sizes are equal to the expected sample sizes
#         observed_sample_sizes = expand_sample_size_matrix(sample_size,
#                                                           variant_indices,
#                                                           phenotype_indices)
#
#         expected_sample_sizes = np.array(
#             [[52, 52, 75, 75, 23, 75, 52, 75],
#              [133, 133, 156, 156, 104, 75, 52, 156],
#              [133, 133, 156, 156, 104, 75, 52, 156],
#              [81, 81, 104, 104, 104, 23, 0, 104]]).T
#
#         self.assertTrue(np.array_equal(observed_sample_sizes, expected_sample_sizes))
#
#     def test_a_inverse_with_missingness(self):
#
#         variant_indices = [np.array([1, 3, 4, 6, -1, 9, 10, 12]),
#                            np.array([-1, -1, 3, 5, 6, 8, -1, 11]),
#                            np.array([1, 2, 3, 5, 7, -1, -1, 4])]
#
#         constant_a = [np.array([[52., 57.],
#                                 [57., 127.]]),
#                       np.array([[23., 30.],
#                                 [30., 68.]]),
#                       np.array([[81., 110.],
#                                 [110., 256.]])]
#
#         variant_dependent_a = np.array([[[48., 50., 78.],
#                                          [58., 64., 98.],
#                                          [55., 61., 101.],
#                                          [49., 47., 75.],
#                                          [48., 52., 72.],
#                                          [52., 55., 88.],
#                                          [53., 70., 87.],
#                                          [52., 60., 86.]],
#                                         [[25., 41., 41.],
#                                          [27., 37., 49.],
#                                          [23., 25., 39.],
#                                          [27., 39., 47.],
#                                          [27., 32., 49.],
#                                          [26., 38., 44.],
#                                          [26., 39., 46.],
#                                          [24., 31., 40.]],
#                                         [[63., 82., 101.],
#                                          [70., 93., 118.],
#                                          [101., 137., 171.],
#                                          [80., 98., 134.],
#                                          [80., 111., 132.],
#                                          [80., 106., 132.],
#                                          [84., 123., 132.],
#                                          [90., 107., 158.]]])
#
#         variant_dependent_a[np.array(variant_indices) == -1] = 0
#         variant_dependent_a = variant_dependent_a.sum(axis=0)
#
#         a_inverse = a_inverse_extended_allow_missingness(variant_dependent_a=variant_dependent_a,
#                                                          constant_a=constant_a,
#                                                          variant_indices=variant_indices)
#
#         a_inverse_expected = np.array(
#             [[[0.0257807, -0.00768435, -0.01032025],
#               [-0.00768435, 0.00579112, 0.00049461],
#               [-0.01032025, 0.00049461, 0.01162156]],
#
#              [[0.02717178, -0.00747424, -0.01066913],
#               [-0.00747424, 0.00577504, 0.00023158],
#               [-0.01066913, 0.00023158, 0.01078375]],
#
#              [[0.02718515, -0.00641103, -0.01104979],
#               [-0.00641103, 0.00494715, 0.00014263],
#               [-0.01104979, 0.00014263, 0.009473]],
#
#              [[0.02608415, -0.00694555, -0.01090292],
#               [-0.00694555, 0.00498668, 0.00064827],
#               [-0.01090292, 0.00064827, 0.01008427]],
#
#              [[0.03821472, -0.01004418, -0.01465556],
#               [-0.01004418, 0.00737881, 0.00010806],
#               [-0.01465556, 0.00010806, 0.01410328]],
#
#              [[0.04766925, -0.01179808, -0.01985591],
#               [-0.01179808, 0.01064337, -0.00052715],
#               [-0.01985591, -0.00052715, 0.0196802]],
#
#              [[0.05953603, -0.0120931, -0.02653899],
#               [-0.0120931, 0.01660502, -0.0059933],
#               [-0.02653899, -0.0059933, 0.03248388]],
#
#              [[0.02640711, -0.00685729, -0.01065435],
#               [-0.00685729, 0.00497599, 0.00053896],
#               [-0.01065435, 0.00053896, 0.00937292]]])
#
#         self.assertTrue(np.allclose(a_inverse, a_inverse_expected))
#

class TestPermutationStrategies(unittest.TestCase):
    number_of_individuals = 81
    number_of_variants = 8

    def test_without_interactions(self):
        """
        We want to test here what a permutation strategy could look like using
        the encoded matrices in hase.
        """
        genotype_matrix = get_genotype_matrix(self.number_of_individuals,
                                              self.number_of_variants)
        covariates_matrix = get_covariates_matrix(self.number_of_individuals, 1)  # One covariate
        covariates_matrix_extended = covariates_matrix.repeat(self.number_of_variants, axis=1)
        base_x = np.stack((covariates_matrix_extended, genotype_matrix.T), axis=1)
        matrix_with_determinants_x = get_x_no_interaction(base_x)

        phenotype_matrix = get_phenotype_matrix_no_interaction(matrix_with_determinants_x)[:, 0:10]

        # Encode the data
        encoder = EncoderCopy(self.number_of_individuals)
        genotype_matrix_encoded = encoder.encode_genotype_matrix(genotype_matrix)
        phenotype_matrix_encoded = encoder.encode_with_inverse(phenotype_matrix)

        indices = np.arange(self.number_of_individuals)
        permutation_indices = np.random.permutation(indices)

        permuted_phenotype_matrix_encoded = encoder.encode_with_inverse(phenotype_matrix[permutation_indices,:])
        B4_new = B4(phenotype_matrix_encoded[permutation_indices,:], genotype_matrix_encoded)
        B4_theoretical = B4(permuted_phenotype_matrix_encoded, genotype_matrix_encoded)


def get_single_study_data(number_of_individuals, number_of_variants):
    genotype_matrix = get_genotype_matrix(number_of_individuals,
                                          number_of_variants)
    covariates_matrix = get_covariates_matrix(number_of_individuals, 1)  # One covariate
    covariates_matrix_extended = covariates_matrix.repeat(number_of_variants, axis=1)
    base_x = np.stack((covariates_matrix_extended, genotype_matrix.T), axis=1)
    matrix_with_determinants_x = get_x_no_interaction(base_x)

    phenotype_matrix = get_phenotype_matrix_no_interaction(matrix_with_determinants_x)[:, 0:2]

    # Encode the data
    encoder = EncoderCopy(number_of_individuals)
    genotype_matrix_encoded = encoder.encode_genotype_matrix(genotype_matrix)
    phenotype_matrix_encoded = encoder.encode_with_inverse(phenotype_matrix)

    a_cov, variant_dependent_a = get_a_parts_alternative(covariates_matrix, genotype_matrix)

    b_cov = B_covariates(covariates_matrix, phenotype_matrix)

    C = C_matrix(phenotype_matrix)

    b4 = B4(phenotype_matrix_encoded, genotype_matrix_encoded)

    return a_cov, variant_dependent_a, b_cov, C, b4, matrix_with_determinants_x, phenotype_matrix


def get_genotype_matrix(number_of_individuals, number_of_variants):
    return np.random.randint(
        0, 3, number_of_variants * number_of_individuals)\
        .reshape(number_of_variants, number_of_individuals)


def get_covariates_matrix(number_of_individuals, number_of_covariates):
    return np.random.randint(0, 4, number_of_individuals * number_of_covariates)\
        .reshape(number_of_individuals, number_of_covariates)


def get_phenotype_matrix_with_interaction(X):
    reshape = np.random.normal(0, 2, X.shape[0] * (X.shape[2])).reshape(X.shape[0], X.shape[2])
    return X[:,5] * 4 + X[:,3] * 0.5 + X[:,1] * 3 + 10 + reshape


def get_phenotype_matrix_no_interaction(X):
    reshape = np.random.normal(0, 2, X.shape[0] * (X.shape[2])).reshape(X.shape[0], X.shape[2])
    return X[:,2] * 4 + X[:,1] * 3 + 10 + reshape


def fit_model(X2, y):
    X2 = X2[:,1:]
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

    p_values = [2 * (1 - stats.t.cdf(np.abs(i), (len(newX) - 1))) for i in ts_b]

    myDF3 = pd.DataFrame()
    myDF3["Coefficients"], myDF3["Standard Errors"], myDF3["t values"], myDF3["Probabilites"] = [params, sd_b, ts_b,
                                                                                                 p_values]
    return myDF3


def calculate_a(covariates_matrix, genotype_matrix):
    a_covariates = A_covariates(covariates_matrix)
    a_tests = A_tests(covariates_matrix, genotype_matrix)
    return A_inverse(a_covariates, a_tests)[0]


def test_hase2_no_interaction(number_of_individuals,
                              number_of_variants):
    genotype_matrix = get_genotype_matrix(number_of_individuals,
                                          number_of_variants)
    covariates_matrix = get_covariates_matrix(number_of_individuals, 1) # One covariate
    covariates_matrix_extended = covariates_matrix.repeat(number_of_variants, axis=1)
    base_x = np.stack((covariates_matrix_extended, genotype_matrix.T), axis=1)
    matrix_with_determinants_x = get_x_no_interaction(base_x)

    phenotype_matrix = get_phenotype_matrix_no_interaction(matrix_with_determinants_x)[:, 0:2]

    # Encode the data
    encoder = EncoderCopy(number_of_individuals)
    genotype_matrix_encoded = encoder.encode_genotype_matrix(genotype_matrix)
    phenotype_matrix_encoded = encoder.encode_with_inverse(phenotype_matrix)

    # Calculate the A inverse the new way
    a_inverse_alternative = calculate_a_alternative(covariates_matrix, genotype_matrix)

    # Calculate the A inverse the old way
    a_inverse = calculate_a(covariates_matrix, genotype_matrix)

    b_cov = B_covariates(covariates_matrix, phenotype_matrix)

    C = C_matrix(phenotype_matrix)

    b4 = B4(phenotype_matrix_encoded, genotype_matrix_encoded)

    b_variable = b4[np.newaxis, ...]

    number_of_variable_terms = b_variable.shape[0]

    N_con = a_inverse_alternative.shape[1] - number_of_variable_terms

    DF = (number_of_individuals - a_inverse_alternative.shape[1])

    t_stat, SE = HASE(b4, a_inverse, b_cov, C, N_con, DF)

    t_stat2, SE2 = hase_supporting_interactions(b_variable, a_inverse_alternative, b_cov, C, N_con, DF)

    assert(np.allclose(SE, SE2))
    assert(np.allclose(t_stat, t_stat2))
    print("Standard error and t-statistics equal between old model and new model")

    model_table = fit_model(matrix_with_determinants_x[..., 0], phenotype_matrix[..., 0])

    assert(np.isclose(model_table.loc[a_inverse_alternative.shape[1] - 1, u"Standard Errors"], SE[0, 0, 0]))
    assert(np.isclose(model_table.loc[a_inverse_alternative.shape[1] - 1, u"t values"], t_stat[0, 0, 0]))
    print("Standard error and t-statistics equal between model and regular regression analysis")

class TestSlicingB4(unittest.TestCase):
    number_of_individuals = 81
    number_of_variants = 10000

    def test_b4_slicing(self):

        genotype_matrix = get_genotype_matrix(self.number_of_individuals,
                                              self.number_of_variants)
        covariates_matrix = get_covariates_matrix(self.number_of_individuals, 1) # One covariate
        covariates_matrix_extended = covariates_matrix.repeat(self.number_of_variants, axis=1)
        base_x = np.stack((covariates_matrix_extended, genotype_matrix.T), axis=1)
        matrix_with_determinants_x = get_x_no_interaction(base_x)

        phenotype_matrix = get_phenotype_matrix_no_interaction(matrix_with_determinants_x)[:, 0:5]

        # Encode the data
        encoder = EncoderCopy(self.number_of_individuals)
        genotype_matrix_encoded = encoder.encode_genotype_matrix(genotype_matrix)
        phenotype_matrix_encoded = encoder.encode_with_inverse(phenotype_matrix)

        # Calculate the A inverse the new way
        a_inverse_alternative = calculate_a_alternative(covariates_matrix, genotype_matrix)

        # Calculate the A inverse the old way
        a_inverse = calculate_a(covariates_matrix, genotype_matrix)

        b_cov = B_covariates(covariates_matrix, phenotype_matrix)

        C = C_matrix(phenotype_matrix)

        b4_actual = B4(phenotype_matrix_encoded, genotype_matrix_encoded)
        b4_expected = B4(phenotype_matrix, genotype_matrix)
        b4_expected_1b = B4(encoder.encode_with_inverse(phenotype_matrix),
                           encoder.encode_genotype_matrix(genotype_matrix))

        self.assertTrue(np.allclose(b4_expected, b4_actual))
        self.assertTrue(np.allclose(b4_expected_1b, b4_actual))

        mask = np.ones(self.number_of_individuals, bool)
        mask[44] = False

        encoder_80 = EncoderCopy(self.number_of_individuals-1)

        b4_expected_2 = B4(phenotype_matrix[mask, :],
                           genotype_matrix[:, mask])
        b4_actual_2 = B4(encoder_80.encode_with_inverse(phenotype_matrix[mask, :]),
                         encoder_80.encode_genotype_matrix(genotype_matrix[:, mask]))

        self.assertTrue(np.allclose(b4_expected_2, b4_actual_2))

        for i in range(self.number_of_individuals):
            mask_for_expected = np.ones(self.number_of_individuals, bool)
            mask_for_expected[i] = False
            b4_expected = B4(encoder_80.encode_with_inverse(phenotype_matrix[mask_for_expected, :]),
                             encoder_80.encode_genotype_matrix(genotype_matrix[:, mask_for_expected]))

            print(i)

            for i2 in range(self.number_of_individuals):
                mask_for_actual = np.ones(self.number_of_individuals, bool)
                mask_for_actual[i2] = False
                b4_actual = B4(phenotype_matrix_encoded[mask_for_actual, :],
                               genotype_matrix_encoded[:, mask_for_actual])

                self.assertFalse(np.allclose(b4_expected, b4_actual))


def get_a_parts_alternative(covariates_matrix, genotype_matrix):
    # Get empty matrix of covariate values
    factor_matrix = np.empty((0, 0))
    # Get the constant part of A
    a_cov = A_covariates(covariates_matrix)
    variant_dependent_a = calculate_variant_dependent_a(genotype_matrix,
                                                        factor_matrix,
                                                        covariates_matrix)
    return a_cov, variant_dependent_a


def calculate_a_alternative(covariates_matrix, genotype_matrix):
    a_cov, variant_dependent_a = get_a_parts_alternative(covariates_matrix, genotype_matrix)
    a_inverse_alternative = get_a_inverse_extended(a_cov, variant_dependent_a)
    return a_inverse_alternative


def get_x_no_interaction(base_x):
    x = np.stack([preprocessing.add_dummy_feature(base_x[..., i]) for i in range(base_x.shape[-1])],
                 axis=2)
    return x


def get_x_with_interaction(base_X):
    poly = PolynomialFeatures(interaction_only=True)
    # indices = list(range(base_X.shape[1]))
    # indices.extend(range(indices[-1]))
    # indices.append(base_X.shape[1])
    # print(indices)
    X = np.stack([poly.fit_transform(base_X[..., i])[:, [0, 1, 2, 5, 6, 3]] for i in range(base_X.shape[-1])],
                             axis=2)
    return X


def test_hase2_with_interaction(number_of_individuals, number_of_variants):
    # Get a genotype matrix
    genotype_matrix = get_genotype_matrix(number_of_individuals,
                                          number_of_variants)
    # Get the matrix of covariates.single
    covariates_matrix = get_covariates_matrix(number_of_individuals, 2) # two covariates.single
    # Get a matrix of covariates.single that is extended to be the same shape as the
    # genotype matrix
    covariates_matrix_extended = np.tile(covariates_matrix[:, np.newaxis, :], (1, number_of_variants, 1))
    # Create the X matrix with base determinants genotype and the covariate
    base_X = np.concatenate((covariates_matrix_extended, genotype_matrix.T[:, :, np.newaxis]), axis=2).transpose((0, 2, 1))
    # Get the X matrix with the
    X = get_x_with_interaction(base_X)
    # Get the phenotype matrix
    phenotype_matrix = get_phenotype_matrix_with_interaction(X)[:, 0:2]

    # Encode the stuff
    encoder = EncoderCopy(number_of_individuals)
    genotype_matrix_encoded = encoder.encode_genotype_matrix(genotype_matrix)
    phenotype_matrix_encoded = encoder.encode_with_inverse(phenotype_matrix)

    encoded_interaction_phenotype_matrices = list()
    # Calculate the phenotypes * interaction covariate
    for covariate_array in covariates_matrix.T:
        interaction_phenotype_matrix = np.einsum('ij,i->ij', phenotype_matrix, covariate_array)
        encoded_interaction_phenotype_matrices.append(encoder.encode_with_inverse(
            interaction_phenotype_matrix))

    # Get the constant part of A
    a_cov = A_covariates(covariates_matrix)
    variant_dependent_A = calculate_variant_dependent_a(genotype_matrix,
                                                        covariates_matrix,
                                                        covariates_matrix)

    # Get the b part that corresponds to the interaction values
    b_interaction_list = list()
    for encoded_interaction_phenotype_matrix in encoded_interaction_phenotype_matrices:
        b_interaction_list.append(np.dot(genotype_matrix_encoded, encoded_interaction_phenotype_matrix))
    b_interaction = np.array(b_interaction_list)

    # Get the constant part of b
    b_cov = B_covariates(covariates_matrix, phenotype_matrix)

    C = C_matrix(phenotype_matrix)

    b4 = B4(phenotype_matrix_encoded, genotype_matrix_encoded)

    a_inv, _ = get_a_inverse_extended(a_cov, variant_dependent_A)

    b_variable = np.append(b_interaction, b4[np.newaxis,...], axis=0)

    number_of_variable_terms = b_variable.shape[0]

    N_con = a_inv.shape[1] - number_of_variable_terms

    DF = (number_of_individuals - a_inv.shape[1])

    t_stat, SE = hase_supporting_interactions(b_variable, a_inv, b_cov, C, N_con, DF)

    model_table = fit_model(X[..., 0], phenotype_matrix[..., 0])

    # Assert whether the HASE2 method returns equal results compared to the regular
    # regression analyses.
    assert(np.isclose(model_table.loc[a_inv.shape[1] - 1, u"Standard Errors"], SE[0, 0, 0]))
    assert(np.isclose(model_table.loc[a_inv.shape[1] - 1, u"t values"], t_stat[0, 0, 0]))
    print("Standard error and t-statistics equal between model and regular regression analysis")


def main(argv=None):
    if argv is None:
        argv = sys.argv

    number_of_individuals = 8
    number_of_variants = 3

    #test_hase2_no_interaction(number_of_individuals, number_of_variants)

    #test_hase2_with_interaction(number_of_individuals, number_of_variants)
    test_b4_slicing(number_of_individuals = 20, number_of_variants = 10)
    return 0


if __name__ == "__main__":
    sys.exit(main())