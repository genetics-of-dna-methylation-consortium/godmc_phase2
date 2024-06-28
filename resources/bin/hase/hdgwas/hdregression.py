from __future__ import print_function
import numpy as np
import os
import sys
from hdgwas.tools import Timer, timer, timing, save_parameters
import scipy.linalg.blas as FB
import h5py
import gc
import tables


# @timing
def A_covariates(covariates, intercept=True):
    '''
    :param covariates: (n_subjects, n_covariates) - only constant covariates.single should be included (age, sex, ICV etc)
    :param intercept: default True, add intercept to model
    :return: matrix (n_cavariates, n_covariates), constant part for the rest of the study
    '''

    S, N = covariates.shape
    if intercept:
        I = np.ones(S).reshape(S, 1)
        covariates = np.hstack((I, covariates))
    a_cov = np.dot(covariates.T, covariates)
    return a_cov


# @timing
def B4(phenotype, genotype):
    b4 = np.tensordot(genotype, phenotype, axes=([1], [0]))
    return b4


def interaction(genotype, factor):
    g = genotype * factor.T
    return g


def calculate_interaction_b(genotype, interaction_values):
    """
    Decodes genotype and interaction values into the sum of Y[i] * COV[i] * SNP[i]
    :param genotype: The encoded genotype matrix
    :param interaction_values: The encoded
    :return:
    """
    b_interactions = np.tensordot(genotype, interaction_values, axes=([1], [0]))
    return b_interactions


# @timing
def A_tests(covariates, genotype, intercept=True):  # TODO (low) extend for any number of tests in model
    '''
    :param covariates: (n_subjects, n_covariates) - only constant covariates.single should be included (age, sex, ICV etc)
    :param genotype: (n_tests, n_subjects) - test could be any kind of quantitative covariance
    :return: (1,n_covariates + intercept)
    '''

    if intercept:
        fst = np.sum(genotype, axis=1).reshape(-1, 1)
        sec = np.dot(genotype, covariates)
        tr = np.sum(np.power(genotype, 2), axis=1).reshape(-1, 1)
        return np.hstack((fst, sec, tr))

    else:
        sec = np.dot(genotype, covariates)
        tr = np.sum(np.power(genotype, 2), axis=1).reshape(-1, 1)
        return np.hstack((sec, tr))


def calculate_variant_dependent_a(genotype, factor_matrix,
                                  covariates, intercept=True):
    """
    Function that calculates the part of A that varies between variants.
    The output is a 2d array for every variant in the genotype matrix, with
    the values representing the sum of the products of the determinant on the row
    and the non-constant determinant on the column for every sample
    For instance, when no interaction stuff is required and intercept is true,
    the variable dependent a for every variant is equivalent to the following.
    [[ sum ( genotypes[variant_index] * intercept value ) ]
     [ sum ( genotypes[variant_index] * covariate values ) ]
     [ sum ( genotypes[variant_index] * genotypes[variant_index] ) ]

    When interaction values are used, the variable dependent a for every variant
    is equivalent to the following, where interaction_values = interaction_factor * genotypes
    [[ sum (                   interaction_values[variant_index] * intercept value ) , sum ( genotypes[variant_index] * intercept value ) ]
     [ sum (                  interaction_values[variant_index] * covariate values ) , sum ( genotypes[variant_index] * covariate values ) ]
     [ sum ( interaction_values[variant_index] * interaction_values[variant_index] ) , sum ( genotypes[variant_index] * interaction_values[variant_index] ) ]
     [ sum (                                                                       0 , sum ( genotypes[variant_index] * gentoypes[variant_index] ) ]]

    :param genotype: A 2d matrix
    :param factor_matrix: A 2d matrix for the interaction values
    :param covariates: A 2d matrix with covariates.single
    :param intercept: If the intercept is true.
    :return: 3d array with variable dependent a (see above.)
    """
    number_of_variable_terms = (factor_matrix.shape[1] + 1)
    number_of_total_terms = (covariates.shape[1] + number_of_variable_terms)
    if intercept:
        number_of_total_terms += 1

    variant_dependent_a = np.zeros(
        number_of_variable_terms * number_of_total_terms * genotype.shape[0]).reshape((
        number_of_variable_terms, number_of_total_terms, genotype.shape[0]))
    variable_term_index = 0
    # Have to add extra columns to the covariates.single

    #covariates = np.tile(covariates, (genotype.shape[0],1,1))

    # Loop through the columns of the factor matrix
    for factor_column in factor_matrix.T:
        raise NotImplementedError("Interactions not yet implemented")
        # Multiply the genotypes with the factor column
        interaction_values = genotype * factor_column
        # This creates a 2d array with columns representing the
        # variants and values representing the covariates.single from individuals

        # Calculate the A values for interaction values and other independent
        # determinants already in covariates.single matrix
        sec = calculate_dot_product_for_variants(covariates, interaction_values)
        tr = np.sum(np.power(interaction_values, 2), axis=1).reshape(-1, 1)

        if intercept:
            fst = np.sum(interaction_values, axis=1).reshape(-1, 1)
            variable_term_values = np.hstack((fst, sec, tr)).T
            variant_dependent_a[variable_term_index, 0:variable_term_values.shape[0]] = variable_term_values
        else:
            variable_term_values = np.hstack((sec, tr)).T
            variant_dependent_a[variable_term_index, 0:variable_term_values.shape[0]] = variable_term_values

        variable_term_index += 1
        covariates = np.dstack((covariates, interaction_values))

    # Calculate the A values for genotypes with the other independent determinants
    sec = calculate_dot_product_for_variants(covariates, genotype)
    tr = np.sum(np.power(genotype, 2), axis=1).reshape(-1, 1)

    if intercept:
        fst = np.sum(genotype, axis=1).reshape(-1, 1)
        variable_term_values = np.hstack((fst, sec, tr)).T
        variant_dependent_a[variable_term_index, 0:variable_term_values.shape[0]] = variable_term_values
    else:
        variable_term_values = np.hstack((sec, tr)).T
        variant_dependent_a[variable_term_index, 0:variable_term_values.shape[0]] = variable_term_values
    if variant_dependent_a.shape[0] == 1:
        variant_dependent_a = np.squeeze(variant_dependent_a, axis=0)
    return variant_dependent_a.T


def calculate_dot_product_for_variants(covariates, other_independent_determinant):
    # In the dot einsum notation, labels represent the following:
    # v: variants
    # i: individuals (dot product of genotypes, covariates.single)
    # k: different covariates.single
    sec = np.einsum('vi,ik->vk', other_independent_determinant, covariates)
    # (Values get summed along individuals)
    return sec


def expand_B_covariates(b_cov, variant_indices):
    pass


def a_inverse_extended_allow_missingness(variant_dependent_a, constant_a, variant_indices, random_effect_intercept):
    variant_indices_arr = np.array(variant_indices).T

    # Obtain the dimensions of an expanded constant_a
    number_of_variants = variant_indices[0].shape[0]

    # Define empty array of dimensions nVariants, dim(constant_a)
    a_covariates_expanded = np.repeat(
        np.array(constant_a)[np.newaxis, ...],
        number_of_variants, axis = 0)

    a_covariates_expanded[variant_indices_arr == -1, :, :] = 0

    # For every variant, sum the A_cov's for every study
    #a_covariates_expanded = np.sum(a_covariates_expanded, axis = 1)

    a_inverse_list = []

    n, m = constant_a[0].shape
    k = n + variant_dependent_a.shape[2]

    for i in range(variant_dependent_a.shape[0]):
        a_complete = np.zeros(k*k).reshape(k,k)
        a_complete[0:n, 0:n] = a_covariates_expanded[i, ...]
        a_complete[:, n:k] = variant_dependent_a[i, :]
        a_complete[n:k, :k-1] = np.maximum(a_complete[n:k, :k-1],
                                    variant_dependent_a[i,:k-1].T)
        try:
            a_inverse_list.append(np.linalg.inv(a_complete))
        except:
            a_inverse_list.append(np.zeros(k * k).reshape(k, k))  # TODO (high) test; check influence on results; warning;
    return np.array(a_inverse_list)


def expand_C_matrix(c_stack, variant_indices):
    pass


def expand_sample_size_matrix(sample_size, variant_indices, phenotype_indices=None):
    """
    Function that expands sample sizes per cohort to sample sizes per variant x phenotype combination.

    :param sample_size: List of sample size per cohort.
    :param variant_indices: List of np.array objects listing variant indices (-1 when missing)
    :param phenotype_indices: List of np.array objects listing phenotype indices (-1 when missing)
    :return: matrix of dimensions nVariants * nPhenotypes, with every combination portraying the number
    of samples available to test that combination.
    """
    # Get the number of variants. (across 0th dimension, should be equal for all cohorts in the list)
    number_of_variants = variant_indices[0].shape[0]

    # Expand the sample size across all variants across a new 1st axis
    sample_size_expanded = np.array(sample_size)[:, np.newaxis].repeat(number_of_variants, 1)

    # Set cohort sample size to 0 for variants that are missing (-1) for that cohort
    sample_size_expanded[np.array(variant_indices) == -1] = 0

    # Expand the sample size across all phenotypes if the indices are not none
    if phenotype_indices is not None:

        # Get the number of phenotypes. (across 0th dimension, should be equal for all cohorts in the list)
        number_of_phenotypes = phenotype_indices[0].shape[0]

        # Expand the sample size across all phenotypes across a new 2nd axis
        sample_size_expanded = sample_size_expanded[:, :, np.newaxis].repeat(number_of_phenotypes, 2)

        # Set cohort sample size to 0 for phenotypes that are missing (-1) for that cohort
        sample_size_expanded[np.array(phenotype_indices)[:, np.newaxis, :].repeat(number_of_variants, 1) == -1] = 0

    return np.sum(sample_size_expanded, 0)


# @timing
def B_covariates(covariates, phenotype, intercept=True):
    S, N = covariates.shape

    b_cov = np.dot(covariates.T, phenotype)
    if intercept:
        b1 = np.sum(phenotype, axis=0).reshape(1, phenotype.shape[1])
        B13 = np.append(b1, b_cov, axis=0)
        return B13
    else:
        return b_cov


def get_a_inverse_extended(a_covariates, variant_dependent_a):
    n,m = a_covariates.shape
    if variant_dependent_a.ndim == 2:
        return A_inverse(a_covariates, variant_dependent_a)
    k = n + variant_dependent_a.shape[2]

    a_complete = np.zeros((variant_dependent_a.shape[0], k, k))
    a_complete[..., 0:n, 0:n] = a_covariates
    a_complete[..., :, n:k] = variant_dependent_a
    a_complete = np.triu(a_complete) + np.triu(a_complete, 1).transpose((0,2,1))

    return(np.linalg.inv(a_complete))


# @timing
def A_inverse(a_covariates, a_test):  # TODO (low) extend for any number of tests in model

    A_inv = []
    n, m = a_covariates.shape
    k = n + 1
    for i in xrange(a_test.shape[0]):  # TODO (low) not in for loop
        inv = np.zeros(k * k).reshape(k, k)
        inv[0:k - 1, 0:k - 1] = a_covariates
        inv[k - 1, :] = a_test[i, :]
        inv[0:k, k - 1] = a_test[i, 0:k]
        try:
            A_inv.append(np.linalg.inv(inv))
        except:
            A_inv.append(np.zeros(k * k).reshape(k, k))  # TODO (high) test; check influence on results; warning;

    return np.array(A_inv)


# @timing
def C_matrix(phenotype):
    C = np.einsum('ij,ji->i', phenotype.T, phenotype)
    return C


# @timing
# @save_parameters
def HASE(b4, A_inverse, b_cov, C, N_con, DF):
    with Timer() as t:
        # These together form the X matrix
        B13 = b_cov
        B4 = b4

        A1_B_constant = np.tensordot(A_inverse[:, :, 0:(N_con)], B13, axes=([2], [0]))

        A1_B_nonconstant = np.einsum('ijk,il->ijl', A_inverse[:, :, N_con:N_con + 1], B4)

        # Combine the inverse of A multiplied with
        # constant part of B and non-constant part of B
        A1_B_full = A1_B_constant + A1_B_nonconstant

        BT_A1B_const = np.einsum('ij,lji->li', B13.T, A1_B_full[:, 0:(N_con), :])

        BT_A1B_nonconst = np.einsum('ijk,ijk->ijk', B4[:, None, :], A1_B_full[:, (N_con):N_con + 1, :])

        BT_A1B_full = BT_A1B_const[:, None, :] + BT_A1B_nonconst

        C_BTA1B = BT_A1B_full - C.reshape(1, -1)

        C_BTA1B = np.abs(C_BTA1B)

        a44_C_BTA1B = C_BTA1B * A_inverse[:, (N_con):N_con + 1, (N_con):N_con + 1]

        a44_C_BTA1B = np.sqrt((a44_C_BTA1B))

        t_stat = np.sqrt(DF) * np.divide(A1_B_full[:, (N_con):N_con + 1, :], a44_C_BTA1B)

        SE = a44_C_BTA1B / np.sqrt(DF)

    print("time to compute GWAS for {} phenotypes and {} SNPs .... {} sec".format(b4.shape[1],
                                                                                  A_inverse.shape[0],
                                                                                  t.secs))
    return t_stat, SE


# @timing
# @save_parameters
def hase_supporting_interactions(b_variable, a_inverse, b_cov, C, number_of_constant_terms, DF):
    with Timer() as t:
        # These together form the X matrix
        B13 = b_cov
        number_of_variable_terms = b_variable.shape[0]
        variant_effect_index = a_inverse.shape[1] - 1

        A1_B_constant = np.tensordot(a_inverse[:, :, 0:(number_of_constant_terms)], B13, axes=([2], [0]))

        # In the einsum notation, the labels represent the following:
        # i: The variant axis
        # j: ...
        # k: The axis with regression terms (interaction, genotype)
        # l: The phenotype axis
        A1_B_nonconstant = np.einsum('ijk,kil->ijl',
                                     a_inverse[:, :, number_of_constant_terms:number_of_constant_terms + number_of_variable_terms],
                                     b_variable)

        # Combine the inverse of A multiplied with
        # constant part of B and non-constant part of B
        A1_B_full = A1_B_constant + A1_B_nonconstant

        # In the einstein summation notation, the labels represent the following:
        # l: The variant axis
        # j: The the axis with the regression terms (interaction, genotype)
        # i: The phenotype axis
        BT_A1B_const = np.einsum('ij,lji->li', B13.T, A1_B_full[:, :(number_of_constant_terms), :])

        # In the einstein summation notation, the labels represent the following
        # i: The variant axis
        # j: Terms
        # k: Phenotype
        BT_A1B_nonconst = np.einsum(
            'ijk,ijk->ik', b_variable.transpose((1, 0, 2)),
            A1_B_full[:, (number_of_constant_terms):number_of_constant_terms + number_of_variable_terms, :])

        # Combine the constant and nonconstant parts of the BT, beta matrix
        BT_A1B_full = BT_A1B_const + BT_A1B_nonconst

        # Get the difference between C (dot product of phenotypes) and the matrix of estimates
        C_BTA1B = BT_A1B_full - C.reshape(1, -1)
        C_BTA1B = np.abs(C_BTA1B)

        # Multiply the far right / lower part of every A_inv matrix
        # (this is the part with the sum of the squares of genotype dosages),
        # with the differences between C and the BT_A1B_full matrix
        a44_C_BTA1B = np.einsum(
            'il,i...->i...l', C_BTA1B,
            a_inverse[:, variant_effect_index, (variant_effect_index):variant_effect_index+1])

        # Square this result
        a44_C_BTA1B = np.sqrt((a44_C_BTA1B))

        # Get the t-statistics
        t_stat = np.sqrt(DF) * np.divide(
            A1_B_full[:, (variant_effect_index):variant_effect_index+1, :], a44_C_BTA1B)
        # The t-statistics are stored in a 3d array with the first dimension
        # representing the variants (main determinants to test), the second
        # dimension being synonymous (or at least as it appears to me) to the
        # third dimension representing phenotypes.
        # Second dimension thus only contains one element (phenotype array)

        # Get the standard error.
        SE = a44_C_BTA1B / np.sqrt(DF)
        # Standard errors are stored in the same ways as the t-statistics:
        # SE[variant][0][phenotype]

    print("time to compute GWAS for {} phenotypes and {} SNPs .... {} sec".format(b_variable.shape[2],
                                                                                  a_inverse.shape[0],
                                                                                  t.secs))
    return t_stat, SE


# @timing
# @save_parameters
def hase_allow_missingness_supporting_interactions(b_variable, a_inverse, b_cov, C, number_of_constant_terms, DF):

    # What do we have:
    # - variable part of the B-matrix (already expanded)
    # - summed A inverse (already expanded)
    # - summed C (already expanded)
    # - summed static part of B

    with Timer() as t:
        # These together form the X matrix
        B13 = b_cov
        number_of_variable_terms = b_variable.shape[0]
        variant_effect_index = a_inverse.shape[1] - 1

        A1_B_constant = np.tensordot(a_inverse[:, :, 0:(number_of_constant_terms)], B13, axes=([2], [0]))

        print(A1_B_constant[3,...])

        # In the einsum notation, the labels represent the following:
        # i: The variant axis
        # j: ...
        # k: The axis with regression terms (interaction, genotype)
        # l: The phenotype axis
        A1_B_nonconstant = np.einsum('ijk,kil->ijl',
                                     a_inverse[:, :, number_of_constant_terms:number_of_constant_terms + number_of_variable_terms],
                                     b_variable)

        # Combine the inverse of A multiplied with
        # constant part of B and non-constant part of B
        A1_B_full = A1_B_constant + A1_B_nonconstant

        # In the einstein summation notation, the labels represent the following:
        # l: The variant axis
        # j: The the axis with the regression terms (interaction, genotype)
        # i: The phenotype axis
        BT_A1B_const = np.einsum('ij,lji->li', B13.T, A1_B_full[:, :(number_of_constant_terms), :])

        # In the einstein summation notation, the labels represent the following
        # i: The variant axis
        # j: Terms
        # k: Phenotype
        BT_A1B_nonconst = np.einsum(
            'ijk,ijk->ik', b_variable.transpose((1, 0, 2)),
            A1_B_full[:, (number_of_constant_terms):number_of_constant_terms + number_of_variable_terms, :])

        # Combine the constant and nonconstant parts of the BT, beta matrix
        BT_A1B_full = BT_A1B_const + BT_A1B_nonconst

        print(BT_A1B_full[3,...])

        # Get the difference between C (dot product of phenotypes) and the matrix of estimates
        C_BTA1B = BT_A1B_full - C.reshape(1, -1)
        C_BTA1B = np.abs(C_BTA1B)

        # Multiply the far right / lower part of every A_inv matrix
        # (this is the part with the sum of the squares of genotype dosages),
        # with the differences between C and the BT_A1B_full matrix
        a44_C_BTA1B = np.einsum(
            'il,i...->i...l', C_BTA1B,
            a_inverse[:, variant_effect_index, (variant_effect_index):variant_effect_index+1])

        # Square this result
        a44_C_BTA1B = np.sqrt((a44_C_BTA1B))

        # Get the t-statistics
        t_stat = np.sqrt(DF) * np.divide(
            A1_B_full[:, (variant_effect_index):variant_effect_index+1, :], a44_C_BTA1B)
        # The t-statistics are stored in a 3d array with the first dimension
        # representing the variants (main determinants to test), the second
        # dimension being synonymous (or at least as it appears to me) to the
        # third dimension representing phenotypes.
        # Second dimension thus only contains one element (phenotype array)

        # Get the standard error.
        SE = a44_C_BTA1B / np.sqrt(DF)
        # Standard errors are stored in the same ways as the t-statistics:
        # SE[variant][0][phenotype]

    print("time to compute GWAS for {} phenotypes and {} SNPs .... {} sec".format(b_variable.shape[2],
                                                                                  a_inverse.shape[0],
                                                                                  t.secs))
    return t_stat, SE