from __future__ import print_function

import sys
import os
import numpy as np

from config import MAPPER_CHUNK_SIZE, basedir, CONVERTER_SPLIT_SIZE, PYTHON_PATH
from hdgwas.meta_classic import CohortAnalyser
from hdgwas.meta_classic import ClassicMetaAnalyser

os.environ['HASEDIR'] = basedir
if PYTHON_PATH is not None:
    for i in PYTHON_PATH: sys.path.insert(0, i)
import h5py
import tables
from hdgwas.tools import Timer, Checker, study_indexes, Mapper, HaseAnalyser, merge_genotype, Reference, timing, \
    check_np, check_converter, get_intersecting_individual_indices, select_identifiers
from hdgwas.converter import GenotypePLINK, GenotypeMINIMAC, GenotypeVCF
from hdgwas.data import Reader, MetaParData, MetaPhenotype
from hdgwas.fake import Encoder
from hdgwas.hdregression import HASE, A_covariates, A_tests, B_covariates, C_matrix, A_inverse, B4, \
    get_a_inverse_extended, hase_supporting_interactions, expand_sample_size_matrix, \
    a_inverse_extended_allow_missingness
import argparse
import gc
from hdgwas.pard import partial_derivatives
from hdgwas.regression import haseregression
import pandas as pd
import time
from hdgwas.protocol import Protocol

__version__ = '1.1.0'

HEAD = "*********************************************************************\n"
HEAD += "* HASE: Framework for efficient high-dimensional association analyses \n"
HEAD += "* Version {V}\n".format(V=__version__)
HEAD += "* (C) 2015-2017 Gennady Roshchupkin and Hieab Adams\n"
HEAD += "* Erasmus MC, Rotterdam /  Department of Medical Informatics, Radiology and Epidemiology \n"
HEAD += "* GNU General Public License v3\n"
HEAD += "*********************************************************************\n"

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    global MAPPER_CHUNK_SIZE

    start = time.time()

    parser = argparse.ArgumentParser(description='Script to use HASE in command line')

    # MAIN SETTINGS
    parser.add_argument("-thr", type=float,
                        help="predefined threshold for t-value, if not defined save everything")

    parser.add_argument("-o", "--out", type=str, required=True, help="path to save result folder")

    parser.add_argument("-mode", required=True, type=str,
                        choices=['regression', 'converting', 'single-meta', 'meta-stage', 'meta-classic', 'encoding'],
                        help='Main argument, specify type of analysis'
                             '*****************'
                             'converting - before start to use any kind of analysis you need to convert data, it could take some time, '
                             'you need to do it ones and then everything will be super fast'
                             '*****************'
                             'single-meta - if you are the part of meta-analysis and will send data to main center, use this mode'
                             '*****************'
                             'encoding - if you are doing HD association, before send data to main center create fake data'
                             '*****************'
                             'meta-stage - if you are the main center of meta-analysis and already got precomputed data from all sites, use it'
                             '*****************'
                             'regression - if you want run everything yourself and do not use precomputed data use it'

                             'If you want to run single site full association analysis run (1) single-meta mode and then (2) meta-stage'

                             '*****************'
                        )

    parser.add_argument("-g", "--genotype", nargs='+', type=str, help="path/paths to genotype data folder")
    parser.add_argument("-ph", "--phenotype", nargs='+', type=str, help="path to phenotype data folder")
    parser.add_argument("-cov", "--covariates", type=str, help="path to covariates data folder")
    parser.add_argument("-ci", "--covariate-indices", type=str,
                        help="path to tab-delimited file with each row representing "
                             "the indices of covariates to include in the classical meta-analysis")

    parser.add_argument("-derivatives", nargs='+', type=str, help="path to derivatives data folder")

    parser.add_argument("-protocol", type=str, help="path to study protocol file")

    parser.add_argument('-study_name', type=str, required=True, nargs='+',
                        help=' Name for saved genotype data, without ext')

    parser.add_argument('-mapper', type=str, help='Mapper data folder')

    parser.add_argument('-ref_name', type=str, default='1000Gp1v3_ref', help='Reference panel name')

    # ADDITIONAL SETTINGS
    parser.add_argument(
        "-snp_id_inc", type=str, nargs='+',
        help="path to file with SNPs id to include to analysis")
    parser.add_argument(
        "-snp_id_exc", type=str, nargs='+',
        help="path to file with SNPs id to exclude from analysis")
    parser.add_argument(
        "-ph_id_inc", type=str, nargs='+',
        help="path to file with phenotype id to exclude from analysis")
    parser.add_argument(
        "-ph_id_exc", type=str, nargs='+',
        help="path to file with phenotype id to exclude from analysis")

    # parser.add_argument("-ind_id_inc", type=str, help="path to file with individuals id to include to analysis") #TODO (low)
    # parser.add_argument("-ind_id_exc", type=str, help="path to file with individuals id to exclude from analysis")#TODO (low)
    # parser.add_argument("-cov_name_inc", type=str, help="path to file with covariates names to include to analysis")#TODO (low)
    # parser.add_argument("-cov_name_exc", type=str, help="path to file with covariates names to exclude from analysis")#TODO (low)

    parser.add_argument('-intercept',
                        type=lambda i: i.lower() in ['y', 'yes'],
                        default=True,
                        help='include intercept to regression, default True {y(es)/...}')

    # parser.add_argument('-intercept', type=str, default='y', choices=['y','n'], help='include intercept to regression, default yes')

    parser.add_argument('-maf', type=float, default=0.0, help='MAF for genetics data')

    parser.add_argument('-encoded', nargs='+', type=int,
                        help='Value per every study, 1 - if encoded, 0 - if not')  # Option not documented
    ###

    # FLAGS
    parser.add_argument('-hdf5', action='store_true', default=True, help='flag for genotype data format')
    parser.add_argument('-id', action='store_true', default=False,
                        help='Flag to convert minimac data to genotype per subject files first (default False)')
    parser.add_argument('-pd_full', action='store_true', default=False, help='For not HD association study')
    parser.add_argument('-effect_intercept', action='store_true', default=False,
                        help='Flag for add study effect to PD regression model')
    parser.add_argument('-permute_ph', action='store_true', default=False, help='Flag for phenotype permutation')
    parser.add_argument('-vcf', action='store_true', default=False, help='Flag for VCF data to convert')
    parser.add_argument('-np', action='store_true', default=True, help='Check BLAS/LAPACK/MKL')
    parser.add_argument('-allow_missingness', action='store_true', default=False,
                        help='Flag for using missingness, instead of using default implementation for meta-analysis.')

    # TODO (low) save genotype after MAF
    ###

    ###CLUSTER SETTING
    parser.add_argument('-cluster', type=str, default='n', choices=['y', 'n'],
                        help=' Is it parallel cluster job, default no')
    parser.add_argument('-node', nargs='+', type=int, help='number of nodes / this node number, example: 10 2 ')
    ###

    # ADVANCED SETTINGS
    parser.add_argument('-interaction', type=str,
                        help='path to file with data for genotype interaction test')  # TODO (low)

    parser.add_argument('-interaction_encoded', type=str, nargs='+',
                        help='path to file with data for genotype interaction test')

    parser.add_argument('-mapper_chunk', type=int, help='Change mapper chunk size from config file')
    ###
    args = parser.parse_args(argv)
    if not args.thr:
        print ('WARNING!!! You did not set threshold for t-value, all results will be saved')
    if args.mapper_chunk:
        MAPPER_CHUNK_SIZE = args.mapper_chunk
    ARG_CHECKER = Checker()
    print(args)
    os.environ['HASEOUT'] = args.out

    if args.cluster == 'y':
        if args.node is not None:
            if args.node[1] > args.node[0]:
                raise ValueError('Node # {} > {} total number of nodes'.format(args.node[1], args.node[0]))

    if not os.path.exists(args.out):
        print("Creating output folder {}".format(args.out))
        os.makedirs(args.out)

    if args.np:
        check_np()

    ################################### CONVERTING ##############################
    if args.mode == 'converting':

        # ARG_CHECKER.check(args,mode='converting')

        R = Reader('genotype')
        R.start(args.genotype[0], vcf=args.vcf)

        with Timer() as t:
            if R.format == 'PLINK':
                G = GenotypePLINK(args.study_name[0], reader=R)
                G.split_size = CONVERTER_SPLIT_SIZE
                G.plink2hdf5(out=args.out)

            elif R.format == 'MINIMAC':
                G = GenotypeMINIMAC(args.study_name[0], reader=R)
                if args.cluster == 'y':
                    G.cluster = True
                G.split_size = CONVERTER_SPLIT_SIZE
                G.MACH2hdf5(args.out, id=args.id)

            elif R.format == 'VCF':
                G = GenotypeVCF(args.study_name[0], reader=R)
                if args.cluster == 'y':
                    G.cluster = True
                G.split_size = CONVERTER_SPLIT_SIZE
                G.VCF2hdf5(args.out)
            else:
                raise ValueError('Genotype data should be in PLINK/MINIMAC/VCF format and alone in folder')

        check_converter(args.out, args.study_name[0])
        print ('Time to convert all data: {} sec'.format(t.secs))

    ################################### ENCODING ##############################

    elif args.mode == 'encoding':

        # ARG_CHECKER.check(args,mode='encoding')
        mapper = Mapper()
        mapper.genotype_names = args.study_name
        mapper.chunk_size = MAPPER_CHUNK_SIZE
        mapper.reference_name = args.ref_name
        mapper.load_flip(args.mapper)
        mapper.load(args.mapper)

        phen = Reader('phenotype')
        phen.start(args.phenotype[0])

        gen = Reader('genotype')
        gen.start(args.genotype[0], hdf5=args.hdf5, study_name=args.study_name[0], ID=False)

        # Try to load a file with values to use for genotype interaction terms if the argument was used.
        interactions = None
        if args.interaction is not None:
            # Read the interactions
            interactions = Reader('interaction')
            interactions.start(args.interaction)

        e = Encoder(args.out)
        e.study_name = args.study_name[0]

        # The following function identifies for every sample what the index of the sample is in each data source.
        # row_index contains the output with in
        # row_index[0] the indices for genotype part of the samples, in
        # row_index[1] the indices for the phenotype part of the samples, and in
        # row_index[2] the indices for the interaction part of the samples if available.
        row_index, intersecting_identifiers = study_indexes(phenotype=phen.folder._data, genotype=gen.folder._data,
                                                            covariates=interactions.folder._data if interactions else None)
        with Timer() as t:

            e.matrix(len(intersecting_identifiers), save=True)
            N_snps_read = 0
            while True:
                with Timer() as t_gen:
                    genotype = gen.get_next()
                    if isinstance(genotype, type(None)):
                        break
                    
                    flip = mapper.flip[args.study_name[0]][N_snps_read:N_snps_read + genotype.shape[0]]
                    N_snps_read += genotype.shape[0]
                    flip_index = (flip == -1)
                    genotype = np.apply_along_axis(lambda x: flip * (x - 2 * flip_index), 0, genotype)
                    genotype = genotype[:, row_index[0]]
                    encode_genotype = e.encode(genotype, data_type='genotype')
                    e.save_hdf5(encode_genotype, os.path.join('encode_genotype'), info=gen.folder,
                                index=row_index[0])
                    encode_genotype = None
                    gc.collect()

                print ('Time to create fake genotype is {}sec'.format(t_gen.secs))

            while True:
                with Timer() as t_phen:
                    phenotype = phen.get_next(
                        index=row_index[1])  # Get the indices of the common ids in the phenotype data (1)
                    if interactions:
                        # Assume that we do not have to get the interaction values in chunks
                        interaction_phenotype_values = interactions.get_next(
                            index=row_index[2])  # use the indices of the common ids in the interaction data (2)
                        if isinstance(interaction_phenotype_values, type(None)):
                            break

                        # For the interaction data, for every factor use for interaction term with genotype data,
                        # multiply the data with the phenotype. After this, encode the data using the same matrix that
                        # was used for the phenotype data.
                        for i in range(interaction_phenotype_values.shape[1]):
                            # np.multiply broadcasts the vector of interaction values into a 2d array with the same
                            # dimensions as the phenotype data. This results in the correct 2d array to encode using
                            # the inverse of F
                            encode_product_of_phenotype_and_interaction_values = e.encode(
                                np.einsum('ij,i->ij', phenotype, interaction_phenotype_values[:, i]),
                                data_type='phenotype')
                            if interactions.folder.format == '.npy':
                                e.save_npy(encode_product_of_phenotype_and_interaction_values,
                                           save_path=os.path.join('encode_interaction',
                                                                  interactions.folder._data.names[i]),
                                           info=phen.folder, index=row_index[2])
                            if interactions.folder.format in ['.csv', '.txt']:
                                e.save_csv(encode_product_of_phenotype_and_interaction_values,
                                           save_path=os.path.join('encode_interaction',
                                                                  interactions.folder._data.names[i]),
                                           info=phen.folder, index=row_index[2])

                    if isinstance(phenotype, type(None)):
                        break

                    encode_phenotype = e.encode(phenotype, data_type='phenotype')

                    if phen.folder.format == '.npy':
                        e.save_npy(encode_phenotype, save_path=os.path.join('encode_phenotype'),
                                   info=phen.folder, index=row_index[1])
                    if phen.folder.format in ['.csv', '.txt']:
                        e.save_csv(encode_phenotype, save_path=os.path.join('encode_phenotype'),
                                   info=phen.folder, index=row_index[1])
                    encode_phenotype = None
                    gc.collect()

            if phen.folder.format == '.npy':
                np.save(os.path.join(os.path.join(args.out, 'encode_phenotype', 'info_dic.npy')), e.phen_info_dic)
            print ('Time to create fake phenotype is {}sec'.format(t_phen.secs))

        print ('Time to encode all data: {} sec'.format(t.secs))

    ################################### SINGLE META STAGE ##############################

    elif args.mode == 'single-meta':

        # ARG_CHECKER.check(args,mode='single-meta')
        mapper = Mapper()
        mapper.genotype_names = args.study_name
        mapper.chunk_size = MAPPER_CHUNK_SIZE
        mapper.reference_name = args.ref_name
        mapper.load_flip(args.mapper)
        mapper.load(args.mapper)
        mapper.cluster = args.cluster
        mapper.node = args.node

        phen = Reader('phenotype')
        phen.start(args.phenotype[0])

        cov = Reader('covariates')
        cov.start(args.covariates)

        if cov.folder.n_files > 1:
            raise ValueError('In covariates folder should be only one file!')

        # Try to load a file with values to use for genotype interaction terms if the argument was used.
        interaction = None
        if args.interaction is not None:
            # Read the interactions
            interaction = Reader('interaction')
            interaction.start(args.interaction)
            if interaction.folder.n_files > 1:
                raise ValueError('In interaction folder should be only one file!')

        gen = Reader('genotype')
        gen.start(args.genotype[0], hdf5=args.hdf5, study_name=args.study_name[0], ID=False)

        with Timer() as t:
            partial_derivatives(save_path=args.out, COV=cov, PHEN=phen, GEN=gen, INTERACTION=interaction,
                                MAP=mapper, MAF=args.maf, R2=None, B4_flag=args.pd_full,
                                study_name=args.study_name[0], intercept=args.intercept)
        print ('Time to compute partial derivatives : {} sec'.format(t.secs))

    ################################### MULTI META STAGE ##############################

    elif args.mode == 'meta-stage':

        # ARG_CHECKER.check(args,mode='meta-stage')

        ##### Init data readers #####
        if args.derivatives is None:
            raise ValueError('For meta-stage analysis partial derivatives data are required!')
        mapper = load_mapper(
            mapper_chunk_size=MAPPER_CHUNK_SIZE, study_name=args.study_name, ref_name=args.ref_name,
            mapper_folder=args.mapper, encoded=args.encoded, cluster=args.cluster, node=args.node,
            snp_id_inc=args.snp_id_inc, snp_id_exc=args.snp_id_exc)

        Analyser = HaseAnalyser()

        # It appears to me that the pard list contains
        # partial derivatives matrices for every study
        partial_derivatives_folders = []

        with Timer() as t:
            for i, j in enumerate(args.derivatives):
                partial_derivatives_folders.append(Reader('partial'))
                partial_derivatives_folders[i].start(j, study_name=args.study_name[i])
                partial_derivatives_folders[i].folder.load()

        print("Time used to load partial derivatives is {}s".format(t.secs))

        # Create a list containing booleans representing the presence
        # of the b4 data in the partial derivatives folders
        b4_presence_per_study = [False if isinstance(i.folder._data.b4, type(None)) else True for i in
                                 partial_derivatives_folders]

        # The if statement roughly translates to:
        # if one study contains b4 data in the partial derivatives folder,
        # all should have b4 data, otherwise an error is raised.
        if np.sum(b4_presence_per_study) != len(partial_derivatives_folders) and np.sum(b4_presence_per_study) != 0:
            raise ValueError('All studies should have b4 data for partial derivatives!')

        if args.protocol is not None:
            protocol = Protocol(args.protocol)
        else:
            protocol = None

        meta_pard = MetaParData(partial_derivatives_folders, args.study_name, protocol=protocol)
        encoded_interactions = None

        is_no_b4_present_in_partial_derivatives = np.sum(b4_presence_per_study) == 0
        if is_no_b4_present_in_partial_derivatives:
            phen = []

            with Timer() as t:
                for i, j in enumerate(args.phenotype):
                    phen.append(Reader('phenotype'))
                    phen[i].start(j)
            print("Time to set pheno {} s".format(t.secs))

            #I edited this - Olalekan
            #meta_phen = MetaPhenotype(phen, include=args.ph_id_inc[0], exclude=args.ph_id_exc[0])

            meta_phen = MetaPhenotype(phen, include=args.ph_id_inc, exclude=args.ph_id_exc)

            N_studies = len(args.genotype)

            gen = []
            with Timer() as t:
                for i, j in enumerate(args.genotype):
                    gen.append(Reader('genotype'))
                    gen[i].start(j, hdf5=args.hdf5, study_name=args.study_name[i], ID=False)
            print("Time to set gen {}s".format(t.secs))

            # Create a dictionary with the datasets for obtaining
            # the indices of shared identifiers
            datasets = {"phenotype": tuple(i.folder._data for i in phen),
                        "genotype": tuple(i.folder._data for i in gen),
                        "partial_derivatives": tuple(i.folder._data.metadata for i in partial_derivatives_folders)}

            # If interaction folders are supplied, read these.
            if args.interaction_encoded:
                encoded_interaction_folders = []
                with Timer() as t:
                    for i, j in enumerate(args.interaction_encoded):
                        encoded_interaction_folders.append(Reader('interaction_folder'))
                        # Have to check if these are the correct arguments for .start(...)
                        encoded_interaction_folders[i].start(j)
                datasets["interaction"] = tuple(i.folder._data for i in encoded_interaction_folders)
                encoded_interactions = MetaPhenotype(encoded_interaction_folders, include=args.ph_id_inc[0],
                                                     exclude=args.ph_id_exc[0])
                print("Time used to load encoded interaction data: {}s".format(t.secs))

            # Get common ids
            row_index, intersecting_identifiers = get_intersecting_individual_indices(datasets)

            # Do something for all covariates
            if row_index["partial_derivatives"].shape[0] != np.sum(
                    [i.folder._data.metadata['id'].shape[0] for i in partial_derivatives_folders]):
                raise ValueError(
                    'Partial Derivatives covariates have different number of subjects {} than genotype and phenotype {}'.format(
                        row_index["covariates"].shape[0],
                        np.sum([i.folder._data.metadata['id'].shape[0] for i in partial_derivatives_folders])))

        # Start looping over all genotype chunks
        while True:
            if mapper.cluster == 'n':
                variant_indices, keys = mapper.get(
                    allow_missingness=args.allow_missingness)
            else:
                ch = mapper.chunk_pop()
                if ch is None:
                    variant_indices = None
                    break
                variant_indices, keys = mapper.get(chunk_number=ch, allow_missingness=args.allow_missingness)

            if isinstance(variant_indices, type(None)):
                break

            Analyser.rsid = keys
            if is_no_b4_present_in_partial_derivatives:
                genotype = np.array([])
                with Timer() as t_g:
                    genotype = merge_genotype(gen, variant_indices, mapper)
                    genotype = genotype[:, row_index["genotype"]]
                print("Time to get G {}s".format(t_g.secs))
            # TODO (low) add interaction

            a_test = np.array([])
            b_cov = np.array([])
            C = np.array([])
            a_cov = np.array([])
            b4 = np.array([])
            b_interaction = np.array([])

            if args.protocol is not None:
                if protocol.enable:
                    regression_model = protocol.regression_model()
            else:
                regression_model = None

            with Timer() as t_pd:
                if is_no_b4_present_in_partial_derivatives:
                    # If we allow for missingness in the variants. We should consider each of the variants separately
                    # In some way.
                    # In the get_expanded method we
                    if args.allow_missingness:
                        a_test, a_cov_stack, b_cov_stack, c_full_stack = meta_pard.get_expanded(
                            variant_indices=variant_indices, regression_model=regression_model,
                            random_effect_intercept=args.effect_intercept)
                    else:
                        a_test, b_cov, C, a_cov = meta_pard.get(variant_indices=variant_indices,
                                                                regression_model=regression_model,
                                                                random_effect_intercept=args.effect_intercept)
                else:
                    if args.allow_missingness:
                        a_test, a_cov_stack, b_cov_stack, c_full_stack, b4 = meta_pard.get_expanded(
                            variant_indices=variant_indices, regression_model=regression_model,
                            random_effect_intercept=args.effect_intercept)
                    else:
                        a_test, b_cov, C, a_cov, b4 = meta_pard.get(variant_indices=variant_indices, B4=True,
                                                                    regression_model=regression_model,
                                                                    random_effect_intercept=args.effect_intercept)

            print("Time to get PD {}s".format(t_pd.secs))

            MAF = meta_pard.maf_pard(SNPs_index=variant_indices)

            sample_size = meta_pard.get_n_id()

            # Check if missingness is allowed.
            if args.allow_missingness:
                # If missingness is defined, define the sample size by a matrix of dimensions nVariants * nPhenotypes,
                # with every combination portraying the number of samples available to test that combination.
                sample_size = expand_sample_size_matrix(
                    sample_size=meta_pard.get_n_per_study(),
                    variant_indices=variant_indices)

            if args.maf != 0:
                filter = (MAF > args.maf) & (MAF < 1 - args.maf) & (MAF != 0.5)
                # This probably handles the large a_test matrices as well.
                a_test = a_test[filter, :]
                Analyser.maf = MAF[filter]
                if is_no_b4_present_in_partial_derivatives:
                    genotype = genotype[filter, :]
                else:
                    b4 = b4[filter, :]
                Analyser.rsid = Analyser.rsid[filter]
                if a_test.shape[0] == 0:
                    print('NO SNPs > MAF')
                    continue
            else:
                Analyser.maf = MAF

            # Use new A_inverse that supports variable number of non-constant A parts
            # a_inv = A_inverse(a_cov, a_test)
            if args.allow_missingness:
                a_inv = a_inverse_extended_allow_missingness(a_test, a_cov_stack, variant_indices)
            else:
                a_inv = get_a_inverse_extended(a_cov, a_test)

            number_of_variable_terms = a_test.shape[2]
            number_of_constant_terms = a_inv.shape[1] - number_of_variable_terms
            print('There are {} subjects in study.'.format(meta_pard.get_n_id()))

            DF = (sample_size - a_inv.shape[1])

            if is_no_b4_present_in_partial_derivatives:

                while True:
                    # Get the next phenotype chunk.
                    phenotype = np.array([])

                    with Timer() as t_ph:
                        phenotype, phen_names = meta_phen.get()
                    print("Time to get PH {}s".format(t_ph.secs))

                    # If the phenotype type is None, the loop is done...
                    if isinstance(phenotype, type(None)):
                        # Reset the processed phenotypes when the loop is done.
                        # With the next chunk of SNPs we need to do these again
                        meta_phen.processed = 0
                        # The encoded interactions are also processed at the
                        # same rate as the phenotypes.
                        # Reset the number of processed values for this as well.
                        if encoded_interactions:
                            encoded_interactions.processed = 0
                        print('All phenotypes processed!')
                        break
                    print ("Merged phenotype shape {}".format(phenotype.shape))
                    # TODO (middle) select phen from protocol
                    phenotype = phenotype[row_index["phenotype"], :]
                    print ("Selected phenotype shape {}".format(phenotype.shape))
                    keys = meta_pard.phen_mapper.dic.keys()
                    phen_ind_dic = {k: i for i, k in enumerate(keys)}
                    phen_ind = np.array([phen_ind_dic.get(i, -1) for i in phen_names])
                    if np.sum(phen_ind == -1) == len(phen_ind):
                        print('There is no common ids in phenotype files and PD data!')
                        break
                    else:
                        print('There are {} common ids in phenotype files and PD data!'.format(np.sum(phen_ind != -1)))
                    C_test = C[phen_ind]
                    b_cov_test = b_cov[:, phen_ind]

                    b4 = B4(phenotype, genotype)
                    b_variable = b4[np.newaxis, ...]

                    # If the encoded interactions are supplied add these to the b_variable part.
                    if encoded_interactions:
                        interaction_phenotype_values, phen_names = encoded_interactions.get()
                        interaction_phenotype_values = interaction_phenotype_values[row_index["interaction"], :]
                        interaction_values = B4(interaction_phenotype_values, genotype)
                        b_variable = np.append(b_variable, [interaction_values], axis=0)

                    print ("B4 shape is {}".format(b4.shape))
                    t_stat, SE = hase_supporting_interactions(b_variable, a_inv, b_cov_test, C_test,
                                                              number_of_constant_terms, DF)

                    if mapper.cluster == 'y':
                        Analyser.cluster = True
                        Analyser.chunk = ch
                        Analyser.node = mapper.node[1]
                    Analyser.t_stat = t_stat
                    Analyser.standard_error = SE
                    Analyser.threshold = args.thr
                    Analyser.out = args.out
                    Analyser.save_result(phen_names[(phen_ind != -1)])

                    t_stat = None
                    Analyser.t_stat = None
                    gc.collect()

            else:
                t_stat, SE = hase_supporting_interactions(b4, a_inv, b_cov, C, number_of_constant_terms, DF)

                Analyser.t_stat = t_stat
                Analyser.standard_error = SE
                Analyser.threshold = args.thr
                Analyser.out = args.out
                Analyser.save_result(np.array(meta_pard.phen_mapper.dic.keys()))

                t_stat = None
                Analyser.t_stat = None
                gc.collect()

    ################################### CLASSIC META STAGE ##############################

    elif args.mode == 'meta-classic':

        # ARG_CHECKER.check(args,mode='meta-stage')

        ##### Init data readers #####
        if args.derivatives is None:
            raise ValueError('For meta-classic analysis partial derivatives data are required!')

        mapper = load_mapper(
            mapper_chunk_size=MAPPER_CHUNK_SIZE, study_name=args.study_name, ref_name=args.ref_name,
            mapper_folder=args.mapper, encoded=args.encoded, cluster=args.cluster, node=args.node,
            snp_id_inc=args.snp_id_inc, snp_id_exc=args.snp_id_exc)

        # It appears to me that the pard list contains
        # partial derivatives matrices for every study
        partial_derivatives_folders = []

        with Timer() as t:
            for i, j in enumerate(args.derivatives):
                partial_derivatives_folders.append(Reader('partial'))
                partial_derivatives_folders[i].start(j, study_name=args.study_name[i])
                partial_derivatives_folders[i].folder.load()

        print("Time used to load partial derivatives is {}s".format(t.secs))

        # Create a list containing booleans representing the presence
        # of the b4 data in the partial derivatives folders
        b4_presence_per_study = [False if isinstance(i.folder._data.b4, type(None)) else True for i in
                                 partial_derivatives_folders]

        # The if statement roughly translates to:
        # if one study contains b4 data in the partial derivatives folder,
        # all should have b4 data, otherwise an error is raised.
        if np.sum(b4_presence_per_study) != len(partial_derivatives_folders) and np.sum(b4_presence_per_study) != 0:
            raise ValueError('All studies should have b4 data for partial derivatives!')

        if args.protocol is not None:
            raise NotImplementedError("Protocols not implemented in classic meta analysis")
            protocol = Protocol(args.protocol)
        else:
            protocol = None

        meta_pard = MetaParData(partial_derivatives_folders, args.study_name,
                                protocol=protocol,
                                allow_missingness=args.allow_missingness)

        is_no_b4_present_in_partial_derivatives = np.sum(b4_presence_per_study) == 0
        if not is_no_b4_present_in_partial_derivatives:
            raise NotImplementedError("Classic meta analysis with b4 present not implemented")

        phen = []

        with Timer() as t:
            for i, j in enumerate(args.phenotype):
                phen.append(Reader('phenotype'))
                phen[i].start(j)
        print("Time to set pheno {} s".format(t.secs))
        meta_phen = MetaPhenotype(phen,
                                  include=args.ph_id_inc,
                                  exclude=args.ph_id_exc,
                                  allow_missingness=args.allow_missingness)

        gen = []
        with Timer() as t:
            for i, j in enumerate(args.genotype):
                gen.append(Reader('genotype'))
                gen[i].start(j, hdf5=args.hdf5, study_name=args.study_name[i], ID=False)
        print("Time to set gen {}s".format(t.secs))

        # Create a dictionary with the datasets for obtaining
        # the indices of shared identifiers
        datasets = {"phenotype": tuple(i.folder._data for i in phen),
                    "genotype": tuple(i.folder._data for i in gen),
                    "partial_derivatives": tuple(i.folder._data.metadata for i in partial_derivatives_folders)}

        # Get common ids
        row_index, intersecting_identifiers = get_intersecting_individual_indices(datasets)

        # Do something for all covariates
        if row_index["partial_derivatives"].shape[0] != np.sum(
                [i.folder._data.metadata['id'].shape[0] for i in partial_derivatives_folders]):
            raise ValueError(
                'Partial Derivatives covariates have different number of subjects {} than genotype and phenotype {}'.format(
                    row_index["covariates"].shape[0],
                    np.sum([i.folder._data.metadata['id'].shape[0] for i in partial_derivatives_folders])))

        covariate_indices = None
        if args.covariate_indices:

            covariate_indices_df = pd.read_csv(
                args.covariate_indices,
                sep="\t", header=None, dtype=str, index_col=0)

            covariate_indices = covariate_indices_df.apply(
                lambda col: col.astype(np.float), axis=1).T.to_dict('list')

            for key, indices in covariate_indices.items():
                if np.min(indices) < 1:
                    raise ValueError(
                        'Lowest covariate index below 1 for cohort {}. '
                        'Covariate indices must start at 1 or above.'.format(
                            key
                        ))

                covariate_indices[key] = np.array(indices).astype(int)

        classic_meta_analyser = ClassicMetaAnalyser(
            meta_phen, meta_pard, intersecting_identifiers, row_index,
            args.study_name, args.out,
            covariate_indices=covariate_indices, maf_threshold=args.maf,
            t_statistic_threshold=args.thr)

        # Start looping over all genotype chunks
        # We use while true, since implementing an iterator is apparently
        # too much work...
        while True:
            # We start with an empty chunk...
            ch = None

            # Then we check whether or not we are working in an environment
            # wherein analysis is split over multiple nodes.

            # This determines how the current chunk should be obtained.
            if mapper.cluster == 'n':
                variant_indices, keys = mapper.get(
                    allow_missingness=args.allow_missingness)
            elif mapper.cluster == 'y':
                ch = mapper.chunk_pop()
                if ch is None:
                    break
                variant_indices, keys = mapper.get(
                    chunk_number=ch, allow_missingness=args.allow_missingness)

            # Now we can check if we have looped through all genotype chunks...
            if isinstance(variant_indices, type(None)):
                break

            # And we can get the genotypes for this chunk...
            if is_no_b4_present_in_partial_derivatives:
                with Timer() as t_g:
                    genotype = merge_genotype(gen, variant_indices, mapper)
                print("Time to get G {}s".format(t_g.secs))

            # Now analyse the genotype chunk
            if mapper.cluster == 'y':
                classic_meta_analyser.analyse_genotype_chunk(
                    genotype, keys, variant_indices,
                    chunk=ch, node=mapper.node[1])
            else:
                classic_meta_analyser.analyse_genotype_chunk(
                    genotype, keys, variant_indices)


    ################################### TO DO EVERYTHING IN ONE GO ##############################

    elif args.mode == 'regression':

        # ARG_CHECKER.check(args,mode='regression')

        print ('START regression mode...')
        if args.mapper is None:
            os.environ['MAPPER'] = "False"  # need to set it here, before genotype Reader start

        phen = Reader('phenotype')
        phen.start(args.phenotype[0])
        phen.permutation = args.permute_ph

        cov = Reader('covariates')
        cov.start(args.covariates)

        interaction = None
        if args.interaction is not None:
            interaction = Reader('interaction')
            interaction.start(args.interaction)

        if (cov.folder.n_files > 1 and cov.folder.format != '.npy') or (
                cov.folder.n_files > 2 and cov.folder.format == '.npy'):  # TODO (middle) test
            raise ValueError('In covariates folder should be only one file!')

        gen = []
        for i, j in enumerate(args.genotype):
            gen.append(Reader('genotype'))
            gen[i].start(j, hdf5=args.hdf5, study_name=args.study_name[i], ID=False)

        if args.mapper is not None:
            mapper = Mapper()
            mapper.chunk_size = MAPPER_CHUNK_SIZE
            mapper.genotype_names = args.study_name
            mapper.reference_name = args.ref_name
            if args.snp_id_inc[0] is not None:
                mapper.include = pd.DataFrame.from_csv(args.snp_id_inc[0], index_col=None)
                print('Include:')
                print(mapper.include.head())
                if 'ID' not in mapper.include.columns and (
                        'CHR' not in mapper.include.columns or 'bp' not in mapper.include.columns):
                    raise ValueError('{} table does not have ID or CHR,bp columns'.format(args.snp_id_inc[0]))
            if args.snp_id_exc[0] is not None:
                mapper.exclude = pd.DataFrame.from_csv(args.snp_id_exc[0], index_col=None)
                print('Exclude:')
                print(mapper.exclude.head())
                if 'ID' not in mapper.exclude.columns and (
                        'CHR' not in mapper.exclude.columns or 'bp' not in mapper.exclude.columns):
                    raise ValueError('{} table does not have ID or CHR,bp columns'.format(args.snp_id_exc[0]))
            mapper.load(args.mapper)
            mapper.load_flip(args.mapper)
            mapper.cluster = args.cluster
            mapper.node = args.node
        else:
            if len(args.genotype) == 1:
                mapper = Mapper()
                mapper.chunk_size = MAPPER_CHUNK_SIZE
                mapper.genotype_names = args.study_name
                mapper.reference_name = args.ref_name
                mapper.cluster = args.cluster
                mapper.node = args.node
                mapper.n_study = 1
                mapper.n_keys = gen[0].folder._data.names.shape[0]
                mapper.keys = np.array(gen[0].folder._data.names.tolist())
                mapper.values = np.array(range(mapper.n_keys)).reshape(-1, 1)
                mapper.flip[args.study_name[0]] = np.array([1] * mapper.n_keys)
                if args.snp_id_exc is not None or args.snp_id_inc is not None:
                    raise ValueError('You can not exclude or include variants to analysis without mapper!')
            else:
                raise ValueError('You can not run regression analysis with several genotype data without mapper!')
        # mapper=None

        Analyser = HaseAnalyser()
        Analyser.threshold = args.thr
        Analyser.out = args.out


        haseregression(phen, gen, cov, mapper, Analyser, args.maf, intercept=args.intercept, interaction=interaction)

    end = time.time()

    print ('experiment finished in {} s'.format((end - start)))


def load_mapper(mapper_chunk_size, study_name, ref_name,
                mapper_folder, encoded, cluster, node,
                snp_id_inc=None, snp_id_exc=None):
    # Create Mapper object. This mapper object makes us able to synchronise across al variants
    # in each of the separate datasets
    mapper = Mapper()
    # We set the chunk size in this mapper. Each chunk is analyzed separately.
    mapper.chunk_size = mapper_chunk_size
    # Assign study and reference names
    mapper.genotype_names = study_name
    mapper.reference_name = ref_name  # Reference dataset
    # Check if we are to exlude, or explicitly include some set of variants
    if snp_id_inc is not None:  # If this is not none the argument contains table(s) of snps to include
        mapper.load_filter_include(snp_id_inc)
    if snp_id_exc is not None:  # If this is not None the argument contains table(s) of snps to exclude
        mapper.load_filter_exclude(snp_id_exc)
    # Load the mapper files for each of the datasets
    mapper.load(mapper_folder)  # Load the mapper files
    # Load the flip files for each of the datasets
    mapper.load_flip(mapper_folder, encode=encoded)  # often args.encoded is is null
    # For cluster stuff:
    mapper.cluster = cluster  # Is n by default
    mapper.node = node
    return mapper


if __name__ == '__main__':
    sys.exit(main())
