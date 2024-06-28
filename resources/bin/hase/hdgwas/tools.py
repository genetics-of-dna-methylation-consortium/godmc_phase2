from __future__ import print_function

import time
import numpy as np
import pandas as pd
import tables
import h5py
import bitarray as ba
import gc
from collections import OrderedDict
import os
from scipy import stats
from hash import *
import glob
import inspect, itertools


def timer(func):
    def f(*args, **kwargs):
        with Timer(verbose=True) as t:
            return func(*args, **kwargs)

    return f


class Timer(object):
    def __init__(self, verbose=False):
        self.verbose = verbose

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        self.msecs = self.secs * 1000  # millisecs
        if self.verbose:
            print('elapsed time: %f ms' % self.msecs)


def timing(f):
    def wrap(*args, **kwargs):
        time1 = time.time()
        ret = f(*args, **kwargs)
        time2 = time.time()
        print('%s function took %0.3f s' % (f.func_name, (time2 - time1)))
        return ret

    return wrap


def save_parameters(f):
    def wrap(*args, **kwargs):
        args_name = inspect.getargspec(f)[0]
        args_dict = dict(itertools.izip(args_name, args))
        np.save('{}.npy'.format(f.func_name), args_dict)
        ret = f(*args, **kwargs)
        return ret

    return wrap


def check_np():
    for i in ['atlas_blas_info', 'openblas_info', 'lapack_opt_info',
              'atlas_info', 'lapack_mkl_info', 'blas_mkl_info',
              'atlas_blas_info', 'mkl_info']:
        try:
            if len(getattr(np.__config__, i)) > 0:
                print('Numpy linked to {}'.format(i))
                return 0
        except:
            continue

    raise ValueError('NO BLAS/LAPACK/MKL found!'
                     'To restart with slow speed use -np False')


########################################################################
########################################################################
########################################################################

class Analyser(object):

    def __init__(self, name):
        self.name = name


class HaseAnalyser(Analyser):

    def __init__(self):
        self.name = None
        self.t_stat = None
        self.p_value = None
        self.betas = None
        self.standard_error = None
        self.DF = None
        self.n_studies = None
        self.probes_path = {}
        self.results = OrderedDict()
        self.split_size = None
        self.threshold = None
        self.result_path = None
        self.result_index = 1
        self.out = None
        self.maf = None
        self.rsid = None
        self.reference = None
        self.cluster = False
        self.permutation = False
        self.rsid_dic = {}
        self.result_folder = None
        self.result_dump_size = 3
        self.file_number = None

    def get_betas(self):
        """
        Returns a 3d array of beta values
        :return: 3d array of beta values.
        -   The 1st dimension of this array comprises each of the included
            phenotypes for the respective cohort or study.
        -   The 2nd dimension of this array does not seem to have a function.
            (could this perhaps have been meant for interaction effects?)
        -   The 3rd dimension of this array represents all of the phenotypes
            that are tested for.
        """
        if self.t_stat.shape == self.standard_error.shape:
            return self.t_stat * self.standard_error

    def summary(self):

        if isinstance(self.result_path, type(None)):
            raise ValueError('Please set result pathway to Analyser!')

        if isinstance(self.DF, type(None)):
            print ('DF is not defined. Forced to use z_score statistics!')

        if self.result_folder is None:
            self.result_folder = [i for i in glob.glob(os.path.join(self.result_path, '*.npy')) if 'RSID' not in i]

        self.results['RSID'] = np.array([])
        self.results['p_value'] = np.array([])
        self.results['t-stat'] = np.array([])
        self.results['phenotype'] = np.array([])
        self.results['SE'] = np.array([])
        self.results['MAF'] = np.array([])
        self.results['BETA'] = np.array([])

        files = []

        if self.file_number is not None:
            try:
                d = np.load(os.path.join(self.result_path, '{}result.npy'.format(self.file_number)), allow_pickle=True).item()
                print(self.file_number)
            except:
                print("Can't read {}".format(i))

            p_value = self.get_p_value(d['t-stat'], df=self.DF)
            self.results['t-stat'] = np.append(self.results['t-stat'], d['t-stat'].flatten())
            self.results['SE'] = np.append(self.results['SE'], d['SE'].flatten())
            self.results['RSID'] = np.append(self.results['RSID'], d['index'])
            self.results['phenotype'] = np.append(self.results['phenotype'], d['phenotype'])
            self.results['MAF'] = np.append(self.results['MAF'], d['MAF'])
            self.results['p_value'] = np.append(self.results['p_value'], p_value.flatten())
            self.results['BETA'] = np.append(self.results['BETA'], d['t-stat'].flatten() * d['SE'].flatten())

        else:

            for i in range(self.result_dump_size):
                try:
                    files.append(self.result_folder.pop())
                except:
                    break
            if len(files) != 0:
                for i in files:
                    print(i)
                    try:
                        d = np.load(i, allow_pickle=True).item()
                    except:
                        print("Can't read {}".format(i))
                        continue
                    p_value = self.get_p_value(d['t-stat'], df=self.DF)
                    self.results['t-stat'] = np.append(self.results['t-stat'], d['t-stat'].flatten())
                    self.results['SE'] = np.append(self.results['SE'], d['SE'].flatten())
                    self.results['RSID'] = np.append(self.results['RSID'], d['index'])
                    self.results['phenotype'] = np.append(self.results['phenotype'], d['phenotype'])
                    self.results['MAF'] = np.append(self.results['MAF'], d['MAF'])
                    self.results['p_value'] = np.append(self.results['p_value'], p_value.flatten())
                    self.results['BETA'] = np.append(self.results['BETA'], d['t-stat'].flatten() * d['SE'].flatten())
            else:
                self.results = None

    def get_p_value(self, t_stat, df=None):
        if df is None:
            return stats.norm.sf(np.abs(t_stat)) * 2
        else:
            return stats.t.sf(np.abs(t_stat), df) * 2

    def save_result(self, phen_names):

        t_threshold = self.threshold
        save_path = self.out

        mask = ([], [], [])

        if t_threshold:
            mask = np.where(np.abs(self.t_stat) > t_threshold)
        else:
            mask = np.where(np.abs(self.t_stat) > 0)

        if (len(mask[0]) != 0):
            print ('Saving results to {}'.format(save_path))
            t_save = self.t_stat[mask[0], mask[1], mask[2]]
            se = self.standard_error[mask[0], mask[1], mask[2]]
            result = {'phenotype': phen_names[mask[2]], 't-stat': t_save, 'index': self.rsid[mask[0]], 'SE': se,
                      'MAF': self.maf[mask[0]]}

            if not self.cluster:
                np.save(os.path.join(save_path, str(self.result_index) + 'result.npy'), result)
            else:
                if self.permutation:
                    result = self.t_stat
                self.rsid_dic[str(self.chunk[0]) + '_' + str(self.chunk[1])] = self.rsid
                np.save(os.path.join(save_path, 'node' + str(self.node) + '_' + str(self.chunk[0]) + '_' + str(
                    self.chunk[1]) + '_' + str(self.result_index) + 'result.npy'), result)
        self.result_index += 1


########################################################################
########################################################################
########################################################################
class Log(object):  # TODO (low) write this class

    def __init__(self):
        pass


########################################################################
########################################################################
########################################################################
class Checker(object):  # TODO (mid) finish or remove this class

    def __init__(self):
        self.mode = None
        self.name = None
        self.args = None

    def check(self, args, mode=None):

        if mode == 'converting':
            return self.converting(args)
        elif mode == 'encoding':
            self.encoding(args)
        elif mode == 'single-meta':
            self.single_meta(args)
        elif mode == 'meta-stage':
            self.meta_stage(args)
        elif mode == 'regression':
            self.regression(args)

    def system_check(self, args):
        pass

    def converting(self, args):

        if len(args.study_name) > 1:
            raise ValueError('Should be only one study name!')

        if len(args.genotype) > 1:
            raise ValueError('Should be only one directory for genotype!')

        if isinstance(args.study_name, type(None)):
            raise ValueError('study name is not defined')

        args.genotype = args.genotype[0]
        if not os.path.isdir(args.genotype):
            raise ValueError('{} is not a directory!'.format(args.genotype))

    def encoding(self, args):

        if len(args.study_name) > 1:
            raise ValueError('Should be only one study name!')

        if len(args.genotype) > 1:
            raise ValueError('Should be only one directory for genotype!')

        if isinstance(args.study_name, type(None)):
            raise ValueError('study name is not defined')

        args.genotype = args.genotype[0]
        args.study_name = args.study_name[0]

        if not os.path.isdir(args.genotype) or not os.path.isdir(args.phenotype) or not os.path.isdir(args.out):
            raise ValueError('{},{} or {} is not a directory'.format(args.genotype, args.phenotype, args.out))

        if not os.path.isdir(os.path.join(args.out, 'encode_genotype')) or not os.path.isdir(
                os.path.join(args.out, 'encode_phenotype')):
            raise ValueError('{} or {} not exist it is up to you to create these directory'.format(
                (os.path.join(args.out, 'encode_genotype')),
                (os.path.join(args.out, 'encode_phenotype'))
            )
            )

        id = np.array(pd.read_hdf(os.path.join(args.genotype, 'individuals', args.study_name + '.h5'),
                                  'individuals').individual.tolist(), dtype=np.int)

    def single_meta(self, args):

        if len(args.study_name) > 1:
            raise ValueError('Should be only one study name!')

        if len(args.genotype) > 1:
            raise ValueError('Should be only one directory for genotype!')

        if isinstance(args.study_name, type(None)):
            raise ValueError('study name is not defined')

        args.genotype = args.genotype[0]
        args.study_name = args.study_name[0]

        if not os.path.isdir(args.genotype) or not os.path.isdir(args.phenotype) or not os.path.isdir(args.covariates):
            raise ValueError(
                '{},{} or {} is nor a directory'.format(os.path.isdir(args.genotype), os.path.isdir(args.phenotype),
                                                        os.path.isdir(args.covariates)))

    def meta_stage(self, args):

        if isinstance(args.mapper, type(None)):
            raise ValueError(
                'You should define mapper folder --mapper, if it does not exist you should first run Mapper script')

        if isinstance(args.mapper_name, type(None)):
            raise ValueError(
                'You should define mapper name --mapper_name, if it does not exist you should first run Mapper script')

        if len(args.genotype) <= 1:
            raise ValueError('Should be only one directory for genotype!')

        g = [os.path.isdir(i) for i in args.genotype]
        d = [os.path.isdir(i) for i in args.derivatives]

        if np.sum(g) != len(args.genotype) or not os.path.isdir(args.phenotype) or np.sum(d) != len(args.derivatives):
            raise ValueError(
                '{},{} or {} is nor a directory'.format(os.path.isdir(args.genotype), os.path.isdir(args.phenotype),
                                                        os.path.isdir(args.derivatives)))

    def regression(self, args):

        g = [os.path.isdir(i) for i in args.genotype]

        if np.sum(g) != len(args.genotype) or not os.path.isdir(args.phenotype[0]) or not os.path.isdir(
                args.covariates):
            raise ValueError(
                '{},{} or {} is not a directory'.format(os.path.isdir(args.genotype), os.path.isdir(args.phenotype),
                                                        os.path.isdir(args.derivatives)))


########################################################################
########################################################################
########################################################################
class Mapper(object):

    def __init__(self):

        self.is_filtered = False
        self.name = None
        self.genotype_names = []
        self.dic = OrderedDict()
        self.n_study = 0
        self.values = np.array([])  # array (N_study, N_keys), where N_keys = number of ID in reference table
        self.keys = None
        self.n_keys = None
        self.processed = 0
        self._limit = None
        self.reference = None
        self.include = None
        self.exclude = None
        self.include_ind = np.array([], dtype=np.int64)
        self.exclude_ind = np.array([], dtype=np.int64)
        # self.hash=HashTable()
        self.cluster = None
        self.node = None
        self.chunk_pool = None
        self.chunk_size = None
        self.column_names = []
        self.flip = {}  # keys - study names; values - array, with length equal to number of probes
        self.probes = None
        self.encoded = {}

    def chunk_pop(self):
        if self.chunk_pool is None and self.cluster is None:
            raise ValueError('Cluster settings are not defined!')

        if self.chunk_pool is None:
            self.chunk_pool = []
            if self.node is None:
                raise ValueError('cluster setting not defined! (number of nodes / this node number) ')
            self.node = [int(i) for i in self.node]
            N = self.n_keys / self.node[0]
            n = N / self.chunk_size
            if N < 1:
                raise ValueError('Too many nodes! Change chunk size in mapper')
            elif N == 1:
                self.chunk_pool = [[self.node[1] - 1, self.node[1]]]
            else:
                if n == 0:
                    self.chunk_pool = [[N * (self.node[1] - 1), N * (self.node[1] - 1) + N]]
                else:
                    self.chunk_pool = []
                    for i in range(n):
                        self.chunk_pool.append([N * (self.node[1] - 1) + self.chunk_size * i,
                                                N * (self.node[1] - 1) + self.chunk_size * (i + 1)])
                    if (N - n * self.chunk_size) != 0:
                        self.chunk_pool.append([N * (self.node[1] - 1) + self.chunk_size * (n),
                                                N * (self.node[1] - 1) + self.chunk_size * (n) + (
                                                        N - n * self.chunk_size)])

            if self.node[1] == self.node[0]:
                if (self.n_keys - N * self.node[0]) != 0:
                    self.chunk_pool.append([self.chunk_pool[-1][1], self.n_keys])
            self.chunk_pool = self.chunk_pool[::-1]

        if len(self.chunk_pool) != 0:
            ch = self.chunk_pool.pop()
            return ch
        else:
            return None

    # @timing
    def fill(self, keys, name, repeats=False, reference=None):  # TODO (middle) remove
        """
		Fill the mapper
		:param keys: The labels for phenotypes/covariates.single that are passed
		:param name: Study name
		:param repeats: Dunno, False is passed afaik
		:param reference: Dunno, False is passed afaik
		"""
        # Assign some stuff...
        # Assign the given reference to the field, this appears to be False most of the time
        self.reference = reference
        self.reference_name = name
        # Add the name of the column
        self.column_names.append(name)

        # Check the number of labels,
        # and the number of unique labels.
        number_of_keys = len(keys)
        unique_keys = np.unique(keys)
        number_of_unique_keys = len(unique_keys)

        if not repeats:
            # The labels should be unique,
            # to check this the lengths of the keys and the unique keys should be equal.
            if number_of_keys != number_of_unique_keys:
                raise ValueError('length of keys {} does not =  {} number of uniq values!!!'.format(
                    number_of_keys, number_of_unique_keys))
        else:
            # When repeats is true, ignore the fact that stuff is not unique
            number_of_keys = number_of_unique_keys
            keys = unique_keys
        # Set n_keys to the number of labels
        self.n_keys = number_of_keys
        if isinstance(self._limit, type(None)):
            # If _limit is not yet set,
            # loop through the labels/keys to fill a dictionary with
            # the labels/keys as keys, and the indices of those labels
            # as values so that the keys and values correspond to each other.
            for i, j in enumerate(keys):
                self.dic[j] = [i]
        else:
            # Do the same as in the previous if clause, untill the index is
            # equal to _limit
            for i, j in enumerate(keys):
                if i == self._limit:
                    break
                self.dic[j] = [i]

    # @timing
    def push(self, new_keys, name=None, new_id=True):  # TODO (middle) remove
        # Add stuff to the mapper
        if not self.reference and len(self.dic) == 0:
            raise ValueError('You should fill mapper first with ref panel or your own rsids!')
        self.n_study += 1
        # The name should not be the same as the names of previous study
        if name is not None:
            try:
                # Exception should be thrown.
                self.genotype_names.index(name)
                raise ValueError('Trying to push the same study to mapper')
            except:
                self.genotype_names.append(name)
                self.column_names.append(name)

        # If the limit is not set
        if isinstance(self._limit, type(None)):
            for i, j in enumerate(new_keys):
                # Fill the dictionary in the same order as the ref
                if self.dic.get(j):
                    self.dic[j] = self.dic[j] + [i]
                else:
                    if new_id:
                        print ('WARNING! You are pushing ID {} which is not present in reference panel!'.format(j))
                        self.dic[j] = [-1] * (self.n_study + 1)
                        self.dic[j][self.n_study] = i
                    else:
                        continue
        else:
            for i, j in enumerate(new_keys):
                if i == self._limit:
                    break
                if self.dic.get(j):
                    self.dic[j] = self.dic[j] + [i]
                else:
                    self.dic[j] = [-1] * (self.n_study + 1)
                    self.dic[j][self.n_study] = i

        for k in self.dic:
            if len(self.dic[k]) < self.n_study + 1:
                self.dic[k] = self.dic[k] + [-1]

    # @timing
    def load_flip(self, folder, encode=None):

        if folder is None:
            raise ValueError('Mapper is not defined!')
        if not os.path.isdir(folder):
            raise ValueError('{} is not a folder'.format(folder))

        if self.reference_name is None:
            raise ValueError('Reference name for mapper is not defined!')

        #This is edited by me-Olalaken
        #if len(glob.glob(folder + 'flip_*')) == 0:
        if len(glob.glob(os.path.join(folder, 'flip_*'))) == 0:
            raise ValueError('There is no flip mapper data in folder {}'.format(folder))

        for j, i in enumerate(self.genotype_names):
            if encode is not None:
                self.encoded[i] = encode[j]
                if encode[j] == 0:
                    self.flip[i] = np.load(os.path.join(folder, 'flip_' + self.reference_name + '_' + i + '.npy'))
                elif encode[j] == 1:
                    flip = np.load(os.path.join(folder, 'flip_' + self.reference_name + '_' + i + '.npy'))
                    self.flip[i] = np.ones(flip.shape[0])
            else:
                self.flip[i] = np.load(os.path.join(folder, 'flip_' + self.reference_name + '_' + i + '.npy'))

    # @timing
    def load(self, folder):
        if folder is None:
            raise ValueError('Mapper is not defined!')
        if not os.path.isdir(folder):
            raise ValueError('{} is not a folder'.format(folder))

        if self.reference_name is None:
            raise ValueError('Reference name for mapper is not defined!')

        #This is edited by me-Olalekan
        #if len(glob.glob(folder + 'keys_*')) == 0:
        if len(glob.glob(os.path.join(folder, 'keys_*'))) == 0:
            raise ValueError('There is no mapper data in folder {}'.format(folder))


        #This is edited by me-Olalekan
        #keys = glob.glob(folder + 'keys_*')
        keys = glob.glob(os.path.join(folder, 'keys_*'))
        if len(keys) > 1:
            raise ValueError('There are more than one reference keys in folder {}'.format(folder))

        self.keys = np.load(os.path.join(keys[0]), allow_pickle=True)  # TODO (middle) not safety to load only one file
        self.n_keys = self.keys.shape[0]

        #This is edited by me-Olalekan
        #values = glob.glob(folder + 'values_*')
        values = glob.glob(os.path.join(folder, 'values_*'))
        if len(values) == 0:
            raise ValueError('There is no mapper data in folder {}'.format(folder))

        if self.genotype_names is None:
            raise ValueError('Genotype names are not defined!')

        self.values = np.zeros((self.n_keys, len(self.genotype_names)))
        # Mapper files should always maintain the same order as the order of study names, so they stay aligned
        # with all other data structures.
        for j, i in enumerate(self.genotype_names):
            self.values[:, j] = np.load(os.path.join(folder, 'values_' + self.reference_name + '_' + i + '.npy'))

        if self.n_keys != self.values.shape[0]:
            raise ValueError(
                'Number of indexes {} in mapper values is different from keys length {}!'.format(self.values.shape[0],
                                                                                                 self.n_keys))

        self.n_study = self.values.shape[1]
        if self.n_keys < self.n_study:
            print ('WARNING!!! Normally should be more test values that studies! Check your saved data for Mapper')

        print ('You loaded values for {} studies and {} test values'.format(self.n_study, self.n_keys))

        if self.include is not None or self.exclude is not None or True:
            self.reference = Reference(self.reference_name)
            self.reference.load_index()

    # @timing
    def get_old(self, chunk_number=None):

        if isinstance(self.keys, type(None)) and isinstance(self.values, type(None)):
            raise ValueError('mapper data is not loaded!')

        if isinstance(self.chunk_size, type(None)):
            raise ValueError('chunk_size and should be defined')

        if self.processed == self.n_keys:
            return None, None

        if chunk_number is not None:
            start = chunk_number[0]
            finish = chunk_number[1]
        else:
            start = self.processed
            finish = self.processed \
                     + self.chunk_size if (self.processed + self.chunk_size) < self.n_keys else self.n_keys
            self.processed = finish

        with Timer() as t_q:
            if self.include is not None:
                if len(self.include_ind) == 0:
                    if 'ID' in self.include.columns:
                        self.include_ind = np.in1d(self.keys, self.include.ID)
                        self.include_ind = np.where(self.include_ind == True)[0]
                    # self.query_include=self.reference.index.select('reference',where='ID=self.include.ID')
                    # self.include_ind=self.query_include.index
                    elif 'CHR' in self.include.columns and 'bp' in self.include.columns:
                        self.query_include = self.reference.index.select('reference',
                                                                         where='CHR=self.include.CHR & bp=self.include.bp')
                        self.include_ind = self.query_include.index

            if self.exclude is not None:
                if len(self.exclude_ind) == 0:
                    if 'ID' in self.exclude.columns:
                        self.exclude_ind = np.in1d(self.keys, self.exclude.ID)
                        self.exclude_ind = np.where(self.exclude_ind == True)[0]
                    # self.query_exclude=self.reference.index.select('reference',where='ID=self.exclude.ID')
                    # self.exclude_ind=self.query_exclude.index
                    elif 'CHR' in self.exclude.columns and 'bp' in self.exclude.columns:
                        self.query_exclude = self.reference.index.select('reference',
                                                                         where='CHR=self.exclude.CHR & bp=self.exclude.bp')
                        self.exclude_ind = self.query_exclude.index

        print("Time to select SNPs {}s".format(t_q.secs))

        if self.exclude is not None or self.include is not None:
            self.include_ind = np.setxor1d(self.include_ind, self.exclude_ind)
            if len(self.include_ind) == 0:
                print('None of included ID found in genotype data!')
                return None, None
            else:
                ind = np.intersect1d(np.arange(start, finish), self.include_ind)
        else:
            ind = np.arange(start, finish)
            self.include_ind = np.arange(len(self.values))

        ind = ind.astype('int')
        indexes = self.values[ind, :]
        r = (indexes == -1).any(axis=1)
        ind = ind[~r]
        indexes = indexes[~r]
        keys = self.keys[ind]

        while len(ind) < self.chunk_size:
            if self.processed == self.n_keys:
                break
            if chunk_number is None:
                start = self.processed
                finish = self.processed + self.chunk_size if (self.processed + self.chunk_size) \
                                                             < self.n_keys else self.n_keys
                self.processed = finish
            else:
                ch = self.chunk_pop()
                if ch is not None:
                    start = ch[0]
                    finish = ch[1]
                else:
                    break
            if len(np.append(ind, np.intersect1d(np.arange(start, finish), self.include_ind))) > 1.5 * self.chunk_size:

                if chunk_number is None:
                    self.processed = self.processed - len(np.arange(start, finish))
                    break
                else:
                    print ('Added back chunk {}'.format(ch))
                    self.chunk_pool.append(ch)
                    break
            else:
                ind = np.append(ind, np.intersect1d(np.arange(start, finish), self.include_ind))
                ind = ind.astype('int')
                indexes = self.values[ind, :]
                r = (indexes == -1).any(axis=1)
                ind = ind[~r]
                indexes = indexes[~r]
                keys = self.keys[ind]

        if len(ind) == 0:
            return None, None

        return [indexes[:, i].astype(np.int64) for i in range(self.n_study)], keys

    # @timing
    def get(self, chunk_number=None, allow_missingness=False):

        if isinstance(self.keys, type(None)) and isinstance(self.values, type(None)):
            raise ValueError('mapper data is not loaded!')

        if isinstance(self.chunk_size, type(None)):
            raise ValueError('chunk_size and should be defined')

        if self.processed == self.n_keys:
            return None, None

        if chunk_number is not None:
            start = chunk_number[0]
            finish = chunk_number[1]
        else:
            start = self.processed
            finish = self.processed \
                     + self.chunk_size if (self.processed + self.chunk_size) < self.n_keys else self.n_keys
            self.processed = finish

        with Timer() as t_q:

            if not self.is_filtered:

                # Now, depending on the presence of a grouping by dataset,
                # either set all unwanted variants in 'values' to -1, or
                # only do this for values unwanted per group.
                #

                # First make a dictionary with keys and indices
                key_index_dictionary = dict(zip(self.keys, range(0, len(self.keys))))
                print("Finished making key:index dictionary of size: {}".format(len(key_index_dictionary)))

                if self.include is not None:
                    if 'ID' in self.include.columns:
                        # Now select the indices for those keys that should be included.
                        print("Subsetting mapper to include specified variants. (n = {})".format(self.include.shape[0]))
                        self.include['indices'] = self.include.ID.apply(lambda x: key_index_dictionary.get(x, (-1)))
                        print("Retrieved all indices of variants that should be included!")
                        self.include = self.include[self.include['indices'] != -1]
                        print("Number of variants that remain after excluding variants not in mapper: {}".format(self.include.shape[0]))
                    elif 'CHR' in self.include.columns and 'bp' in self.include.columns:
                        raise NotImplementedError()
                        # TODO: check first if this works as expected

                        # self.query_include = self.reference.index.select(
                        #     'reference', where='CHR=self.include.CHR & bp=self.include.bp')
                        # self.include['indices'] = self.query_include.index
                    if "cohort" not in self.include.index.names:
                        self.include = pd.concat(
                            [self.include]*self.values.shape[1],
                            keys=np.arange(self.values.shape[1]),
                            names=["cohort", "row"])

                    new_values = np.full_like(self.values, -1)
                    new_values[
                        self.include.indices,
                        self.include.index.get_level_values("cohort").values] = (
                        self.values[
                            self.include.indices,
                            self.include.index.get_level_values("cohort").values])

                    self.values = new_values


                if self.exclude is not None:
                    if 'ID' in self.exclude.columns:
                        # Now select the indices for those keys that should be included.
                        self.exclude['indices'] = self.exclude.ID.apply(lambda x: key_index_dictionary.get(x, (-1)))
                        self.exclude = self.exclude[self.exclude['indices'] != -1]
                    elif 'CHR' in self.exclude.columns and 'bp' in self.exclude.columns:
                        raise NotImplementedError()
                        # TODO: check first if this works as expected

                        # self.query_exclude = self.reference.index.select(
                        #     'reference', where='CHR=self.include.CHR & bp=self.include.bp')
                        # self.exclude['indices'] = self.query_exclude.index
                    if "cohort" not in self.exclude.index.names:
                        self.exclude = pd.concat(
                            [self.exclude]*self.values.shape[1],
                            keys=np.arange(self.values.shape[1]),
                            names=["cohort", "row"])

                    self.values[self.exclude.indices, self.exclude.index.get_level_values("cohort").values] = -1
                self.is_filtered = True

        print("Time to select SNPs {}s".format(t_q.secs))

        if allow_missingness:
            # Get the indices for those values that are present in ANY studies
            # Values that are greater than -1 are present
            all_indices_for_values_with_all_studies_available = np.where(
                (self.values[:, :] > -1).any(axis=1))[0]
        else:
            # Get the indices for those values that are present in ALL studies
            # Values that are greater than -1 are present
            all_indices_for_values_with_all_studies_available = np.where(
                (self.values[:, :] > -1).all(axis=1))[0]

        all_indices_for_values_with_all_studies_available = all_indices_for_values_with_all_studies_available[
            all_indices_for_values_with_all_studies_available >= start]

        # Loop over chunks, popping them, while the current chunk does not contain more than self.chunk_size values.

        while np.sum(all_indices_for_values_with_all_studies_available < finish) < self.chunk_size:
            if self.processed == self.n_keys:
                break

            # Get the new finish
            if chunk_number is None:
                start = self.processed
                finish = self.processed + self.chunk_size if (self.processed + self.chunk_size) \
                                                             < self.n_keys else self.n_keys
                self.processed = finish
            else:
                ch = self.chunk_pop()
                if ch is not None:
                    start = ch[0]
                    finish = ch[1]
                else:
                    break

        # Give back the latest chunk if there will be over 1.5 times self.chunk in the values

        if np.sum(all_indices_for_values_with_all_studies_available < finish) > 1.5 * self.chunk_size:
            if chunk_number is None:
                self.processed = self.processed = start
                finish = start
            else:
                print ('Added back chunk {}'.format(ch))
                self.chunk_pool.append(ch)
                finish = ch[0]

        indices_for_values_with_all_studies_available = all_indices_for_values_with_all_studies_available[
            all_indices_for_values_with_all_studies_available < finish]

        keys = self.keys[indices_for_values_with_all_studies_available]

        indexes = self.values[indices_for_values_with_all_studies_available, :]

        if len(indices_for_values_with_all_studies_available) == 0:
            return None, None

        return [indexes[:, i].astype(np.int64) for i in range(self.n_study)], keys

    def get_all(self, name, nonempty=True):
        ind = self.genotype_names.index(name)
        indexes = self.values[:, ind]
        if nonempty:
            r = (self.values == -1).any(axis=1)
            if len(r) == self.values.shape[0]:
                raise ValueError('There is no common names between studies')
            indexes = indexes[~r]
        return indexes

    def load_filter_exclude(self, paths):
        self.exclude = self.load_variant_filter_dataframes(paths, (("ID",), ("CHR", "bp")))

    def load_filter_include(self, paths):
        self.include = self.load_variant_filter_dataframes(paths, (("ID",), ("CHR", "bp")))

    @staticmethod
    def load_variant_filter_dataframes(paths, column_sets=None):
        filter_dataframes = []
        for variant_filter_file in paths:
            filter_dataframe = pd.read_csv(variant_filter_file, index_col=None)
            if column_sets is not None:
                matching_column_sets = np.array(
                    [pd.Series(column_set).isin(filter_dataframe.columns).all() for column_set in column_sets]).sum()
                if matching_column_sets == 0:
                    raise ValueError('{} table should have one of the following column combinations:{}{}'.format(
                        variant_filter_file, os.linesep, (
                            " OR:{}".format(os.linesep)
                                .join(["'{}'".format("', '".join(column_set)) for column_set in column_sets]))))
                if matching_column_sets > 1:
                    raise ValueError('{} table columns matches more than one possible column combination:{}{}'.format(
                        variant_filter_file, os.linesep, (
                            " OR:{}".format(os.linesep)
                                .join(["'{}'".format("', '".join(column_set)) for column_set in column_sets]))))
            filter_dataframes.append(filter_dataframe)
        if len(filter_dataframes) == 1:
            mapper_exclude = filter_dataframes[0]
        else:
            mapper_exclude = pd.concat(
                filter_dataframes, keys=np.arange(len(paths)), names=["cohort", "row"])
        return mapper_exclude




########################################################################
########################################################################
########################################################################

class Reference(object):

    def __init__(self, name="1000Gp1v3_ref", path=None):
        self.name = name
        if path is None:
            self.path_default = os.path.join(os.environ['HASEDIR'], 'data')
        else:
            self.path_default=path
        self.path = {
            self.name: {
                'table': os.path.join(self.path_default, '{}.ref.gz'.format(self.name.rstrip("_ref"))),
                'index': os.path.join(self.path_default, '{}.ref_info.h5'.format(self.name.rstrip("_ref")))}
        }
        self.dataframe = None
        self.loaded = False
        self.chunk = 10000
        self.read = 0
        self.columns = ['ID', 'allele1', 'allele2', 'CHR', 'bp']
        self.index = None

    def load(self):
        if self.name is not None:
            if self.path.get(self.name) is not None:
                try:
                    self.dataframe = pd.read_csv(self.path[self.name]['table'], compression='gzip', sep=' ',
                                                 chunksize=self.chunk)
                except:
                    self.dataframe = pd.read_csv(self.path[self.name]['table'], sep=' ', chunksize=self.chunk)
                self.loaded = True
            else:
                if os.path.isfile(os.path.join(self.path_default, self.name)):
                    try:
                        self.dataframe = pd.read_csv(os.path.join(self.path_default, self.name), compression='gzip',
                                                     sep=' ', chunksize=self.chunk)
                    except:
                        self.dataframe = pd.read_csv(os.path.join(self.path_default, self.name), sep=' ',
                                                     chunksize=self.chunk)
                    self.loaded = True
                    self.name = os.path.basename(self.name)
                else:
                    raise ValueError('Unknown reference {}!'.format((self.name)))
        else:
            raise ValueError('Reference name is not define!')

    def load_index(self):
        try:
            self.index = pd.HDFStore(os.path.join(self.path[self.name]['index']), 'r')
        except:
            if os.path.isfile(os.path.join(self.path_default, self.name)):
                self.index = pd.HDFStore(os.path.join(self.path_default, self.name), 'r')
            else:
                raise ValueError('There is {} no index file {}'.format(self.path_default, self.name))

    def next(self):
        df = self.dataframe.get_chunk()
        self.read += df.shape[0]
        return df


########################################################################
########################################################################
########################################################################

def _get_id(object):
    if isinstance(object, type(None)):
        # If the given object is none, return an empty array
        id_p = np.array([], dtype=np.str)
    else:
        # If collect the ids from the given object
        if isinstance(object, tuple):
            id_p = np.array([], dtype=np.str)
            for i in object:
                if isinstance(i, dict):
                    id_p = np.append(id_p, i['id'])
                else:
                    id_p = np.append(id_p, i.get_id())
        else:
            id_p = object.get_id()
    return id_p.astype(np.str)


def get_identifier_array(object):
    """
    Method that returns a 2d array of samples, with the first column
    representing the sample id, and the second column representing the
    second column.
    :param object:
    :return:
    """
    if isinstance(object, type(None)):
        # If the given object is none, return an empty array
        id_p = np.array([], dtype=np.str)
    else:
        # If collect the ids from the given object
        if isinstance(object, tuple):
            id_p = np.empty((0,2), dtype=np.str)
            # For each study / cohort, the ids are appended.
            # AFAIK, no care is given to possible duplicates.
            # Here, we take care of possible duplicates by prepending
            # ids for a specific study with the study name
            for study_index, study_data in enumerate(object):
                if isinstance(study_data, dict):
                    sample_ids_in_current_study = study_data['id']
                else:
                    sample_ids_in_current_study = study_data.get_id()
                new_sample_ids = np.vstack((
                    sample_ids_in_current_study,
                    np.array([study_index] * len(sample_ids_in_current_study))
                )).T

                id_p = np.append(id_p, new_sample_ids, 0)
        else:
            id_p = object.get_id()
    return id_p.astype(np.str)


def get_intersecting_individual_indices(datasets):
    # Make sure that the length of dictionary datasets is larger than 1
    if len(datasets) < 2:
        raise ValueError("Number of datasets should be above 1")
    # Indices
    indices = dict()
    # Get the identifiers from the datasets
    identifiers = {key: get_identifier_array(dataset) for key, dataset in datasets.items()}
    # Get the keys from the identifiers
    keys = list(identifiers.keys())
    # Start the intersections
    intersection = set([tuple(ids) for ids in identifiers[keys[0]]])

    # Loop through the datasets and get the preliminary intersection
    for key in keys[1:]:
        new_set = set([tuple(ids) for ids in identifiers[key]])
        intersection = set([x for x in set(tuple(x) for x in new_set) & set(tuple(x) for x in intersection)])

    intersection = np.array(list(intersection))

    for key in keys:
        # Loop through the datasets and get the indices of the values in the
        # intersection within the original list of identifiers.
        identifier = identifiers[key][:, np.newaxis]
        intersection_indices = (identifier == intersection).all(-1)
        indices[key] = intersection_indices.argmax(axis=0)

    print ('There are {} common ids'.format(len(intersection)))
    # Could possibly save the results...
    return indices, intersection


def select_identifiers(study_index, row_index, intersection):
    # Go through all datasets, picking all samples from the intersection
    # where the study index column is equal to the requested study index.
    indices_of_selected_samples = dict()

    # Get the identifiers from the datasets
    for key, indices in row_index.items():
        indices_of_selected_samples[key] = indices[np.where(
            intersection[:,1] == str(study_index)
        )]

    return indices_of_selected_samples


def study_indexes(args=None,
                  genotype=None,
                  phenotype=None,
                  covariates=None,
                  interaction=None):
    """
	Method that matches sample identifiers between genotype, phenotype and
	covariates.single datasets.
	:param args: Maybe supposed to be the HASE output path? Currently does not do a whole lot.
	:param genotype: The genotype data (gen.folder._data)
	:param phenotype: The phenotype data (phen.folder._data)
	:param covariates: The covariate data
	:param interaction: The interaction data
	:return: A list with, for every dataset, the indices of common indices in
	the corresponding dataset.
	And, the common ids.
	"""

    id_i = _get_id(interaction)
    id_c = _get_id(covariates)
    id_p = _get_id(phenotype)
    id_g = _get_id(genotype)

    index_g = np.array([])
    index_p = np.array([])
    index_c = np.array([])
    index_i = np.array([])

    if args is not None:
        if not isinstance(args.ind_id_inc, type(None)) or not isinstance(args.ind_id_exc, type(None)):
            pass
        if not isinstance(args.cov_name_inc, type(None)) or not isinstance(args.cov_name_exc, type(None)):
            pass

    if len(id_g) != 0 and len(id_p) != 0 and len(id_c) != 0 and len(id_i) != 0:
        common_id = (np.intersect1d(np.intersect1d(np.intersect1d(id_g, id_p), id_c), id_i))

        index_g = np.array([np.where(id_g == i)[0][0] for i in common_id])
        index_p = np.array([np.where(id_p == i)[0][0] for i in common_id])
        index_c = np.array([np.where(id_c == i)[0][0] for i in common_id])
        index_i = np.array([np.where(id_i == i)[0][0] for i in common_id])

    elif len(id_i) != 0:
        raise ValueError
    elif len(id_g) != 0 and len(id_p) != 0 and len(id_c) != 0:
        common_id = (np.intersect1d(np.intersect1d(id_g, id_p), id_c))

        index_g = np.array([np.where(id_g == i)[0][0] for i in common_id])
        index_p = np.array([np.where(id_p == i)[0][0] for i in common_id])
        index_c = np.array([np.where(id_c == i)[0][0] for i in common_id])

    elif len(id_g) == 0:
        common_id = (np.intersect1d(id_p, id_c))
        index_p = np.array([np.where(id_p == i)[0][0] for i in common_id])
        index_c = np.array([np.where(id_c == i)[0][0] for i in common_id])

    elif len(id_c) == 0:
        common_id = (np.intersect1d(id_p, id_g))
        index_g = np.array([np.where(id_g == i)[0][0] for i in common_id])
        index_p = np.array([np.where(id_p == i)[0][0] for i in common_id])

    else:
        common_id = (np.intersect1d(id_g, id_c))
        index_g = np.array([np.where(id_g == i)[0][0] for i in common_id])
        index_c = np.array([np.where(id_c == i)[0][0] for i in common_id])

    print ('There are {} common ids'.format(len(common_id)))
    np.savetxt(os.path.join(os.environ['HASEOUT'], 'study_common_id.txt'), common_id, fmt='%s')
    np.savetxt(os.path.join(os.environ['HASEOUT'], 'gen_id.txt'), index_g, fmt='%s')
    np.savetxt(os.path.join(os.environ['HASEOUT'], 'phen_id.txt'), index_p, fmt='%s')
    np.savetxt(os.path.join(os.environ['HASEOUT'], 'cov_id.txt'), index_c, fmt='%s')
    np.savetxt(os.path.join(os.environ['HASEOUT'], 'interaction_id.txt'), index_i, fmt='%s')
    if len(common_id) == 0:
        exit(0)

    return [index_g, index_p, index_c, index_i], np.array(common_id)


# @timing
def merge_genotype(genotype, variant_indices, mapper, flip_flag=True):
    if variant_indices is None:
        gen = genotype[0].get_next()
        if gen is not None:
            for i in range(1, len(genotype)):
                gen = np.hstack((gen, genotype[i].get_next()))
        return gen
    else:
        if len(genotype) != len(variant_indices):
            raise ValueError('There are not equal number of genotypes and SNPs indexes {}!={}'.format(len(genotype),
                                                                                                      len(variant_indices)))
        gen = np.zeros((len(variant_indices[0]), np.sum([i.folder._data.id.shape[0] for i in genotype])))
        a = 0
        b = 0
        for i in range(0, len(genotype)):
            g = genotype[i].get(variant_indices[i], impute=False)
            print(gen.shape)
            b += g.shape[1]
            if flip_flag:
                if not mapper.encoded.get(genotype[i].folder.name, 0):
                    flip = mapper.flip[genotype[i].folder.name][variant_indices[i]]
                    flip_index = (flip == -1)
                    g = np.apply_along_axis(lambda x: flip * (x - 2 * flip_index), 0, g)

            g[np.where(variant_indices[i] == -1), :] = 0
            gen[:, a:b] = g
            a = b

        return gen


# @timing
def get_genotype(genotype, study_index, variant_indices, mapper, flip_flag=True):
    if variant_indices is None:
        gen = genotype[study_index].get_next()
        return gen
    else:
        if len(genotype) != len(variant_indices):
            raise ValueError('There are not equal number of genotypes and SNPs indexes {}!={}'.format(len(genotype),
                                                                                                      len(variant_indices)))

        # Now get the genotype for this study and perform flipping if necessary
        g = genotype[study_index].get(variant_indices[study_index], impute=False)
        print(g.shape)
        if flip_flag:
            if not mapper.encoded.get(genotype[study_index].folder.name, 0):
                flip = mapper.flip[genotype[study_index].folder.name][variant_indices[study_index]]
                flip_index = (flip == -1)
                g = np.apply_along_axis(lambda x: flip * (x - 2 * flip_index), 0, g)

        # Set all genotypes to 0 if the variant in question is not present
        # in this dataset.
        g[np.where(variant_indices[study_index] == -1), :] = 0

        return g


def check_converter(converted_folder, study_name):
    l_ind = os.listdir(os.path.join(converted_folder, 'individuals'))
    l_probes = os.listdir(os.path.join(converted_folder, 'probes'))
    l_genotype = os.listdir(os.path.join(converted_folder, 'genotype'))

    if os.path.isdir(os.path.join(converted_folder, 'tmp_files')):
        l_tmp = os.listdir(os.path.join(converted_folder, 'tmp_files'))
        if len(l_tmp) == 0:
            print ('There is no tmp files!')
        else:
            try:
                df_tmp = pd.read_csv(os.path.join(converted_folder, 'tmp_files', 'info.txt'), header=None)
            except:
                df_tmp = None

    else:
        df_tmp = None

    if len(l_ind) != 1 or '.h5' not in l_ind[0]:
        print ('There should be 1 hdf5 file in individuals folder!')

    else:
        df_ind = pd.read_hdf(os.path.join(converted_folder, 'individuals', study_name + '.h5'), 'individuals')
        if df_tmp is not None and int(df_tmp[0][0]) != df_ind.shape[0]:
            print ('Converted number of IDs {} != {} original number'.format(df_ind.shape[0], df_tmp[0][0]))
        else:
            print ('Number of individuals {} '.format(df_ind.shape[0]))
            print (df_ind.head())

    if len(l_probes) != 2:
        print ('There should be 2 files in probes folder, you have {}!'.format(len(l_probes)))
        probes = pd.HDFStore(os.path.join(converted_folder, 'probes', study_name + '.h5'), 'r')
        N_probs = probes.get_storer('probes').nrows

    else:
        probes = pd.HDFStore(os.path.join(converted_folder, 'probes', study_name + '.h5'), 'r')
        N_probs = probes.get_storer('probes').nrows

    if df_tmp is not None and int(df_tmp[0][1]) != N_probs:
        print ("Converted number of variants {} diff from  original number {}".format(N_probs, df_tmp[0][1]))

    else:
        print ('Converted number of variants {}'.format(N_probs))

        print (probes.select('probes', start=0, stop=10))
        probes.close()

    if len(l_genotype) == 0:
        print('There is no genotype data converted!')
    else:
        index = 0
        for i in l_genotype:
            n, m = h5py.File(os.path.join(converted_folder, 'genotype', i), 'r')['genotype'].shape
            index += n

        if index != N_probs:
            print('Number of variants in genotype folder {} != {} number in probes file'.format(index, N_probs))

        else:
            print ('Number of variants in genotype folder {}'.format(index))

        print (h5py.File(os.path.join(converted_folder, 'genotype', l_genotype[0]), 'r')['genotype'][...])


if __name__ == '__main__':
    print ('tools')
