#!/usr/bin/env python3

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import with_statement

import errno
import os
import sys

import numpy as np
import pandas as pd
import tables


class GenotypeDataException(Exception):
    pass


class GenotypeData:
    _genotype_matrix = None  # type: np.ndarray(dtype=np.float) # shape = genotype,individual
    _alleles = None  # type: list[list[str]]
    _variant_data_frame = None  # type: pd.DataFrame['ID', "CHR", "bp", "allele1", "allele2"]

    def __init__(self, genotype_matrix, variant_data_frame, individual_ids):
        self.number_of_variants = genotype_matrix.shape[0]
        self.number_of_individuals = genotype_matrix.shape[1]
        self._genotype_matrix = genotype_matrix
        self._variant_data_frame = variant_data_frame
        self._individuals_data = pd.DataFrame(
            individual_ids,
            columns=["individual"])

    def get_individuals_data_frame(self):
        return self._individuals_data

    def get_variant_chunks(self):
        # Initialize counter for the number of snps
        if self.number_of_variants != self._variant_data_frame.shape[0]:
            raise GenotypeDataException(
                "Variant data frame does not have the expected length ({} vs {} respectively)"
                    .format(self._variant_data_frame.shape[0],
                            self.number_of_variants))

        return [self._variant_data_frame]

    def get_alleles(self, variant_index):
        return self._variant_data_frame.iloc[variant_index][['str_allele1', 'str_allele2']].values.flatten().tolist()

    def get_dosages(self, variant_index):
        return self._genotype_matrix[variant_index, :]

    def to_csv(self, path):
        pd.DataFrame(data=self._genotype_matrix,
                     index=self._variant_data_frame.ID,
                     columns=self._individuals_data.individual).to_csv(path, sep="\t")


class HaseHDF5WriterException(Exception):
    pass


class FileExistsError(OSError):
    def __init__(self, msg):
        super(FileExistsError, self).__init__(errno.EEXIST, msg)


class HaseHDF5Writer(object):
    def __init__(self, path, chunk_size, study_name):
        self.chunk_size = chunk_size
        self.abs_path = os.path.realpath(os.path.expanduser(path))
        self.study_name = study_name
        self.bad_variant_indices = list()

        self.probes_directory_path = os.path.join(self.abs_path, "probes")
        try:
            print("Creating probes folder at {}...".format(self.probes_directory_path))
            os.makedirs(self.probes_directory_path)
        except FileExistsError as e:
            raise HaseHDF5WriterException("Directory '{}' already exists"
                                          .format(self.probes_directory_path), e)

        self.individuals_directory_path = os.path.join(self.abs_path, "individuals")
        try:
            print("Creating individuals folder at {}...".format(self.individuals_directory_path))
            os.mkdir(self.individuals_directory_path)
        except FileExistsError as e:
            raise HaseHDF5WriterException("Directory '{}' already exists"
                                          .format(self.individuals_directory_path), e)

        self.genotype_directory_path = os.path.join(self.abs_path, "genotype")
        try:
            print("Creating genotype folder at {}...".format(self.genotype_directory_path))
            os.mkdir(self.genotype_directory_path)
        except FileExistsError as e:
            raise HaseHDF5WriterException("Directory '{}' already exists"
                                          .format(self.genotype_directory_path), e)

    def write_genotype_matrix(self, genotype_data):
        """

        :type genotype_data: GenotypeData
        """
        self._write_probes(genotype_data)
        print("Completed probes file")
        self._write_individuals(genotype_data)
        print("Completed individuals file")
        self._write_genotype(genotype_data)
        print("Completed genotype file")

    def _write_probes(self, genotype_data):
        probes_path = os.path.join(self.probes_directory_path, self.study_name + '.h5')
        if os.path.isfile(probes_path):
            raise HaseHDF5WriterException("File '{}' already exists".format(probes_path))

        hash_table = {'keys': np.array([], dtype=np.int), 'allele': np.array([])}

        variant_index = 0
        for i, variant_chunk in enumerate(genotype_data.get_variant_chunks()):

            alleles1 = list()
            alleles2 = list()
            chunked_variant_index = 0
            variant_chunk_length = len(variant_chunk)
            bad_variant_indices_probes_chunk = list()
            while chunked_variant_index < variant_chunk_length:
                alleles = genotype_data.get_alleles(variant_index)
                if len(alleles) != 2:
                    print("Not found 2 alleles for variant '{}': discarding variant...".format(
                        variant_chunk["ID"][chunked_variant_index]), file=sys.stderr)
                    print(variant_index, chunked_variant_index)
                    self.bad_variant_indices.append(variant_index)
                    bad_variant_indices_probes_chunk.append(chunked_variant_index)

                    # if len(alleles) == 1:
                    #     print("len alleles is 1", file=sys.stderr)
                    #     alleles.append(alleles[0])
                    # elif len(alleles) == 0:
                    #     print("len alleles is 0", file=sys.stderr)
                    #     alleles = ["N", "N"]
                else:
                    alleles1.append(alleles[0])
                    alleles2.append(alleles[1])

                variant_index += 1
                chunked_variant_index += 1

            # Drop every variant that does not have exactly 2 alleles
            filtered_variant_chunk = variant_chunk.drop(
                bad_variant_indices_probes_chunk)

            filtered_variant_chunk = filtered_variant_chunk.reset_index(drop=True)

            alleles1_series = pd.Series(alleles1)
            alleles2_series = pd.Series(alleles2)

            hash_1 = alleles1_series.apply(hash)
            hash_2 = alleles2_series.apply(hash)

            k, indices = np.unique(np.append(hash_1, hash_2), return_index=True)
            s = np.append(alleles1_series, alleles2_series)[indices]
            ind = np.invert(np.in1d(k, hash_table['keys']))
            hash_table['keys'] = np.append(hash_table['keys'], k[ind])
            hash_table['allele'] = np.append(hash_table['allele'], s[ind])
            filtered_variant_chunk["allele1"] = hash_1
            filtered_variant_chunk["allele2"] = hash_2

            filtered_variant_chunk.to_hdf(probes_path,
                                          data_columns=["CHR", "bp", "ID", 'allele1', 'allele2'],
                                          key='probes', format='table', append=True,
                                          min_itemsize=25, complib='zlib', complevel=9, dropna=True)

            print("Wrote {} variants to probes file".format(variant_index))

        pd.DataFrame.from_dict(hash_table).to_csv(
            os.path.join(self.probes_directory_path, self.study_name + '_hash_table.csv.gz'),
            index=False, compression='gzip', sep='\t')

    def _write_individuals(self, genotype_data):
        individuals_path = os.path.join(self.individuals_directory_path, self.study_name + '.h5')
        if os.path.isfile(individuals_path):
            raise HaseHDF5WriterException("File '{}' already exists".format(individuals_path))

        genotype_data.get_individuals_data_frame().to_hdf(
            individuals_path, key='individuals', format='table',
            min_itemsize=25, complib='zlib', complevel=9)

    def _write_genotype(self, genotype_data):

        number_of_chunks = (genotype_data.number_of_variants // self.chunk_size) + 1
        for chunk_index in xrange(number_of_chunks):
            start = chunk_index * self.chunk_size
            end = min((chunk_index + 1) * self.chunk_size, genotype_data.number_of_variants)
            dosage_matrix = np.empty((end - start, genotype_data.number_of_individuals))

            print("Loading {}-{} variants to write to chunk {} out of {} total chunks".format(
                start, end, chunk_index, number_of_chunks))

            if dosage_matrix.shape[0] == 0:
                print("Empty chunk. Skipping...")
                continue

            # Get the dosages for every variant within the
            for chunked_variant_index, variant_index in enumerate(xrange(start, end)):
                dosage_matrix[chunked_variant_index,] = genotype_data.get_dosages(variant_index)

            # Drop every variant that did not have two alleles.
            bad_variant_indices_chunk = list()
            for bad_variant_index in self.bad_variant_indices:
                if start <= bad_variant_index < end:
                    bad_variant_indices_chunk.append(bad_variant_index - start)
            print(bad_variant_indices_chunk)
            dosage_matrix = np.delete(dosage_matrix, bad_variant_indices_chunk, 0)

            h5_gen_file = tables.open_file(
                os.path.join(self.genotype_directory_path, str(chunk_index) + '_' + str(self.study_name) + '.h5'), 'w',
                title=self.study_name)

            atom = tables.Float16Atom()
            genotype = h5_gen_file.create_carray(h5_gen_file.root, 'genotype', atom,
                                                 (dosage_matrix.shape),
                                                 title='Genotype',
                                                 filters=tables.Filters(complevel=9, complib='zlib'))

            genotype[:] = dosage_matrix
            h5_gen_file.close()
        print("Discarded {} variants that did not have two alleles".format(
            len(self.bad_variant_indices)), file=sys.stderr)
