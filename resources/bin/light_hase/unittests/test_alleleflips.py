#!/usr/bin/env python3

"""
Test script for the classic meta analysis
"""

import unittest
import os

from tools import mapper


class TestAlleleFlips(unittest.TestCase):

    def test_snp_mapping(self):
        genotype_directory = os.path.join(os.path.dirname(__file__), "resources", "alleleflips")
        reference_directory = os.path.join(os.path.dirname(__file__), "resources", "alleleflips", "hase_reference")

        out_dir = os.path.join(os.path.dirname(__file__), "unittests", "alleleflips", "out")

        mapper.main(["-g", genotype_directory,
                     "-study_name", "GTEx_2017-06-05_v8_EUR",
                     "-ref_name", "1000G-30x",
                     "-ref_path", reference_directory,
                     "-o", out_dir])