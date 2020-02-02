"""
This script is to test the pickle file generated with 'pickle_for_sqlite_data.py'
"""

import pickle
from itertools import product

import unittest

pickle_filename = "population_genotype_map.p"

with open(pickle_filename, "rb") as pickle_file:
    population_coverage = pickle.load(pickle_file)
    country_ethnicity = pickle.load(pickle_file)
    ethnicity = pickle.load(pickle_file)

class TestPickleFile(unittest.TestCase):

    def test_none_result(self):
        _mhc_class = ['I']
        _population = ['Algeria']

        class_population_combo = list(product(_mhc_class, _population))
        for mhc_class, population in class_population_combo:
            result = population_coverage.get(mhc_class).get(population)
            self.assertIsNone(result)

    def test_basic(self):
        _mhc_class = ['I']
        _population = ['Japan']

        class_population_combo = list(product(_mhc_class, _population))
        for mhc_class, population in class_population_combo:
            result = population_coverage.get(mhc_class).get(population)
            self.assertIsNotNone(result)
            self.assertIsInstance(result, dict)
            self.assertEqual(result.keys(), ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-G'])