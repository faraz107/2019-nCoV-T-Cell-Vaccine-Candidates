"""
Created on: 03/06/2017

@author: Dorjee Gyaltsen
@brief: calculates population coverage given list of population(s), mhc class option(s) and
        a file containing eptiope-allele combination - to calculate_coverage() method
        file syntax: <epitope name 1> \t <allele name 1>,<allele name 2>,...
"""
from __future__ import print_function

import copy
from itertools import product, permutations
from collections import OrderedDict
from util import *

from population_coverage_pickle import population_coverage


class PopulationCoverage(object):

    def __init__(self):
        self.total_genotype_map = OrderedDict()
        self.total_genotype_list = []
        self.allele_tuple_list = None
        self.input_epitope_allele_list = None

    def set_file_content(self, filename):
        self.input_epitope_allele_list = zip(*read_input_file(filename))

    @staticmethod
    def validate_input_file(population_ethnicity=None):
        # all_input_alleles = [allele_names for epitope, allele_names in self.input_epitope_allele_list]
        # flatten list of tuples
        # all_input_alleles = flatten_list(all_input_alleles)

        population_ethnicity_list = []
        # all_unique_mhc_i_ii_alleles = []
        for population_data in population_coverage.values():
            population_ethnicity_list.append(population_data.keys())

        population_ethnicity_list = flatten_list(population_ethnicity_list)
        if population_ethnicity:
            if population_ethnicity not in population_ethnicity_list:
                print("population/ethnicity doesn't exist!")
                sys.exit(0)

    @staticmethod
    def dictionarized_allele_map(allele_map):
        return OrderedDict([(ag_tuple[0], ag_tuple[1]) for ag_tuple in allele_map])

    def calculate_coverage(self, query_by="area_country_ethnicity", population=None, mhc_class=None, filename=None):
        """ @brief: Returns 'result' - a list of dictionary(ies) keyed by individual result name and
                            'negative' - if any, a list of invalid input(s)
            @note: This is start method for population_coverage commandline. population, mhc_class and filename
                   are the required required parameters
        """
        result = []
        negative = []

        #  list of available mhc_class options - hard-coded
        available_mhc_class = ["I", "II", "combined"]

        for _mhc_class in mhc_class:
            if _mhc_class not in available_mhc_class:
                print("-c '{}' value is not an option".format(_mhc_class))
                sys.exit(0)

        class_population_combo = list(product(mhc_class, population))

        self.set_file_content(filename)
        self.validate_input_file()

        population_map = self.get_population_map(class_population_combo)
        for population, class_map in population_map.iteritems():
            for mhc_class, coverage_data in class_map.iteritems():
                frequency = self.get_frequency(coverage_data)
                result_map = {}
                if frequency:
                    result_map = frequency[0]

                for frequency_map in frequency[1:]:
                    merged_locus = self.merge_loci(result_map, frequency_map)
                    result_map = merged_locus

                if len(result_map) > 1:
                    frequency_map = self.compute_graph_frequency(result_map)
                    if frequency_map:
                        frequency_map.update({"mhc_class": mhc_class})
                        frequency_map.update({"population": population})
                    result.append(frequency_map)

        _mhc_population_values = [(rd.get("mhc_class"), rd.get("population")) for rd in result]
        for _negative_class, _negative_population in list(set(class_population_combo) - set(_mhc_population_values)):
            negative.append({
                'cumulative_coverage': [],
                'mhc_class': _negative_class,
                'average_hit': 0.0,
                'coverage': 0.0,
                'pc90': 0.0,
                'epitope_hits': [],
                'percent_individuals': [],
                'population': _negative_population
            })

        return result, negative

    def get_frequency(self, coverage_data):
        """

        """
        adjusted_genotype = self.get_adjusted_genotype(coverage_data)

        # locus_map is a dictionary of list of dictionaries where outer dictionary is key by locus names, and
        # values are a list of key-value pairs of each allele_name and genotype frequency
        #   syntax: {'locus': [{'allele_name': 'genotype'}, {...}, ...]}
        locus_map = OrderedDict()
        for locus, tupleized_allele_map in adjusted_genotype.iteritems():
            allele_map = self.dictionarized_allele_map(tupleized_allele_map)
            locus_map[locus] = allele_map

        total_hits = self.count_hits(locus_map)
        frequency = self.compute_frequency(locus_map, total_hits)
        return frequency

    def get_population_map(self, class_population_combo):
        d = {}
        for mhc_class, population in class_population_combo:
            d.setdefault(population, []).append(mhc_class)

        population_map = {}
        for _population, _mhc_list in d.iteritems():
            class_map = {}
            for _mhc in _mhc_list:
                if _mhc != "combined":
                    coverage_data = self.get_coverage_data(_mhc, _population)
                    class_map[_mhc] = coverage_data
                else:
                    # if any(_m in class_map for _m in ("I", "II")):
                    #     _class_map = copy.deepcopy(class_map)
                    #     # print("deep copy: ", _class_map)
                    #     class_map[_mhc] = filter(None, _class_map.values())
                    # else:
                    _combined_coverage_data = []
                    _combined_mhc = ["I", "II"]
                    for _cmhc in _combined_mhc:
                        coverage_data = self.get_coverage_data(_cmhc, _population)
                        _combined_coverage_data.append(coverage_data)
                    _combined_coverage_data = filter(None, _combined_coverage_data)

                    _coverage_data = {}
                    for cd in _combined_coverage_data:
                        _coverage_data.update(cd)
                    class_map[_mhc] = OrderedDict(sorted(_coverage_data.items()))

            class_map = OrderedDict((k, v) for k, v in class_map.items() if v is not None)

            # for _mhc in _mhc_list:
            #     if _mhc not in class_map:
            #         class_map.update({
            #             'cumulative_coverage': [],
            #             'mhc_class': _mhc,
            #             'average_hit': 0.0,
            #             'coverage': 0.0,
            #             'pc90': 0.0,
            #             'epitope_hits': [],
            #             'percent_individuals': [],
            #             'population': _population
            #         })

            population_map[_population] = class_map

        return population_map

    def get_coverage_data(self, _mhc, _population):
        coverage_data = population_coverage.get(_mhc).get(_population)
        return coverage_data

    @staticmethod
    def merge_loci(result_map, frequency_map):
        merged_locus = {}
        for hit1, genotype1 in result_map.iteritems():
            for hit2, genotype2 in frequency_map.iteritems():
                total_hit = hit1 + hit2
                total_frequency = genotype1 * genotype2
                if total_hit not in merged_locus:
                    merged_locus[total_hit] = total_frequency
                else:
                    value = merged_locus.get(total_hit)
                    merged_locus[total_hit] = total_frequency + value
        return merged_locus

    def count_hits(self, locus_map):
        total_hits = OrderedDict()
        for locus, ag_map in locus_map.iteritems():
            allele_list = sorted(ag_map.keys())
            for allele_name in allele_list:
                number_hit = 0
                for epitope, tupleized_allele in self.input_epitope_allele_list:
                    if tupleized_allele.count(allele_name):
                        number_hit += 1
                total_hits[allele_name] = number_hit
        return total_hits

    def compute_graph_frequency(self, merged_data):
        graph_output = {}

        epitope_hits = merged_data.keys()
        y = merged_data.values()

        last_epitope_hit = len(epitope_hits)

        if epitope_hits:
            last_epitope_hit = epitope_hits[-1]

        _x = []
        _y = []

        for i in range(last_epitope_hit + 1):
            _x.append(i)
            if i in epitope_hits:
                _y.append(merged_data[i])
            else:
                _y.append(0.0)
        _merged_data = dict(zip(_x, _y))

        # TODO: code below must change according to above _x, _y values
        _merged_data = dict((k, v * 100) for k, v in _merged_data.iteritems())
        percent_individuals = _merged_data.values()

        cumulative_coverage = []
        for k, v in sorted(_merged_data.iteritems(), reverse=True):
            if not cumulative_coverage:
                cumulative_coverage.append(v)
            else:
                new_v = v + cumulative_coverage[-1]
                cumulative_coverage.append(new_v)

        pc90 = self.calculate_pc90(cumulative_coverage)
        coverage = self.calculate_frequency_coverage(_x, _y)
        average_hit = self.compute_average_epitope_hit(_x, _y)

        graph_output.update({
            "epitope_hits": _x,
            "percent_individuals": _y,
            "cumulative_coverage": cumulative_coverage[::-1],
            "pc90": pc90,
            "coverage": coverage,
            "average_hit": average_hit
        })

        return graph_output

    @staticmethod
    def calculate_pc90(data):
        data = data[::-1]
        index1 = 0
        index2 = len(data) - 1

        if data[len(data) - 1] >= 90:
            pc90 = len(data) - 1

        for i in range(len(data) - 2, -1, -1):
            if data[i] == 90:
                pc90 = i

            if data[i] > 90:
                index1 = i
                break

        index2 = index1 + 1

        if index2 > len(data) - 1:
            index2 = len(data) - 1

        cur_value = data[index2]

        for i in range(index2, len(data)):
            if data[i] == cur_value:
                index2 = i

            if data[i] < cur_value:
                break

        # TODO: check the value of index1 and index2
        _range = data[index1] - data[index2]
        _fraction = (90 - data[index2]) / _range
        pc90 = index2 - (_fraction * (index2 - index1))

        return pc90

    @staticmethod
    def calculate_frequency_coverage(x, y):
        frequency = 0
        for i in x:
            if x[i] != 0:
                frequency += y[i]
        return frequency

    @staticmethod
    def compute_average_epitope_hit(x, y):
        n = 0
        sum = 0
        for i in x:
            sum += x[i] * y[i]
            n += y[i]

        if n == 0:
            return 0.0
        return sum/n

    def get_adjusted_genotype(self, coverage_data):
        """ """
        adjusted_coverage_data = OrderedDict()

        for locus, allele_genotype_map in coverage_data.iteritems():
            allele_list = []
            genotype_list = []
            for allele, genotype in allele_genotype_map:
                allele_list.append(allele)
                genotype_list.append(genotype)

            total_genotype = sum(genotype_list)
            self.total_genotype_map[locus] = total_genotype
            self.total_genotype_list.append((locus, total_genotype))

            if total_genotype > 1:
                genotype_list = [g/sum(genotype_list) for g in genotype_list]

            if total_genotype < 1:
                diff = 1 - sum(genotype_list)
                allele_list.append("UNKNOWN")
                genotype_list.append(diff)

            adjusted_coverage_data[locus] = zip(allele_list, genotype_list)
        return adjusted_coverage_data

    def compute_frequency(self, locus_map, total_hits):

        # 'total_hits' is an ordered dictionary of total number of hits keyed by allele name
        # for a particular locus (eg: 'HLA-A') in question (syntax: [{'allele': 'number of hits'}, ...]

        # 'self.ag_map_list' is a list of ordered dictionary(ies) of genotype frequency values
        # keyed by allele name for a locus (syntax: [{'allele': 'genotype'}, ...]

        total_locus_frequency = []
        # for locus_map in self.ag_map_list:
        for locus, ag_map in locus_map.iteritems():
            allele_list = sorted(ag_map.keys())
            locus_frequency = {}
            for i in range(len(allele_list)):
                hit1 = total_hits[allele_list[i]]
                genotype1 = ag_map[allele_list[i]]

                for j in range(len(allele_list)):
                    hit2 = total_hits[allele_list[j]]
                    genotype2 = ag_map[allele_list[j]]

                    # if it's homozygous (same allele)
                    total_hit = hit1

                    if i != j:
                        total_hit = hit1 + hit2

                    total_frequency = genotype1 * genotype2

                    if total_hit not in locus_frequency:
                        locus_frequency[total_hit] = total_frequency
                    else:
                        value = locus_frequency.get(total_hit)
                        locus_frequency[total_hit] = total_frequency + value

            total_locus_frequency.append(locus_frequency)
        return total_locus_frequency
