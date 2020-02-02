import os
import tempfile
from util import *
from population_calculation import PopulationCoverage


def calculate(epitope=None, allele=None, population=None, mhc_class=None):
    tmpfile = tempfile.NamedTemporaryFile(suffix='.txt', delete=False)
    try:
        for ea in zip(epitope, allele):
            tmpfile.write("{}\t{}\n".format(ea[0], ",".join(ea[1])))
        tmpfile.seek(0)
    finally:
        tmpfile.close()

    pc = PopulationCoverage()
    result = pc.calculate_coverage(population=population, mhc_class=mhc_class, filename=tmpfile.name)

    if not result:
        msg = 'Error calling population coverage standalone:\n{}, {}, {}'.format(population, mhc_class, tmpfile.name)
        raise Exception(msg)

    # delete the temp file
    os.remove(tmpfile.name)

    return result


def get_population_coverage_info():
    return get_population_coverage_dict()
