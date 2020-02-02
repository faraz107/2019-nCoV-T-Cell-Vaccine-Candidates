from __future__ import print_function
import os
import sys
import re

# adding dependencies to the Python path
project_dir = os.path.dirname(os.path.realpath(__file__))
deps = os.path.join(project_dir, "deps")
sys.path.append(os.path.join(deps, "population-coverage-pickle"))

from population_coverage_pickle import population_coverage, country_ethnicity, ethnicity as population_by_ethnicity


def read_input_file(filename):
    epitope_list = []
    allele_tuple_list = []
    try:
        for input_row in open(filename).readlines():
            epitope_list.append(input_row.strip().split()[0])
            allele_tuple_list.append(tuple(input_row.strip().split()[1].split(",")))
    except (IOError, TypeError):
        raise IOError("Oops! input file error.")
    return epitope_list, allele_tuple_list


def get_column_header(type):
    header = None

    if type == "chart_table":
        header = ("population/area", "epitope_hits", "percent_individuals", "cumulative_coverage")

    if type == "calculation_table":
        header = ("population/area", "coverage", "average_hit", "pc90")

    return header


def flatten_list(lis):
    return [item for sublis in lis for item in sublis]


def mean(lst):
    """calculates mean of the lst"""
    return float(sum(lst)) / max(len(lst), 1)


def stddev(lst):
    """returns the standard deviation of lst"""
    from math import sqrt
    mn = mean(lst)
    variance = sum([(e-mn)**2 for e in lst]) / len(lst)
    return sqrt(variance)


def into_percentage(n):
    return n*100


def print_chart_table(result=None):
    import itertools

    if not result:
        print("* No result found! *")
    else:
        for mhc_key, result in itertools.groupby(result, lambda item: item["mhc_class"]):
            coverage_list = []
            hit_list = []
            pc90_list = []

            print("class {}".format(mhc_key))

            # print summary table for all individual coverage result for each mhc_class
            calculation_header = get_column_header("calculation_table")
            print("\t".join(calculation_header))
            individual_class_result = list(result)

            for rd in individual_class_result:
                population = rd.get("population")
                coverage = rd.get("coverage")
                coverage_list.append(coverage)
                hit = rd.get("average_hit")
                hit_list.append(hit)
                pc90 = rd.get("pc90")
                pc90_list.append(pc90)

                print("{}\t{}%\t{}\t{}".format(
                    population,
                    round(into_percentage(coverage), 2),
                    round(hit, 2),
                    round(pc90, 2))
                )

            print("{}\t{}%\t{}\t{}".format(
                "average",
                round(into_percentage(mean(coverage_list)), 2),
                round(mean(hit_list), 2),
                round(mean(pc90_list), 2))
            )

            print("{}\t{}%\t{}\t{}\n".format(
                "standard_deviation",
                round(into_percentage(stddev(coverage_list)), 2),
                round(stddev(hit_list), 2),
                round(stddev(pc90_list), 2))
            )

            # print individual population coverage result
            chart_header = get_column_header("chart_table")
            print("\t".join(chart_header))
            for rd in individual_class_result:
                population = rd.get("population")
                epitope_hits = rd.get("epitope_hits")
                percent_individuals = [round(pi*100, 2) for pi in rd.get("percent_individuals")]
                cumulative_coverage = [round(cc, 2) for cc in rd.get("cumulative_coverage")]
                for row in zip(epitope_hits, percent_individuals, cumulative_coverage):
                    print("{}\t{}\t{}\t{}".format(population, row[0], row[1], row[2]))
            print("\n")


def replace_w_char(s):
    s = re.sub(r";|\s|-|,", "_", s)
    return re.sub(r"_+", "_", s)


def generate_plot(result=None, plot_path=None):
    import numpy as np
    import matplotlib.cbook
    matplotlib.use('Agg')
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib.figure import Figure

    if not os.path.isdir(plot_path):
        print("plot directory doesn't exist!")
        sys.exit(0)

    for rd in result:
        fig = Figure(figsize=(6.5, 4.0), facecolor="#F5F5F5")
        ax = fig.add_subplot(1, 1, 1)

        population = rd.get("population")
        mhc_class = rd.get("mhc_class")
        epitope_hits = rd.get("epitope_hits")
        percent_individuals = rd.get("percent_individuals")
        percent_individuals = [pi*100 for pi in percent_individuals]
        cumulative_coverage = rd.get("cumulative_coverage")

        ax.bar(epitope_hits, percent_individuals, color="#6698FF", align="center")
        ax.set_title("{} - Class {} Coverage".format(population, mhc_class), fontsize=11)
        ax.set_xlabel("Number of epitope hits/HLA combination recognized", fontsize=9)
        ax.set_ylabel("Percent of individuals", fontsize=9)

        ax.set_xticks(np.arange(len(epitope_hits)))
        ax.set_xticklabels(epitope_hits)

        ax.set_ylim(ymax=max(percent_individuals) + 15)
        ax.xaxis.set_tick_params(labelsize=9)
        ax.yaxis.set_tick_params(labelsize=9)
        ax.grid(True)

        ax2 = ax.twinx()
        ax2.set_xlim(left=-0.55)
        ax2.set_ylim(0, 100)
        ax2.set_ylabel("Cumulative percent of population coverage", fontsize=9)
        ax2.yaxis.set_tick_params(labelsize=9)
        majorLocator = MultipleLocator(10)
        majorFormatter = FormatStrFormatter("%d")
        minorLocator = MultipleLocator(5)
        ax2.yaxis.set_major_locator(majorLocator)
        ax2.yaxis.set_major_formatter(majorFormatter)
        ax2.axhline(y=90, color="r", linewidth=1, label="Threshold")
        ax2.plot(ax.get_xticks(), cumulative_coverage, linestyle="-", marker="o", markerfacecolor="yellow", linewidth=1.5)

        canvas = FigureCanvas(fig)
        canvas.print_figure("{}/popcov_{}_{}.png".format(plot_path, replace_w_char(population.lower()), replace_w_char(mhc_class.lower())))

    print("* A plot has been generated in '{}' directory with <population>_<mhc_class> suffixed.".format(plot_path))


def get_available_allele_names(population=None, mhc_class=None):
    class_population_alleles = population_coverage.get(mhc_class).get(population)
    for locus, allele_genotype in class_population_alleles.iteritems():
        for allele, genotype in allele_genotype:
            print("{}\t{}".format(allele, genotype))


def get_population_list():
    print("List of population by area:\n---")
    for region, population_ethnicity in country_ethnicity.iteritems():
        print("{}".format(region))
        for population, ethnicity in population_ethnicity.iteritems():
            print("\t{}".format(population))
            for e in ethnicity:
                print("\t\t{}".format(e))
    print("\n")
    print("List of ethnicity:\n---")
    for _ethnicity in population_by_ethnicity:
        print("{}".format(_ethnicity))


def get_all_alleles():
    """ returns a dictionary of alleles keyed by mhc class """
    all_class_alleles = {}
    for mhc_class, population_map in population_coverage.iteritems():
        all_alleles = []
        for population, locus_map in population_map.iteritems():
            for lm_list in locus_map.values():
                all_alleles.extend([lm[0] for lm in lm_list])
        all_class_alleles[mhc_class] = list(set(all_alleles))
    return all_class_alleles


def get_population_coverage_dict():
    return country_ethnicity


def locus_by_allele_name(mhc_class, population, query_alleles):
    for lkey, agvalue in population_coverage.get(mhc_class).get(population).iteritems():
        if [tup for tup in agvalue for qa in query_alleles if qa in tup]:
            return lkey

