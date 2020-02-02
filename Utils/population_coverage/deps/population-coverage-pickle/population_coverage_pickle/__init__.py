import pickle
from pkg_resources import resource_filename  # @UnresolvedImport

package_name = "population_coverage_pickle"
pickle_filename = "population_genotype_map.p"
pickle_file_path = resource_filename(package_name, pickle_filename)

with open(pickle_file_path, "rb") as pickle_file:
    population_coverage = pickle.load(pickle_file)
    country_ethnicity = pickle.load(pickle_file)
    ethnicity = pickle.load(pickle_file)
