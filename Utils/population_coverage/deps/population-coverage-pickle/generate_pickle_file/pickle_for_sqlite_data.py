import pickle
import sqlite3 as lite
from collections import OrderedDict

DB_PATH = './population_coverage.sqlite'
conn = lite.connect(DB_PATH)
cur = conn.cursor()
cur.execute("SELECT class, population FROM population_genotype ORDER BY id;")
rows = cur.fetchall()

d = {}
for row in rows:
    d.setdefault(row[0], []).append(row[-1])
filtered_d = {k: list(set(v)) for k, v in d.items()}
print("Finished creating class-population dictionary!")

population_genotype_map = {}
print("Creating population genotype map...")
for mhc_class, populations in filtered_d.iteritems():
    print("{}:".format(mhc_class))

    population_map = {}
    for population in populations:
        print("\t-{}".format(population))
        cur.execute("""SELECT locus, allele, genotype FROM population_genotype where class = ? and population = ? ORDER BY id;""", (mhc_class, population,))
        lag_rows = cur.fetchall()

        locus_map = {}
        if lag_rows:
            for lag in lag_rows:
                locus_map.setdefault(lag[0], []).append(lag[1:])
        locus_map = OrderedDict(sorted(locus_map.items()))
        population_map[population] = locus_map

    population_map = OrderedDict(sorted(population_map.items()))
    population_genotype_map.update({mhc_class: population_map})

# print(population_genotype_map)

cur.execute("""SELECT distinct(ethnicity), country FROM population_genotype where population_type="Country_ethnicity" order by country;""")
rows2 = cur.fetchall()

d2 = OrderedDict()
for row in rows2:
    d2.setdefault(row[1], []).append("{} {}".format(row[1], row[0]))
# print(d2)

cur.execute("""SELECT distinct(ethnicity), country, area FROM population_genotype where population_type="Country_ethnicity" order by country;""")
ec_rows = cur.fetchall()
# print(ec_rows)

country_ethnicity_map = OrderedDict()
for row in ec_rows:
    country_ethnicity_map.setdefault(row[2], {}).update({row[1]: d2[row[1]]})
# print(country_ethnicity_map)

sortorder = ["East Asia", "Northeast Asia", "South Asia", "Southeast Asia", "Southwest Asia",
             "Europe", "East Africa", "West Africa", "Central Africa", "North Africa", "South Africa",
             "West Indies", "North America", "Central America", "South America", "Oceania"]
sorted_country_ethnicity_map = OrderedDict([(k, country_ethnicity_map[k]) for k in sortorder])
# print(sorted_country_ethnicity_map)

ethnicity_list = ["Amerindian", "Arab", "Asian", "Australian Aborigines", "Austronesian", "Berber",
                  "Black", "Caucasoid", "Hispanic", "Jew", "Kurd", "Melanesian", "Mestizo", "Micronesian",
                  "Mixed", "Mulatto", "Oriental", "Persian", "Polynesian", "Siberian"]

with open("population_genotype_map.p", "wb") as pickle_file:
    pickle.dump(population_genotype_map, pickle_file)
    pickle.dump(sorted_country_ethnicity_map, pickle_file)
    pickle.dump(ethnicity_list, pickle_file)
print "All done!"