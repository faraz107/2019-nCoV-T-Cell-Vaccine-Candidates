Provides population coverage data in the form of a dictionary: {
        <mhc_class>: {
            <population>: {
                <locus>: [("<allele_name>", <genotype>), (...)]
            }
        }
}

Example use:
    from population_coverage_pickle import population_coverage, country_ethnicity, ethnicity
