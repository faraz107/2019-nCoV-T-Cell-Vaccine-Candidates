import os
from setuptools import setup

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

with open(os.path.join(os.path.dirname(__file__), "README.rst")) as readme:
    README = readme.read()

setup(
    name="population-coverage-pickle",
    version="1.1.2",
    packages=["population_coverage_pickle"],
    package_data={
        "population_coverage_pickle": ["population_genotype_map.p"]
    },
    description="PyPA package providing data for population coverage tool.",
    long_description=README,
    classifiers=[
        "Intended Audience :: Developers",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 2",
    ]
)
