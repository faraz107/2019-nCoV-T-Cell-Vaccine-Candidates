# This is a standalone code for population coverage calculation

Prerequisites:
-------------
+ Python 2.7 or higher
  * http://www.python.org/
    - Under ubuntu: sudo apt-get install python2.7

+ numpy
  * https://pypi.python.org/pypi/numpy
    - Under ubuntu: pip install numpy

+ matplotlib (tested on version 2.0.0)
  * https://matplotlib.org/downloads.html
    - Under ubuntu: pip install matplotlib==2.0.0


Installation:
------------
Unpack the tar.gz files (IEDB_Population_Coverage-1.0.tar.gz)
Run the 'configure' script to install necessary packages and requirements file.

Specifically for population_coverage
  $ tar -zxvf IEDB_Population_Coverage-1.0.tar.gz
  $ cd population_coverage
  $ ./configure


Usage:
-----
* Deplay usage
$ python calculate_population_coverage.py --help

* List all populations
$ python calculate_population_coverage.py --list

* Calculate population coverage for a given file containing epitope sequence and a list of alleles.
$ python calculate_population_coverage.py -p <population name> -c <mhc class> -f <input file path>
Example: $ python calculate_population_coverage.py -p Japan -c I -f ./test/mhci_alleles.txt

* Calculate population coverage and generate plots.
$ python calculate_population_coverage.py -p <population name> -c <mhc class> -f <input file path> --plot <output plot path>
Example: $ python calculate_population_coverage.py -p Japan -c II -f ./test/mhcii_alleles.txt --plot /home/dorjee/
