Analyzer of periodic repetitions in sequences of DNA
====================================================

Use aprdna to look for periodic repetitions of dinucleotides in sequences of DNA, in the GenBank format.

Install
-------

Use pip to install the program:

.. code-block:: bash

  $ pip install aprdna

Usage
-----

usage: aprdna [-h] [-a NAME] [-p PERIOD] [-t TOLERANCE] [-m LENGTH] [-r THROTTLE] path

Search for periodic repetitions of dinucleotides. Given a GenBank sequence file,
it will search it for repetitions of each of all possible dinucleotides with a
periodicity of PERIOD +/- TOLERANCE. DNA regions in which there are LENGTH or more
consecutive repetitions are considered a match. The output is the percentage of the
sequence corresponding to matching regions for each dinucleotide, and PDFs comparing
the matching regions with other features published with the given sequence files.


positional arguments:
  path                  Path to GenBank Sequence file

optional arguments:
  -h, --help            show help message and exit
  -p PERIOD, --period PERIOD
                        Periodicity of the repetition
  -t TOLERANCE, --tolerance TOLERANCE
                        Tolerance window for the repetition
  -m LENGTH, --length LENGTH
                        Minimal repetition length
  -r THROTTLE, --throttle THROTTLE
                        Stop after so many fragments
