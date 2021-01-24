Analyzer of periodic repetitions in sequences of DNA
====================================================

Use aprdna to look for periodic repetitions in sequences of DNA, in the GenBank format.

Install
-------

Use pip to install the program:

.. code-block:: bash

  $ pip install aprdna

Usage
-----

usage: aprdna [-h] [-a NAME] [-p PERIOD] [-t TOLERANCE] [-m LENGTH] [-r THROTTLE] path

Search for repetitions of nucleotides.

positional arguments:
  path                  Path to Sequence file

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
