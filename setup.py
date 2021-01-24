from setuptools import setup, find_packages


version = '0.1.0'


requires = [
    'biopython == 1.78',
    'reportlab == 3.5.60',
]

long_description = open('README.rst').read()

setup(name='aprdna',
      version=version,
      description="Analizer of periodic repetitions in DNA",
      long_description=long_description,
      classifiers=[
          "Programming Language :: Python",
      ],
      keywords='',
      author='Enrique Pérez',
      author_email='enriquepablo@gmail.com',
      license='gpl',
      packages=find_packages('src'),
      package_dir={'': 'src'},
      include_package_data=True,
      zip_safe=False,
      install_requires=requires,
      tests_require=[],
      extras_require={},
      entry_points={
          'console_scripts': ['aprdna=aprdna.script:main'],
      }
      )
