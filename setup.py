#! /usr/bin/env python
from setuptools import setup

setup(name='seq-annot',
      version='0.6.0',
      packages=['seq_annot',],
      description='Tools that fascilitate the annotation and functional '
          'comparison of metagenomes',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Software Development :: Libraries :: Python Modules'
      ],
      keywords='bioinformatics sequence profiling annotation',
      url='https://github.com/Brazelton-Lab/seq-annot/',
      download_url = 'https://github.com/Brazelton-Lab/seq-annot/archive/v0.6.0.tar.gz',
      author='Christopher Thornton',
      author_email='christopher.thornton@utah.edu',
      license='GPLv3',
      include_package_data=True,
      zip_safe=False,
      install_requires=['bio_utils', 'HTSeq', 'arandomness'],
      entry_points={
          'console_scripts': [
              'compare_features = seq_annot.compare:main',
              'count_features = seq_annot.count:main',
              'annotate_features = seq_annot.annotate:main',
              'screen_features = seq_annot.screen:main',
              'combine_features = seq_annot.combine:main',
              'reldbs = seq_annot.reldbs:main',
              'colocate_features = seq_annot.colocate:main'
          ]
      }
      )
