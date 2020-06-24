#!/usr/bin/env python
from setuptools import setup, find_packages

entry_points = {
    'console_scripts': [
        'robo-rx=robotools.rx:main',
        'robo-phot=robotools.phot:main',
    ]
}


if __name__ == "__main__":
    setup(name='robotools',
          version='0.2.0',
          description=('Lowell Robo 31" data analysis tools.'),
          author="Michael S. P. Kelley",
          author_email="msk@astro.umd.edu",
          url="https://github.com/mkelley/robotools",
          packages=find_packages(),
          requires=['numpy', 'astropy', 'mskpy',
                    'ccdproc', 'photutils', 'sep'],
          entry_points=entry_points,
          license='BSD',
          )
