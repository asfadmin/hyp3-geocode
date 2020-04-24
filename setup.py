import os

from setuptools import find_packages, setup

_HERE = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(_HERE, 'README.md'), 'r') as f:
    long_desc = f.read()

setup(
    name='hyp3_geocode',
    use_scm_version=True,
    description='HyP3 plugin to create geocoded products from Sentinel-1 granules',
    long_description=long_desc,
    long_description_content_type='text/markdown',

    url='https://github.com/asfadmin/hyp3-geocode',

    author='ASF APD/Tools Team',
    author_email='uaf-asf-apd@alaska.edu',

    license='BSD',
    include_package_data=True,

    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        ],

    python_requires='~=3.6',

    install_requires=[
        'hyp3lib',
        'hyp3proclib',
        'importlib_metadata',
        'numpy',
    ],

    extras_require={
        'develop': [
            'pytest',
            'pytest-cov',
            'pytest-console-scripts',
            'tox',
            'tox-conda',
        ]
    },

    packages=find_packages(),

    entry_points={'console_scripts': [
            'hyp3_geocode = hyp3_geocode.__main__:main',
            'geocode_sentinel.py = hyp3_geocode.sentinel:main',
        ]
    },

    zip_safe=False,
)
