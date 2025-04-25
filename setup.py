"""
Setup script for installing the Evo2_Clinical package.
"""

from setuptools import setup, find_packages
import os
import re

# Read the version from __init__.py
with open(os.path.join('evo2_clinical', '__init__.py'), 'r') as f:
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", f.read(), re.M)
    if version_match:
        version = version_match.group(1)
    else:
        version = '0.1.0'  # Default version if not found

# Read the long description from README.md
with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='evo2_clinical',
    version=version,
    description='Computational framework for investigating endothelial influence on pulmonary disease',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Shabab Khan, Jeffrey Man, Aninda Dibya Saha',
    author_email='shababkhan@example.com',  # Replace with actual email
    url='https://github.com/shababkhan/Evo2_Clinical',  # Replace with actual URL
    packages=find_packages(),
    install_requires=[
        'pandas>=1.3.0',
        'numpy>=1.20.0',
        'pyyaml>=5.1',
        # Add other dependencies as needed
    ],
    entry_points={
        'console_scripts': [
            'evo2-clinical=evo2_clinical.cli:main',
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.7',
)