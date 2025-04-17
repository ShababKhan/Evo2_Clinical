from setuptools import setup, find_packages
import os

# Function to read requirements from requirements.txt
def parse_requirements(filename):
    """Load requirements from a pip requirements file."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    # Filter out comments and empty lines
    requirements = [line.strip() for line in lines if line.strip() and not line.startswith('#')]
    return requirements

# Get the long description from the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Define package requirements
requirements = parse_requirements('requirements.txt')

setup(
    name='evo2_clinical',
    version='0.1.0',
    author='Shabab Khan',
    author_email='your.email@example.com',  # Replace with actual email
    description='Python package for clinical data analysis using Evo2 pipeline',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/your_username/Evo2_Clinical',  # Replace with actual URL
    packages=find_packages(exclude=['tests*', 'examples*', 'notebooks*']),
    include_package_data=True,  # Include non-code files specified in MANIFEST.in
    install_requires=requirements,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',  # Choose an appropriate license
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.8',
    entry_points={
        'console_scripts': [
            'evo2_pipeline=evo2_pipeline:main',
            'evo2_app=app:main', # Assuming app.py has a main function for streamlit
        ],
    },
    # Add data files if needed, e.g., configuration templates
    # data_files=[('config', ['config/default.yaml'])],
)

# Note: The original setup.py contained os.system calls to install packages
# and import checks. This is not standard practice for setup.py.
# Dependencies should be listed in install_requires or requirements.txt.
# The environment setup should be handled by pip during installation.
# The import checks can be part of the test suite.
