from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="transnet",
    version="0.1.0",
    author="Matthias Anagho-Mattanovich",
    author_email="matthias.mattanovich@sund.ku.dk",
    description="A package for trans-omics data integration and network analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mmattano/transnet",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.2.0",
        "networkx>=2.5",
        "matplotlib>=3.3.0",
        "scipy>=1.6.0",
        "scikit-learn>=0.24.0",
        "statsmodels>=0.12.0",
        "bioservices>=1.9.0",
        "mygene>=3.2.0",
        "requests>=2.25.0",
        "zeep>=4.0.0",
        "biopython>=1.78",
        "pyensembl>=2.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.10",
            "black>=21.0",
            "flake8>=3.8",
            "sphinx>=4.0",
            "nbsphinx>=0.8",
        ],
        "viz": [
            "seaborn>=0.11.0",
            "plotly>=5.0.0",
        ]
    }
)