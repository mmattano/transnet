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
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "pandas",
        "networkx",
        "matplotlib",
        "scipy",
        "bioservices",
        "mygene",
        "python-louvain",  # community detection
        "requests",
        "biopython"
    ],
    extras_require={
        "dev": [
            "pytest",
            "pytest-cov",
            "black",
            "flake8"
        ]
    }
)
