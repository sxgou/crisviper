from setuptools import setup, find_packages

setup(
    name="lineage_tracer_amplicon",
    version="1.0.0",
    description="CRISPR lineage tracing amplicon analysis software",
    author="Developer",
    author_email="developer@example.com",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "biopython>=1.81",
        "pandas>=2.0",
        "numpy>=1.24",
        "scipy>=1.10",
        "matplotlib>=3.7",
        "pyyaml>=6.0",
        "edlib>=1.3.9",
        "tqdm>=4.65",
        "pysam>=0.21",
    ],
    extras_require={
        "dev": [
            "pytest>=7.4",
            "pytest-cov>=4.1",
            "black>=23.7",
            "mypy>=1.5",
            "flake8>=6.1",
        ],
    },
    entry_points={
        "console_scripts": [
            "lineage-tracer=lineage_tracer_amplicon.cli:main",
        ],
    },
    python_requires=">=3.9",
)