from setuptools import setup, find_packages

setup(
    name="mindscape",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "ruamel.yaml>=0.15.0",  # For YAML configuration file handling
        "numpy>=1.18.5,<2.0.0",  # Core numerical library
        "pandas>=1.0.1",  # Data manipulation and analysis
        "scipy>=1.9",  # Advanced numerical computations (optional but useful)
        "matplotlib>=3.3,<3.9",  # Plotting library (optional but useful)
        "tqdm>=4.0.0",  # Progress bar library (useful for tracking progress in loops)
        "click>=8.0.0",  # For CLI functionality
    ],
    entry_points={
        "console_scripts": [
            "mindscape=mindscape.cli:main",  # CLI entry point
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",  # Minimum Python version required
)