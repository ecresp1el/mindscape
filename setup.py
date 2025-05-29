from setuptools import setup, find_packages

setup(
    name="mindscape",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "click",  # Add other dependencies here
        "ruamel.yaml",
        "numpy",
        "pandas",
    ],
    entry_points={
        "console_scripts": [
            "mindscape=mindscape.cli:main",  # CLI entry point
        ],
    },
)