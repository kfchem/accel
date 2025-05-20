from pathlib import Path

from setuptools import find_packages, setup

with Path("accel").joinpath("__init__.py").open("r") as f:
    _version = f.readline().split()[2].replace("'", "").replace('"', "")

setup(
    name="accel",
    version=_version,
    description="ACCeL is a Python package for managing and filtering conformers in computational chemistry.",
    long_description=open("README.rst", encoding="utf-8").read(),
    long_description_content_type="text/x-rst",
    author="Keisuke Fukaya",
    author_email="kfukaya@pu-toyama.ac.jp",
    url="https://github.com/kfchem/accel",
    license="MIT License",
    install_requires=["numpy"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    packages=find_packages(exclude=("tests", "docs", "application")),
    include_package_data=True,
)
