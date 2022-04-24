from pathlib import Path

from setuptools import find_packages, setup

with open("README.rst") as f:
    readme = f.read()

with open("LICENSE") as f:
    license = f.read()

with Path("accel").joinpath("__init__.py").open("r") as f:
    _version = f.readline().split()[2]

setup(
    name="accel",
    version=_version,
    description="python package for manimuplating multiple conformers",
    long_description=readme,
    author="Keisuke Fukaya",
    author_email="kfukaya@pu-toyama.ac.jp",
    url="httpx://github.com/kfchem/accel",
    license=license,
    install_requires=["numpy"],
    packages=find_packages(exclude=("tests", "docs", "application")),
)
"""
change __version__
python setup.py bdist_wheel
python setup.py clean --all
"""
