# setup.py
from setuptools import setup, find_packages

setup(
    name="telomerehunter",
    version="1.2.0",
    description="Tool to estimate telomere content from WGS data",
    long_description=open("telomerehunter-1.2.0.dist-info/DESCRIPTION.rst").read(),
    long_description_content_type="text/x-rst",
    author="Lars Feuerbach, Philip Ginsbach, Lina Sieverling, Christian J Bouwens (fork)",
    author_email="l.sieverling@dkfz-heidelberg.de",
    url="https://www.dkfz.de/en/applied-bioinformatics/telomerehunter/telomerehunter.html",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "PyPDF2",
        "pysam==0.20.0",
        "numpy",
    ],
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    keywords=["telomere", "content", "read", "NGS", "WGS", "tumor", "control"],
)