# Copyright (c) 2018 Dana Wyman
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from setuptools import find_packages, setup

with open("README.md", "r") as readme_file:
    LONG_DESCRIPTION = readme_file.read()

setup(
    name="talon",
    version="5.0",
    description="TALON is a Python program for identifying known and novel "
                "genes/isoforms in long read transcriptome data sets",
    author="Dana Wyman",
    author_email="dwyman@uci.edu",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    license="MIT",
    keywords="talon long read transcriptome novel gene isoform analysis",
    zip_safe=False,
    packages=find_packages('src'),
    package_dir={'': 'src'},
    package_data={'talon.post': ["r_scripts/*.R"]},
    url="https://github.com/dewyman/TALON",
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires=">=3.6",
    install_requires=[
        "pandas",
        "pysam>=0.15.4",
        "pybedtools",
        "pyfaidx"
    ],
    entry_points={
        "console_scripts": [
            'talon=talon.talon:main',
            'talon_label_reads=talon.talon_label_reads:main',
            'talon_initialize_database=talon.initialize_talon_database:main',
            'talon_filter_transcripts=talon.post.filter_talon_transcripts:main',
            'talon_abundance=talon.post.create_abundance_file_from_database:main',
            'talon_create_GTF=talon.post.create_GTF_from_database:main',
            'talon_reformat_gtf=talon.reformat_gtf:main',
            'talon_generate_report=talon.post.generate_talon_report:main',
            'talon_summarize=talon.post.summarize_datasets:main',
            'talon_fetch_reads=talon.post.get_read_annotations:main',
            'talon_get_sjs=talon.post.get_transcript_sjs:main',
            'talon_longest_end=talon.post.call_longest_ends:main'
        ]
    }
)
