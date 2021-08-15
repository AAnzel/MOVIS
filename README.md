![movis_logo](./Source/images/movis_logo_transparent.png)

[![Build status](https://github.com/AAnzel/MOVIS/actions/workflows/main.yml/badge.svg)](https://github.com/AAnzel/MOVIS/actions/workflows/main.yml)

---
# MOVIS

**Exploratory data analysis and visualization tool for time-series multi-omics data sets.**


## Manuscript

This tool is created for the following paper:

**MOVIS: A Multi-Omics Software Solution for Multi-modal Time-Series Clustering and Embedding Tasks**
Aleksandar Anžel, Dominik Heider, and Georges Hattab

**Paper badge placeholder, link to the PDF placeholder**

Please cite the paper as:
``` Bibtex citation placeholder
```

---
Abstract:

> **Motivation**:
Thanks to recent advances in sequencing and computational technologies, many researchers with biological and/or medical backgrounds are now producing multiple data sets with an embedded temporal dimension. Motivated by the exploration of these time-series data to uncover anomalies, find patterns and gain insights from the data, this step has remained a complex challenge that often requires expertise in computer science and bioinformatics.

>**Results**:
Here we present MOVIS, a modular time-series multi-omics exploration tool with a user-friendly web interface that integrates, analyzes, and visualizes different time-series omics data sets simultaneously. The resulting visualizations are task-specific, reproducible, and publication-ready. MOVIS is built on open-source software and is easily extendable to accommodate different analytical tasks.

>**Availability**:
We provide MOVIS as a web service at INSERTURL. Source code, help and documentation can be found at https://github.com/AAnzel/MOVIS.

>**Contact**: [aleksandar.anzel@uni-marburg.de](mailto:aleksandar.anzel@uni-marburg.de)

>**Supplementary information**: Supplementary data are available at [Bioinformatics]()
online.

**Paper image placeholder**

## Dependancy

The code is written in Python 3.8.10 and tested on Linux with the following libraries installed:

|Library|Version|
|---|---|
|altair|4.1.0|
|biopython|1.78|
|gensim|4.0.1|
|numpy|1.20.3|
|pandas|1.2.5|
|scikit-learn|0.24.2|
|scipy|1.6.2|
|streamlit|0.83.0|


## Data
The data used in the Example 1 comes from the following paper:

> **Integration of time-series meta-omics data reveals how microbial ecosystems respond to disturbance**, Herold, M., Martínez Arbas, S., Narayanasamy, S. et al. Nat Commun 11, 5281(2020).
https://doi.org/10.1038/s41467-020-19006-2.

It is stored at [Data/cached/example_1/](./Data/cached/example_1) in either a raw format or as a [pickle](https://docs.python.org/3/library/pickle.html) object.

## Code
|Script|Description|
|---|---|
|[Source/](./Source/)|contains all scripts necessary to run the tool.
|[Source/main.py](./Source/main.py)|contains the code that builds the main layout and connects all pages.
|[Source/home.py](./Source/home.py)|contains the code that builds the home page.
|[Source/example_1.py](./Source/example_1.py)|contains the code that builds the example page.
|[Source/upload.py](./Source/upload.py)|contains the code that builds the upload page.
|[Source/common.py](./Source/common.py)|contains the code with functions shared by all pages.
|[Source/visualize.py](./Source/visualize.py)|contains the code with functions that create various visualizations present in this tool.

## Getting started
Some text on what a user can expect

## Installation & Running
How to run the tool


## License

Licensed under the GNU General Public License, Version 3.0 ([LICENSE](./LICENSE) or https://www.gnu.org/licenses/gpl-3.0.en.html)

### Contribution

Any contribution intentionally submitted for inclusion in the work by you, shall be licensed under the GNU GPLv3.
