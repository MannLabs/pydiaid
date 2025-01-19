![Pip installation](https://github.com/MannLabs/pydiaid/workflows/Default%20installation%20and%20tests/badge.svg)
![GUI and PyPi releases](https://github.com/MannLabs/pydiaid/workflows/Publish%20on%20PyPi%20and%20release%20on%20GitHub/badge.svg)

# py_diAID
The [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) developed py_diAID, a tool that facilitates the generation of data-independent acquisition (DIA) methods with an optimal window design. py_diAID, an abbreviation for Automated Isolation Design for Data-Independent Acquisition, is available as an open-source Python and Graphical User Interface (GUI) package. To access all the hyperlinks in this document, please view it on [GitHub](https://github.com/MannLabs/pydiaid).

* [**About**](#about)
  * [**py_diAID**](#py_diaid)
  * [**Supported Scan Modes**](#supported-scan-modes)
  * [**Key Features**](#key-features)
* [**License**](#license)
* [**Installation**](#installation)
  * [**One-click GUI**](#one-click-gui)
  * [**Pip installer**](#pip)
  * [**Developer installer**](#developer)
* [**Usage**](#usage)
  * [**GUI**](#gui)
  * [**Python and jupyter notebooks**](#python-and-jupyter-notebooks)
* [**Troubleshooting**](#troubleshooting)
* [**Citations**](#citations)
* [**How to contribute**](#how-to-contribute)
* [**Changelog**](#changelog)

---
## About

### py_diAID
py_diAID (Automated Isolation Design for Data-Independent Acquisition) is a Python package that streamlines the generation of data-independent acquisition (DIA) methods for mass spectrometry-based proteomics.

### Supported Scan Modes
The package supports several advanced scan modes. 
- __dia-PASEF__ offers comprehensive proteome coverage with a high degree of quantitative reproducibility. It uses a much larger ion beam proportion than conventional DIA methods and is particularly advantageous for studying deep proteomes from cell lines, clinical samples with regular and very low sample input, as well as for exploring post-translational modifications such as phosphorylation.
- As the successor to dia-PASEF, __synchro-PASEF__ more efficiently samples precursors leading to shorter cycle times than dia-PASEF, thereby improving quantitative reproducibility. It also allows linking fragment signals with precursor masses increasing specificity. This unites the benefits of both DDA and DIA approaches.
- __Orbitrap Astral DIA__ provides comprehensive proteome coverage with high quantitative reproducibility. Users can generate window schemes with either fixed or dynamic window widths, and the system ensures optimal window border placement within the forbidden zone.

### Key Features
py_diAID is an open-source Python package and also offers a Graphical User Interface (GUI). It was developed by the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann). py_diAID automatically generates optimal isolation windows and supports variable/dynamic isolation widths. These widths are aligned to the precursor density, enabling short acquisition cycles while covering virtually the entire m/z (-ion mobility) range. The package facilitates quality control of measured samples through precursor distribution visualization and enables evaluation of existing DIA methods. 

---
## License

py_diAID was developed by the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and is freely available with an [Apache License 2.0](LICENSE.txt). External Python packages (available in the [requirements](requirements) folder) have their own licenses, which can be consulted on their respective websites.

---
## Installation

py_diAID can be installed and used on the Windows operating system.
There are three different types of installation possible:

* [**One-click GUI installer:**](#one-click-gui) Choose this installation if you only want the GUI and/or keep things simple.
* [**Pip installer:**](#pip) Choose this installation if you want to use py_diAID as a Python package in an existing Python 3.8 environment (e.g. a Jupyter notebook). If needed, the GUI can be installed with pip.
* [**Developer installer:**](#developer) Choose this installation if you are familiar with CLI tools, [conda](https://docs.conda.io/en/latest/), and Python. This installation allows access to all available features of py_diAID and even allows to modify its source code directly. Generally, the developer version of py_diAID outperforms the precompiled versions making this the installation of choice for high-throughput experiments.

### One-click GUI

The GUI of py_diAID is a stand-alone tool that requires no knowledge of Python or CLI tools. Click on the link below to download the latest release for:

* [**Windows**](https://github.com/MannLabs/pydiaid/releases/latest/download/pydiaid_gui_installer_windows.exe)

Older releases remain available on the [release page](https://github.com/MannLabs/pydiaid/releases), but no backward compatibility is guaranteed.

**IMPORTANT: Please refer to the [GUI manual](https://github.com/MannLabs/pydiaid/blob/development/pydiaid/docs/manual.pdf) for detailed instructions on installing and using the stand-alone py_diAID GUI.**


### Pip

py_diAID can be installed in an existing Python 3.8 environment with a single `bash` command. *This `bash` command can also be run directly from within a Jupyter notebook by prepending it with a `!`*:

```bash
pip install pydiaid
```

Installing py_diAID like this avoids conflicts when integrating it in other tools, as this does not enforce strict versioning of dependencies. However, if new versions of dependencies are released, they are not guaranteed to be fully compatible with py_diAID. While this should only occur in rare cases where dependencies are not backward compatible, you can always force py_diAID to use dependency versions that are known to be compatible with:

```bash
pip install "pydiaid[stable]"
```

NOTE: You might need to run `pip install pip==21.0` before installing py_diAID like this. Also, note the double quotes `"`.

For those who are adventurous, it is also possible to directly install any branch (e.g. `@development`) with any extras (e.g. `#egg=pydiaid[stable,development-stable]`) from GitHub with e.g.

```bash
pip install "git+https://github.com/MannLabs/pydiaid.git@development#egg=pydiaid[stable,development-stable]"
```

### Developer

py_diAID can also be installed in editable (i.e. developer) mode with a few `bash` commands. This allows to fully customize the software and even modify the source code to your specific needs. When an editable Python package is installed, its source code is stored in a transparent location of your choice. While optional, it is advised to first (create and) navigate to e.g. a general software folder:

```bash
mkdir ~/folder/where/to/install/software
cd ~/folder/where/to/install/software
```

***The following commands assume you do not perform any additional `cd` commands anymore***.

Next, download the py_diAID repository from GitHub either directly or with a `git` command. This creates a new py_diAID subfolder in your current directory.

```bash
git clone https://github.com/MannLabs/pydiaid.git
```

For any Python package, it is highly recommended to use a separate [conda virtual environment](https://docs.conda.io/en/latest/), as otherwise *dependency conflicts can occur with already existing packages*.

```bash
conda create --name pydiaid python=3.8 -y
conda activate pydiaid
```

Finally, py_diAID and all its [dependencies](requirements) need to be installed. To take advantage of all features and allow development (with the `-e` flag), this is best done by also installing the [development dependencies](requirements/requirements_development.txt) instead of only the [core dependencies](requirements/requirements.txt):

```bash
pip install -e "./pydiaid[development]"
```

By default this installs loose dependancies (no explicit versioning), although it is also possible to use stable dependencies (e.g. `pip install -e "./pydiaid[stable,development-stable]"`).

***By using the editable flag `-e`, all modifications to the [py_diAID source code folder](pydiaid) are directly reflected when running py_diAID. Note that the py_diAID folder cannot be moved and/or renamed if an editable version is installed.***

---
## Usage

There are two ways to use py_diAID:

* [**GUI**](#gui)
* [**Python**](#python-and-jupyter-notebooks)

NOTE: The first time you use a fresh installation of py_diAID, it is often relatively slow because some functions might still need compilation on your local operating system and architecture. Subsequent executions should be a lot faster.

### GUI

If the GUI was not installed through a one-click GUI installer, it could be activated with the following `bash` command:

```bash
pydiaid gui
```

Note that this needs to be prepended with a `!` when you want to run this from within a Jupyter notebook. When the command is run directly from the command-line, make sure you use the right environment (activate it with e.g. `conda activate pydiaid` or set an alias to the binary executable (can be obtained with `where pydiaid` or `which pydiaid`)).

**IMPORTANT: Please refer to the [GUI manual](https://github.com/MannLabs/pydiaid/blob/development/pydiaid/docs/manual.pdf) for detailed instructions on installing, troubleshooting, and using the stand-alone py_diAID GUI.**

### Python and Jupyter notebooks

py_diAID can be imported as a Python package into any Python script or notebook with the command `import pydiaid`.

An ‘nbs’ folder in the GitHub repository contains Jupyter Notebooks as tutorials on using py_diAID as a Python package.  

---
## Troubleshooting

In case of issues, check out the following links:

* [FAQ](https://github.com/MannLabs/pydiaid#faq): This section provides answers to issues of general interest.
* [Issues](https://github.com/MannLabs/pydiaid/issues): Try a few different search terms to find out if a similar problem has been encountered before.
* [Discussions](https://github.com/MannLabs/pydiaid/discussions): Check if your problem or feature request has been discussed earlier.

---
## FAQ
- Where to find test libraries? The py_diAID package includes test libraries for quick workflow testing. These can be found at: py_diAID installation directory\pydiaid\diapasef\static\AlphaPept_results.csv for dia-PASEF and py_diAID installation directory\pydiaid\synchropasef\static\evidence_MaxQuant_270223.txt for synchro-PASEF.
- What is the best input for py_diAID method generation? In general, the best input for py_diAID method generation is dda-PASEF acquired with a wide ion mobility range, for instance, from 0.6-1.6. It provides a complete and unbiased view of the precursor cloud in m/z and the ion mobility plane. In contrast, the data collected with dia-PASEF will present a precursor cloud that is influenced by the position of their isolation windows. The dda-PASEF runs may be an analysis of a single-run representative of the study, or a fractionated peptide library. Both these approaches have yielded comparable isolation window schemes. Regardless of the strategy used, the most critical aspect is a precise ion mobility calibration.
- Using DIA-NN results as input for py_diAID: We have now included an option to load DIA-NN libraries and single-runs.
- How to specify multiple PTMs? The initial versions of py_diAID could only process one PTM or string input at a time. We have now updated it to allow for filtering of the input library for multiple PTMs. To do this, all PTMs need to be specified in a list of strings, for instance ["STY", "GlyGly"].

---
## Citations

Check out the [optimal dia-PASEF](https://doi.org/10.1016/j.mcpro.2022.100279), [synchro-PASEF](https://doi.org/10.1016/j.mcpro.2022.100489) and [PASEF workflows and py_diAID](https://doi.org/10.1038/s41596-024-01104-w) publications.

---
## How to contribute

If you like this software, you can give us a [star](https://github.com/MannLabs/pydiaid/stargazers) to boost our visibility! All direct contributions are also welcome. Feel free to post a new [issue](https://github.com/MannLabs/pydiaid/issues) or clone the repository and create a [pull request](https://github.com/MannLabs/pydiaid/pulls) with a new branch. For even more interactive participation, check out the [discussions](https://github.com/MannLabs/pydiaid/discussions) and [the Contributors License Agreement](misc/CLA.md).


---
## Changelog

See the [HISTORY.md](HISTORY.md) for a complete overview of the changes made in each version.
