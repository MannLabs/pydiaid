![Pip installation](https://github.com/MannLabs/pydiaid/workflows/Default%20installation%20and%20tests/badge.svg)
![GUI and PyPi releases](https://github.com/MannLabs/pydiaid/workflows/Publish%20on%20PyPi%20and%20release%20on%20GitHub/badge.svg)

# py_diAID
The [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) developed py_diAID generating dia-PASEF methods with an optimal window design. py_diAID stands for an open-source Python package for Data-Independent Acquisition with an Automated Isolation Design. To enable all hyperlinks in this document, please view it on [GitHub](https://github.com/MannLabs/pydiaid).

* [**About**](#about)
* [**License**](#license)
* [**Installation**](#installation)
  * [**One-click GUI**](#one-click-gui)
  * [**Pip installer**](#pip)
  * [**Developer installer**](#developer)
* [**Usage**](#usage)
  * [**GUI**](#gui)
  * [**CLI**](#cli)
  * [**Python and jupyter notebooks**](#python-and-jupyter-notebooks)
* [**Troubleshooting**](#troubleshooting)
* [**Citations**](#citations)
* [**How to contribute**](#how-to-contribute)
* [**Changelog**](#changelog)

---
## About

Data-independent acquisition coupled with parallel accumulation – serial fragmentation (dia-PASEF) has gained increasing attention from proteomics researchers over the last years. dia-PASEF offers comprehensive proteome coverage, a high degree of reproducibility, and quantitative accuracy while using a much higher ion beam proportion than conventional DIA methods. Previous tools generated dia-PASEF methods with equidistant isolation widths and required a fitting of the window design to the precursor density cloud by hand. We present py_diAID, a Python-based package for Data-Independent Acquisition offering an Automated Isolation Design. py_diAID generates optimal dia-PASEF methods with variable isolation widths adjusted to the precursor density in m/z and automatically, optimally placed in the m/z – ion mobility (IM) plane. Variable isolation widths enable short acquisition cycles while covering the complete m/z-IM-range essentially. We found dia-PASEF methods, generated with py_diAID, beneficial for optimizing proteomics workflows based on cell lines (HeLa) or clinical samples such as CSF and Plasma and for studying post-translational modifications such as phosphorylation.

We offer py_diAID as a Python module, command-line tool, and graphical user interface on all major operating systems under an Apache 2.0 license. py_diAID generates dia-PASEF methods with an optimal window design. It also allows for quality control of the precursors’ distribution of a dataset in the m/z-ion mobility plane and evaluating the suitability of already existing dia-PASEF methods for the individual experiment.

py_diAID is an open-source Python package from the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann).

---
## License

py_diAID was developed by the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and is freely available with an [Apache License 2.0](LICENSE.txt). External Python packages (available in the [requirements](requirements) folder) have their own licenses, which can be consulted on their respective websites.

---
## Installation

py_diAID can be installed and used on all major operating systems (Windows, macOS, and Linux).
There are three different types of installation possible:

* [**One-click GUI installer:**](#one-click-gui) Choose this installation if you only want the GUI and/or keep things simple.
* [**Pip installer:**](#pip) Choose this installation if you want to use py_diAID as a Python package in an existing Python 3.8 environment (e.g. a Jupyter notebook). If needed, the GUI and CLI can be installed with pip.
* [**Developer installer:**](#developer) Choose this installation if you are familiar with CLI tools, [conda](https://docs.conda.io/en/latest/), and Python. This installation allows access to all available features of py_diAID and even allows to modify its source code directly. Generally, the developer version of py_diAID outperforms the precompiled versions making this the installation of choice for high-throughput experiments.

### One-click GUI

The GUI of py_diAID is a stand-alone tool that requires no knowledge of Python or CLI tools. Click on one of the links below to download the latest release for:

* [**Windows**](https://github.com/MannLabs/pydiaid/releases/latest/download/pydiaid_gui_installer_windows.exe)
* [**macOS**](https://github.com/MannLabs/pydiaid/releases/latest/download/pydiaid_gui_installer_macos.pkg)
* [**Linux**](https://github.com/MannLabs/pydiaid/releases/latest/download/pydiaid_gui_installer_linux.deb)

Older releases remain available on the [release page](https://github.com/MannLabs/pydiaid/releases), but no backward compatibility is guaranteed.

**IMPORTANT: Please refer to the [GUI manual](https://github.com/MannLabs/pydiaid/blob/development/pydiaid/docs/manual.pdf) for detailed instructions on installing, troubleshooting, and using the stand-alone py_diAID GUI.**

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

There are three ways to use py_diAID:

* [**GUI**](#gui)
* [**CLI**](#cli)
* [**Python**](#python-and-jupyter-notebooks)

NOTE: The first time you use a fresh installation of py_diAID, it is often relatively slow because some functions might still need compilation on your local operating system and architecture. Subsequent use should be a lot faster.

### GUI

If the GUI was not installed through a one-click GUI installer, it could be activated with the following `bash` command:

```bash
pydiaid gui
```

Note that this needs to be prepended with a `!` when you want to run this from within a Jupyter notebook. When the command is run directly from the command-line, make sure you use the right environment (activate it with e.g. `conda activate pydiaid` or set an alias to the binary executable (can be obtained with `where pydiaid` or `which pydiaid`)).

**IMPORTANT: Please refer to the [GUI manual](https://github.com/MannLabs/pydiaid/blob/development/pydiaid/docs/manual.pdf) for detailed instructions on installing, troubleshooting, and using the stand-alone py_diAID GUI.**

### CLI

The CLI can be run with the following command (after activating the `conda` environment with `conda activate pydiaid` or if an alias was set to the pydiaid executable):

```bash
pydiaid -h
```

It is possible to get help with each function and its (required) parameters by using the `-h` flag. For instance, the command ```pydiaid optimize -h``` will produce the following output:

```
******************
* py_diAID 0.0.9 *
******************
Usage: pydiaid optimize [OPTIONS]

  Optimize a dia-PASEF method.

Options:
  -p TEXT     Parameter file (check out
              d:\pydiaid\pydiaid\lib\default_parameters.json for an example)
              [required]
  -h, --help  Show this message and exit.
```

py_diAID provides several options:
- charge: Evaluate a dia-PASEF method for multiply charged precursors.
- create: Create a specific dia-PASEF method.
- evaluate: Evaluate a dia-PASEF method.
- gui: Start graphical user interface.
- optimize: Optimize a dia-PASEF method.

All options can be executed with ```pydiaid [option] -p [Text]```. The parameters are saved in a .json parameter file and have to be adjusted in this file. For instance, the command ```pydiaid optimize -p "d:\pydiaid\pydiaid\lib\default_parameters.json"``` will execute one complete optimization process. py_diAID will create a folder at the location specified in the .json parameter file with all generated information and give the following result in the terminal window:

```
******************
* py_diAID 0.0.9 *
******************
Using parameter file d:\pydiaid\pydiaid\lib\default_parameters.json
{'precursors within m/z-range [%]': 97.59}
RUN WITH: [0.7435820751492209, 0.9789579174732773, 1.73455196349983, 1.5945652845606708] | RESULT: 10823.0
RUN WITH: [0.7468616079634967, 0.8864876711590428, 1.6545030790686948, 1.6033036772717069] | RESULT: 11168.0
RUN WITH: [0.7346248237666703, 0.8734794138726382, 1.5853169071713897, 1.6400014924345845] | RESULT: 11138.0
RUN WITH: [0.8118354704659128, 0.9634373141944189, 1.7636724593721786, 1.5662199385632] | RESULT: 10794.0
RUN WITH: [0.8150009805616482, 1.0092137143383832, 1.274205264719582, 1.5310500464391847] | RESULT: 11918.0
########
BEST RESULT
INPUT: [0.8150009805616482, 1.0092137143383832, 1.274205264719582, 1.5310500464391847]
OUTPUT: 11918.0
########
```

### Python and Jupyter notebooks

py_diAID can be imported as a Python package into any Python script or notebook with the command `import pydiaid`.

An ‘nbs’ folder in the GitHub repository contains several Jupyter Notebooks as tutorials on using py_diAID as a Python package.  

---
## Troubleshooting

In case of issues, check out the following:

* [Issues](https://github.com/MannLabs/pydiaid/issues): Try a few different search terms to find out if a similar problem has been encountered before
* [Discussions](https://github.com/MannLabs/pydiaid/discussions): Check if your problem or feature request has been discussed earlier.

---
## Citations

A manuscript is currently in preparation.

---
## How to contribute

If you like this software, you can give us a [star](https://github.com/MannLabs/pydiaid/stargazers) to boost our visibility! All direct contributions are also welcome. Feel free to post a new [issue](https://github.com/MannLabs/pydiaid/issues) or clone the repository and create a [pull request](https://github.com/MannLabs/pydiaid/pulls) with a new branch. For even more interactive participation, check out the [discussions](https://github.com/MannLabs/pydiaid/discussions) and [the Contributors License Agreement](misc/CLA.md).

---
## Changelog

See the [HISTORY.md](HISTORY.md) for a complete overview of the changes made in each version.
