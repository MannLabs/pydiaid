#!bash

# Initial cleanup
rm -rf dist
rm -rf build
cd ../..
rm -rf dist
rm -rf build

conda init
conda config --remove channels defaults
conda config --add channels conda-forge

# Creating a conda environment
conda create -n pydiaid_installer python=3.8 -y
conda activate pydiaid_installer

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_windows_gui
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "../../dist/pydiaid-0.0.37-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.10

conda install -y scikit-optimize=0.9.0 numpy=1.23.5 openblas=0.3.23 scipy=1.10.1 pandas=1.3.4 matplotlib=3.6.0

pyinstaller ../pyinstaller/pydiaid.spec -y --clean
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../pydiaid/data/*.fasta dist/pydiaid/data

# Wrapping the pyinstaller folder in a .exe package
"C:\Program Files (x86)\Inno Setup 6\ISCC.exe" pydiaid_innoinstaller.iss
# WARNING: this assumes a static location for innosetup
