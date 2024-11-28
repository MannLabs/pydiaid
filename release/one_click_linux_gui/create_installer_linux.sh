#!bash

# Initial cleanup
rm -rf dist
rm -rf build
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n pydiaid_installer python=3.8 -y
conda activate pydiaid_installer

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_linux_gui
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "../../dist/pydiaid-0.0.28-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.10
pyinstaller ../pyinstaller/pydiaid.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../pydiaid/data/*.fasta dist/pydiaid/data
# WARNING: this probably does not work!!!!

# Wrapping the pyinstaller folder in a .deb package
mkdir -p dist/pydiaid_gui_installer_linux/usr/local/bin
mv dist/pydiaid dist/pydiaid_gui_installer_linux/usr/local/bin/pydiaid
mkdir dist/pydiaid_gui_installer_linux/DEBIAN
cp control dist/pydiaid_gui_installer_linux/DEBIAN
dpkg-deb --build --root-owner-group dist/pydiaid_gui_installer_linux/
