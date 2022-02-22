#!bash

# Initial cleanup
rm -rf dist
rm -rf build
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n py_diaid_installer python=3.8 -y
conda activate py_diaid_installer

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_linux_gui
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "../../dist/py_diaid-0.0.3-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.2
pyinstaller ../pyinstaller/py_diaid.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../py_diaid/data/*.fasta dist/py_diaid/data
# WARNING: this probably does not work!!!!

# Wrapping the pyinstaller folder in a .deb package
mkdir -p dist/py_diaid_gui_installer_linux/usr/local/bin
mv dist/py_diaid dist/py_diaid_gui_installer_linux/usr/local/bin/py_diaid
mkdir dist/py_diaid_gui_installer_linux/DEBIAN
cp control dist/py_diaid_gui_installer_linux/DEBIAN
dpkg-deb --build --root-owner-group dist/py_diaid_gui_installer_linux/
