#!bash

# Initial cleanup
rm -rf dist
rm -rf build
FILE=py_diAID.pkg
if test -f "$FILE"; then
  rm py_diAID.pkg
fi
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n pydiaidinstaller python=3.8 -y
conda activate pydiaidinstaller

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_macos_gui
pip install "../../dist/pydiaid-0.0.21-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.7
pyinstaller ../pyinstaller/pydiaid.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../pydiaid/data/*.fasta dist/pydiaid/data

# Wrapping the pyinstaller folder in a .pkg package
mkdir -p dist/pydiaid/Contents/Resources
cp ../logos/alpha_logo.icns dist/pydiaid/Contents/Resources
mv dist/pydiaid_gui dist/pydiaid/Contents/MacOS
cp Info.plist dist/pydiaid/Contents
cp pydiaid_terminal dist/pydiaid/Contents/MacOS
cp ../../LICENSE.txt Resources/LICENSE.txt
cp ../logos/alpha_logo.png Resources/alpha_logo.png
chmod 777 scripts/*

pkgbuild --root dist/pydiaid --identifier de.mpg.biochem.pydiaid.app --version 0.3.0 --install-location /Applications/py_diAID.app --scripts scripts py_diAID.pkg
productbuild --distribution distribution.xml --resources Resources --package-path py_diAID.pkg dist/pydiaid_gui_installer_macos.pkg
