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
conda create -n py_diaidinstaller python=3.8 -y
conda activate py_diaidinstaller

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_macos_gui
pip install "../../dist/py_diaid-0.0.3-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.7
pyinstaller ../pyinstaller/py_diaid.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../py_diaid/data/*.fasta dist/py_diaid/data

# Wrapping the pyinstaller folder in a .pkg package
mkdir -p dist/py_diaid/Contents/Resources
cp ../logos/alpha_logo.icns dist/py_diaid/Contents/Resources
mv dist/py_diaid_gui dist/py_diaid/Contents/MacOS
cp Info.plist dist/py_diaid/Contents
cp py_diaid_terminal dist/py_diaid/Contents/MacOS
cp ../../LICENSE.txt Resources/LICENSE.txt
cp ../logos/alpha_logo.png Resources/alpha_logo.png
chmod 777 scripts/*

pkgbuild --root dist/py_diaid --identifier de.mpg.biochem.py_diaid.app --version 0.3.0 --install-location /Applications/py_diAID.app --scripts scripts py_diAID.pkg
productbuild --distribution distribution.xml --resources Resources --package-path py_diAID.pkg dist/py_diaid_gui_installer_macos.pkg
