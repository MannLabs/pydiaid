#!bash

# Initial cleanup
rm -rf dist
rm -rf build
FILE=diAID_PASEF.pkg
if test -f "$FILE"; then
  rm diAID_PASEF.pkg
fi
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n diaid_pasefinstaller python=3.8 -y
conda activate diaid_pasefinstaller

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_macos_gui
pip install "../../dist/diaid_pasef-0.0.2-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.7
pyinstaller ../pyinstaller/diaid_pasef.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../diaid_pasef/data/*.fasta dist/diaid_pasef/data

# Wrapping the pyinstaller folder in a .pkg package
mkdir -p dist/diaid_pasef/Contents/Resources
cp ../logos/alpha_logo.icns dist/diaid_pasef/Contents/Resources
mv dist/diaid_pasef_gui dist/diaid_pasef/Contents/MacOS
cp Info.plist dist/diaid_pasef/Contents
cp diaid_pasef_terminal dist/diaid_pasef/Contents/MacOS
cp ../../LICENSE.txt Resources/LICENSE.txt
cp ../logos/alpha_logo.png Resources/alpha_logo.png
chmod 777 scripts/*

pkgbuild --root dist/diaid_pasef --identifier de.mpg.biochem.diaid_pasef.app --version 0.3.0 --install-location /Applications/diAID_PASEF.app --scripts scripts diAID_PASEF.pkg
productbuild --distribution distribution.xml --resources Resources --package-path diAID_PASEF.pkg dist/diaid_pasef_gui_installer_macos.pkg
