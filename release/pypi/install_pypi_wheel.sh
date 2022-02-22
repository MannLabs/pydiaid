conda create -n pydiaid_pip_test python=3.8 -y
conda activate pydiaid_pip_test
pip install "pydiaid[stable]"
pydiaid
conda deactivate
