conda create -n py_diaid_pip_test python=3.8 -y
conda activate py_diaid_pip_test
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple "py_diaid[stable]"
py_diaid
conda deactivate
