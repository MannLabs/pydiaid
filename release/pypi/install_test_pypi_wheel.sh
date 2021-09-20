conda create -n diaid_pasef_pip_test python=3.8 -y
conda activate diaid_pasef_pip_test
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple "diaid_pasef[stable]"
diaid_pasef
conda deactivate
