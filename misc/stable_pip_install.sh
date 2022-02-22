conda create -n py_diaid python=3.8 -y
conda activate py_diaid
pip install -e '../.[stable,development-stable]'
py_diaid
conda deactivate
