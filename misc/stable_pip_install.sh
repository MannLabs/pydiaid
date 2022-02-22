conda create -n pydiaid python=3.8 -y
conda activate pydiaid
pip install -e '../.[stable,development-stable]'
pydiaid
conda deactivate
