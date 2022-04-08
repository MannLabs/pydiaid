conda create -n pydiaid python=3.8 -y
conda activate pydiaid
pip install -e '../.[development]'
pydiaid
conda deactivate
