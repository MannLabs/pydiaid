conda create -n diaid_pasef python=3.8 -y
conda activate diaid_pasef
pip install -e '../.[stable,development-stable]'
diaid_pasef
conda deactivate
