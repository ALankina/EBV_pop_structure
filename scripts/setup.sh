
# make a pyenv
/share/apps/python-3.9.0-shared/bin/python3 -m venv /SAN/breuerlab/pathseq1/oc/envs/pybio

# activate the env
source /SAN/breuerlab/pathseq1/oc/envs/pybio/bin/activate

# install required things
python3 -m pip install Bio pandas numpy pycountry-convert