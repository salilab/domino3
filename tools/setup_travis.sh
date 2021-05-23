#!/bin/bash -e

# Set up an environment to run tests under Travis CI (see ../.travis.yml)

if [ $# -ne 1 ]; then
  echo "Usage: $0 python_version"
  exit 1
fi

cur_dir=$(pwd)
python_version=$1
temp_dir=$(mktemp -d)

cd ${temp_dir}

conda config --remove channels defaults  # get conda-forge, not main, packages
conda create --yes -q -n python${python_version} -c salilab -c conda-forge python=${python_version} pip nose imp-nightly gxx_linux-64 eigen swig cmake
eval "$(conda shell.bash hook)"
conda activate python${python_version}
pip install coverage

# IMP tests use sys.argv[0] to determine their location, which won't work if
# we use nosetests, so add a workaround
ln -sf $(which nosetests) ${cur_dir}/test/

cd ${cur_dir}

rm -rf ${temp_dir}
