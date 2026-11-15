#!/bin/bash


[ -z "$1" ] && echo "usage $0 /path/to/conda" && exit

export PATH="/home/losafe/users/atarasen/miniconda/miniconda3/bin:$PATH"

[ ! -f $1/bin/activate ] && conda create -y -p $1 python=2.7


. activate $(readlink -f $1)

if [ -z "$CONDA_PREFIX" ] ; then
	echo "failed to activate $1" >&2
	exit 1
fi


conda install -y numpy=1.11 cython krb5 pyproj basemap scipy shapely netcdf4=1.2.2 matplotlib=2.0.2 requests Pillow pandas GDAL tensorflow xerces-c=3.1.2 libnetcdf

conda install -c https://conda.binstar.org/jgomezdans proj4
conda install -c  https://conda.binstar.org/jgomezdans hdf4
conda install -c https://conda.binstar.org/jgomezdans gdal
conda install geos=3.4.2 opencv

# Install cerbere
pip install git+https://git.cersat.fr/cerbere/cerbere.git@develop
pip install git+https://git.cersat.fr/cerbere/cerinterp.git@master
pip install git+https://git.cersat.fr/cerbere/ceraux.git@master
pip install git+https://git.cersat.fr/cerbere/cerform.git@master

#pip install git+https://git.cersat.fr/cerbere/cerbere.git@feature/python-hdf4


# Install sar package
pip install http://git.oceandatalab.com/downloads/sar-0.0.1.tar.gz 

#if ! pip install git+https://git.oceandatalab.com/gilles.guitton/sar.git@develop; then
#    pip install http://git.oceandatalab.com/downloads/sar-0.0.1.tar.gz 
#fi

