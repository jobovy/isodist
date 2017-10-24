# isodist

## AUTHOR

Jo Bovy - bovy at astro dot utoronto dot ca

If you find this code useful in your research, please let me know and
consider citing this package. Thanks!

## INSTALLATION

Standard python setup.py build/install

## DEPENDENCIES

This package requires [NumPy](http://numpy.scipy.org/), [Scipy](http://www.scipy.org/), and [Matplotlib](http://matplotlib.sourceforge.net/).

## DATA DIRECTORY

``isodist`` requires sets of isochrones, which can be downloaded from various isochrone libraries' websites. These should be contained in directories under a general data directory for the ``isodist`` package. The general directory should be referenced by an environment variable ``ISODIST_DATA``, such that it can be accessed as ``$ISODIST_DATA`` (for example, ``ls $ISODIST_DATA`` should give the contents of this directory. Environment variables can be defined as ``export ISODIST_DATA=/path/to/isodist/data/directory`` in bash-style shells and ``setenv ISODIST_DATA /path/to/isodist/data/directory`` in C-style shells. Individual isochrones should be stored in sub-directories of the general directory. For example, the ``parsec-2mass-spitzer-wise`` directory should contain ``PARSEC`` isochrones with the ``2mass-spitzer-wise`` filter set (these names can then be used in the ``isodist`` code and will be found).

Please contact the author of this package for help with obtaining the isochrones.