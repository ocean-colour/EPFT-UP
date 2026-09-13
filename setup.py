
# Standard imports
import glob, os
from setuptools import setup, find_packages


# Begin setup
setup_keywords = dict()
setup_keywords['name'] = 'epft_up'
setup_keywords['description'] = 'Empirical Phytoplankton Functional Types with Uncertainty and Provenance'
setup_keywords['author'] = 'J. Xavier Prochaska'
setup_keywords['author_email'] = 'jxp@ucsc.edu'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/ocean-colour/EPFT-UP'
setup_keywords['version'] = '0.0.dev0'
# Use README.md as long_description.
setup_keywords['long_description'] = ''
if os.path.exists('README.md'):
    with open('README.md') as readme:
        setup_keywords['long_description'] = readme.read()
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['python_requires'] = '>=3.12'
setup_keywords['install_requires'] = [
    'numpy', 'scipy', 'pandas', 'pyarrow', 'matplotlib', 'seaborn',
    'scikit-learn', 'tqdm', 'IPython', 'pytest',
    # Oceanography / gridded data I/O
    'xarray', 'h5netcdf', 'h5py', 'netcdf4', 'cftime',
    # NASA Earthdata access (MOANA track; skip-guarded)
    'earthaccess']
# Sibling packages (e.g. ocpy) are not on PyPI; install them from
# source / GitHub via requirements.txt (git+https://github.com/ocean-colour/...).
setup_keywords['zip_safe'] = False
setup_keywords['packages'] = find_packages()
# Vendored NASA MOANA reference tables (sha256-pinned; see epft_up/data/moana/README.md)
setup_keywords['package_data'] = {'epft_up': ['data/moana/*']}
setup_keywords['include_package_data'] = True

if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
                                 if not os.path.basename(fname).endswith('.rst')]

setup(**setup_keywords)
