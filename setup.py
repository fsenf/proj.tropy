

# from setuptools import setup
from numpy.distutils.core import setup, Extension

#f2py -c -m --f90flags='-O3 -fPIC' _f90_bit_conversion_altair bit_conversion.f90
e1 = Extension( '_f90_bit_conversion', 
                sources = ['tropy/io_tools/bit_conversion.f90',],
                extra_f90_compile_args = ['-O3 -fPIC'], )

setup(name='tropy',
      version='0.1',
      description='Collection of IO and Analysis Tools for TROPOS work',
      author='Fabian Senf',
      author_email='senf@tropos.de',
      license='GPL',
      packages=['tropy', 'tropy.analysis_tools', 'tropy.io_tools', 'tropy.plotting_tools', 'tropy.l15_msevi'],
      ext_package = 'tropy.io_tools',
      ext_modules = [e1, ], 
      zip_safe=False)

