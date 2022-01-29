from setuptools import setup

setup(name='nmrformd',
      version='v0.0.3',
      description='Calculate NMR relaxation time from molecular dynamics trajectory file',
      long_description=open('pypi/README.rst').read(),
      long_description_content_type='text/x-rst',
      url='https://github.com/simongravelle/nmrformd',
      download_url='https://github.com/simongravelle/nmrformd/archive/refs/tags/v0.0.3.tar.gz',
      author='Simon Gravelle',
      author_email='simon.gravelle@live.fr',
      license='GNU GENERAL PUBLIC LICENSE',
      packages=['nmrformd'],
      zip_safe=False,
         install_requires=[
       "mdanalysis",
       "pytest",
       "numpy",
       "scipy",
      ]      
      )
