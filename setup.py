from setuptools import setup

setup(name='nmrformd',
      version='0.0.0',
      description='Calculate NMR relaxation time from molecular dynamics trajectory file',
      long_description=open('README').read(),
      url='https://github.com/simongravelle/nmrformd',
      download_url='https://github.com/simongravelle/nmrformd/archive/refs/tags/v0.0.0.tar.gz',
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
