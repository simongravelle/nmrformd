from setuptools import setup

setup(name='nmrformd',
      version='0.1',
      description='Calculate NMR relaxation time from molecular dynamics trajectory file',
      long_description=open('README.md').read(),
      url='https://github.com/simongravelle/nmrformd',
      author='Simon Gravelle',
      author_email='simon.gravelle@live.fr',
      license='MIT',
      packages=['nmrformd'],
      zip_safe=False,
         install_requires=[
       "mdanalysis",
       "pytest",
       "numpy",
       "random",
       "scipy",
      ]      
      )
