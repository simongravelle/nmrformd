# NMRforMD - developer section

<a href="webp">
  <img src="../docs/source/figures/logo/logo-b.png" align="right" width="30%"/>
</a>

## How to

If you intend to make modification to the code, please raise an issue or send me an email
first. Then, fork the repository, apply your changes, then make a pull request
that will be reviewed.

## Build the documentation

Clone the repository as follow:

``` bash
   git clone https://github.com/simongravelle/nmrformd.git --recurse-submodule
```

Build the documentation locally from the [docs](docs/) folder by typing:

``` bash
    pip3 install -r requirements.txt
    make clean
    make html
```

## Publish new pip version (1)

Publish a new pip version by following those
[instructions](https://gist.github.com/arsho/fc651bfadd8a0f42be72156fd21bd8a9).

1 - if necessary, update *docs/source/conf.py*, *CITATION.cff*, and *setup.py*

2 - Create source distribution using

``` bash
    python3 setup.py sdist
```

3 - Create a new release on Github using the generated tar.gz file located in dist/

4 - Update the link in setup.py

5 - Create wheel using:

``` bash
    python3 setup.py bdist_wheel
```

6 - Upload to pypi using (with the appropriate number):

``` bash
    twine upload dist/nmrformd-0.1.0*
```

![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)
![readthedoc](https://readthedocs.org/projects/nmrformd/badge/?version=latest)

## Publish new pip version (2)

This instruction are from [this page](https://gist.github.com/arsho/fc651bfadd8a0f42be72156fd21bd8a9).

Install last version using 

```bash
    pip install -e .
```

1 - Update version in conf.py

2 - Update CITATION.cff

3 - Update setup.py

4 - Update CHANGELOG

5 - Create source distribution using python3 setup.py sdist

6 - Create a new release on Github using the generated tar.gz file

7 - Update the link in setup.py

8 - Create wheel using python3 setup.py bdist_wheel

9 - Upload to pypi using twine upload dist/nmrformd-0.1.0* \


