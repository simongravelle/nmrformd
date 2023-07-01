Installation
============

Install the last published version
----------------------------------

Simply type in a terminal:

.. code-block:: bash

	pip3 install nmrformd


Install the last (development) version
--------------------------------------

To get the last version of NMRforMD, clone the repository from `Github`_ on your computer
and use pip3 from the main directory:

.. _`Github`: https://github.com/simongravelle/nmrformd

.. code-block:: bash

	git clone https://github.com/simongravelle/nmrformd.git
	
	cd nmrformd/

	pip install .
	
You can run the tests using pytest:
	
.. code-block:: bash	
	
	cd tests
	pytest .