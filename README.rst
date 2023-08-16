====================
elmer_circuitbuilder
====================


.. image:: https://img.shields.io/pypi/v/elmer_circuitbuilder.svg
        :target: https://pypi.python.org/pypi/elmer_circuitbuilder

.. image:: https://readthedocs.org/projects/elmer-circuitbuilder/badge/?version=latest
        :target: https://github.com/ElmerCSC/elmer_circuitbuilder
        :alt: Documentation Status


Python module for creating the circuit simulation definitions for Elmer FEM. The circuit definitions enable easy setup of coils (e.g. massive, stranded, and foil) in 2D and 3D for magnetodynamics applications.


* Free software: GNU Lesser General Public License v3
* Documentation: https://github.com/ElmerCSC/elmer_circuitbuilder.
* Use Examples: https://github.com/ElmerCSC/elmer-elmag/tree/main/CircuitBuilder


Features
--------
* Assembles stiffness and damping matrices of electrical circuit networks
* Available electrical components: Resistors, Capacitors, Inductors, Ideal Current Source, Ideal Voltage Source
* Enables coupling to ElmerFEM to perform circuit-field simulations by setting up coils (massive, stranded, foil) in 2D and 3D

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
