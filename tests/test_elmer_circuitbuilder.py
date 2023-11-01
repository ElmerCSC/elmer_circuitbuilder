#!/usr/bin/env python

"""Tests for `elmer_circuitbuilder` package."""

import pytest
from elmer_circuitbuilder.elmer_circuitbuilder import *

'''
@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string
'''


def test_helloworld_no_params():
    assert say_hello() == "Hello, World!"

def test_helloworld_with_param():
    assert say_hello("Everyone") == "Hello, Everyone!"


if __name__ == '__main__':
    output_file = "harmonic_open3Dmassive_circuit.definition"

    # initialize circuits: number of circuits - do not remove
    c = number_of_circuits(1)

    # ------------------ Circuit 1 (Current Source - Harmonic)---------------------

    # reference/ground node needed - do not remove.
    c[1].ref_node = 1

    # Components

    I1 = I("source2", 1, 2, 10000)
    V1 = V("source3", 1, 1, 10000)
    V2 = V("source4", 1, 1, 10000)
    V3 = V("source5", 1, 1, 10000)
    Wire1 = ElmerComponent("Winding1", 2, 3, 1)
    Wire1.is3D()
    Wire1.isOpen(2, 1)
    C1 = C("Capacitor1", 2, 1, 100)
    Wire2 = ElmerComponent("Winding2", 3, 1, 2)
    Wire2.is3D()
    Wire2.isOpen(101, 1)

    # store components in array components = [comp1, comp2,...] - do not remove
    c[1].components.append([I1, Wire1, Wire2, C1])
    #c[2].components.append([V1, Wire1])

    # --------------------------------------------------

    # generate elmer circuit.definitions - do not remove / do not edit
    write_elmer_circuits(c, output_file)