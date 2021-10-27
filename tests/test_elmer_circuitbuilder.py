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
