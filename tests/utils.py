""" Common utilities for the MIRaGe test suite.

Authors
-------
    - Lauren Chambers

"""

import os
import yaml

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

def parametrized_data():
    """Load parametrized data from file.

    Returns
    -------
    test_data : dict
        Dictionary containing parametrized test data
    """
    parametrized_data_file = os.path.join(__location__, 'test_data', 'parametrized_test_data.yaml')
    with open(parametrized_data_file) as f:
        test_data = yaml.safe_load(f.read())

    return test_data