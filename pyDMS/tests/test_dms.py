import numpy as np
import pytest

import pyDMS.dms as dms
# data to use for the tests
gas = dms.Gas()
gas.c = np.array([[15,41,56,75,87,111],[8,28,40,57,71,87],[1,5,11,19,28,38],[0.3,2,4,8,11,17]])
gas.cerr = np.array([[0.1,0.2,0.3,0.4,0.5,0.6],[0.1,0.2,0.3,0.4,0.5,0.6],[0.1,0.2,0.3,0.4,0.5,0.6],[0.1,0.2,0.3,0.4,0.5,0.6]])
gas.p = np.array([[0.1,0.5,1,2,3,5],[0.1,0.5,1,2,3,5],[0.1,0.5,1,2,3,5],[0.1,0.5,1,2,3,5],[0.1,0.5,1,2,3,5]])
gas.T = np.array([308,328,398,423])

gas.LFER.out = np.array([0.1, 0.2, 0.3, 0.4]) # *get reall data
    
x_LFER = np.array([1.0, 0.1, 0.5, 0.05, 10, 10, 10, 10]) # * get real data
x_vH = np.array([1.0, 0.1, 10, 10, 10, 10]) # * get real data
# def test_LFER():

# test_vH():

# test_analysis():

# test_Gas():

def test_LFER_loss():

    result = dms.LFER_loss(x_LFER, gas)

    assert isinstance(result, float)
    assert result >= 0.0
    assert result < 1000.0 # * get precise result

def test_vH_loss():
    
    result = dms.vH_loss(x_vH, gas)

    assert isinstance(result, float)
    assert result >= 0.0
    assert result < 1E20 # * get precise result

def test_is_outlier():

    array = np.array([1,2,3,12])

    result = dms.is_outlier(array)

    assert isinstance(result, np.ndarray)
    assert np.array_equal(result, np.array([False, False, False, True]))

#def test_hess():

#def test_calc_LFEs():

# test_calc_params():

# test_calc_error():
if __name__ == "__main__":
    pytest.main()