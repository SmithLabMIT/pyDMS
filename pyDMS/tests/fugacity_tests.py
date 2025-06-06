import pytest
import numpy as np

import pyDMS.dms as dms


# Pytest unit tests
def test_fugacity_supplied(capfd):
    gas = dms.Gas()
    gas.f = np.array([1,2,3],[4,5,6])
    gas.temp = np.array([300,400])
    dms.calculate_fugacity(gas)
    out, _ = capfd.readouterr()
    assert "Fugacity data supplied in Gas.f" in out

def test_user_supplied_virial(capfd):
    gas = dms.Gas()
    gas.p = np.array([1,2,3])
    gas.temp = np.array(300)
    gas.virial_coeff = {'B0': 1, 'B1': 1}
    dms.calculate_fugacity(gas)
    out, _ = capfd.readouterr()
    assert "Calculating fugacity using user-supplied Virial coefficients" in out

def test_user_supplied_pr(capfd):
    gas = dms.Gas()
    gas.p = np.array([1,2,3])
    gas.temp = np.array(300)
    gas.pr_coeff = {'Tc': 100, 'Pc': 1, 'omega': 0.1}
    dms.calculate_fugacity(gas)
    out, _ = capfd.readouterr()
    assert "Calculating fugacity using user-supplied Peng-Robinson coefficients" in out

def test_builtin_virial(capfd):
    gas = dms.Gas()
    gas.p = np.array([1,2,3])
    gas.temp = np.array(300)
    gas.formula = 'CO2'
    dms.calculate_fugacity(gas)
    out, _ = capfd.readouterr()
    assert "Calculating fugacity using built-in Virial EoS" in out

def test_builtin_pr(capfd):
    gas = dms.Gas()
    gas.p = np.array([1,2,3])
    gas.temp = np.array(300)
    gas.formula = 'H2S'
    dms.calculate_fugacity(gas)
    out, _ = capfd.readouterr()
    assert "Calculating fugacity using built-in Peng-Robinson EoS" in out

def test_formula_not_found(capfd):
    gas = dms.Gas()
    gas.p = np.array([1,2,3])
    gas.temp = np.array(300)
    gas.formula = 'XYZ'
    with pytest.warns(UserWarning, match="Fugacity data unspecified"):
        dms.calculate_fugacity(gas)