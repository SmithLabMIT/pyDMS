import pytest
import numpy as np
import pyDMS.dms as dms


# Pytest unit tests
def test_fugacity_supplied(capfd):
    gas = dms.Gas()
    gas.f = np.array([[1, 2, 3], [4, 5, 6]])
    gas.temp = np.array([300, 400])
    dms.calculate_fugacity(gas)
    out, _ = capfd.readouterr()

    assert "Fugacity data supplied in Gas.f" in out


def test_user_supplied_virial(capfd):
    gas = dms.Gas()
    gas.p = np.array([[1, 2, 3], [4, 5, 6]])
    gas.temp = np.array([300, 400])
    gas.virial_coeff = {"B0": 1, "B1": 1}
    dms.calculate_fugacity(gas)
    out, _ = capfd.readouterr()

    assert "Calculating fugacity using user-supplied Virial" in out


def test_user_supplied_pr(capfd):
    gas = dms.Gas()
    gas.p = np.array([[1, 2, 3], [4, 5, 6]])
    gas.temp = np.array([300, 400])
    gas.pr_coeff = {"Tc": 100, "Pc": 1, "omega": 0.1}
    dms.calculate_fugacity(gas)
    out, _ = capfd.readouterr()

    assert "Calculating fugacity using user-supplied Peng-Robinson" in out


def test_builtin_virial(capfd):
    gas = dms.Gas()
    gas.p = np.array([[1, 2, 3], [4, 5, 6]])
    gas.temp = np.array([300, 400])
    gas.formula = "CO2"
    dms.calculate_fugacity(gas)
    out, _ = capfd.readouterr()

    assert "Calculating fugacity using built-in Virial EoS" in out


def test_builtin_pr(capfd):
    gas = dms.Gas()
    gas.p = np.array([[1, 2, 3], [4, 5, 6]])
    gas.temp = np.array([300, 400])
    gas.formula = "H2S"
    dms.calculate_fugacity(gas)
    out, _ = capfd.readouterr()

    assert "Calculating fugacity using built-in Peng-Robinson EoS" in out


def test_builtin_virial_values():
    gas = dms.Gas()
    gas.p = np.array([[1, 2, 3], [4, 5, 6]])
    gas.temp = np.array([300, 400])
    gas.formula = "CO2"
    dms.calculate_fugacity(gas)

    # Check all fugacity values are within 1 atm of pressure
    assert np.allclose(
        gas.f, gas.p, atol=1.0
    ), "Fugacity values deviate by more than 1 atm from pressure"


def test_builtin_pr_values():
    gas = dms.Gas()
    gas.p = np.array([[1, 2, 3], [4, 5, 6]])
    gas.temp = np.array([300, 400])
    gas.formula = "H2S"
    dms.calculate_fugacity(gas)

    # Check all fugacity values are within 1 atm of pressure
    assert np.allclose(
        gas.f, gas.p, atol=1.0
    ), "Fugacity values deviate by more than 1 atm from pressure"


def test_formula_not_found(capfd):
    gas = dms.Gas()
    gas.p = np.array([[1, 2, 3], [4, 5, 5]])
    gas.temp = np.array([300, 400])
    gas.formula = "XYZ"
    with pytest.warns(UserWarning, match="Fugacity data unspecified"):
        dms.calculate_fugacity(gas)


if __name__ == "__main__":
    pytest.main()
