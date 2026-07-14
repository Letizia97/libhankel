import math

import libhankel
import numpy as np
import pytest

INPUT_X_ARR = np.array(
    [
        15.0,
        18.54166667,
        22.08333333,
        25.625,
        29.16666667,
        32.70833333,
        36.25,
        39.79166667,
        43.33333333,
        46.875,
        50.41666667,
        53.95833333,
        57.5,
        61.04166667,
        64.58333333,
        68.125,
        71.66666667,
        75.20833333,
        78.75,
        82.29166667,
        85.83333333,
        89.375,
        92.91666667,
        96.45833333,
        100.0,
    ]
)

QWE_p_dict = {"n_eval": 250, "eps_rel": 1e-9}

QWE_p_dict_missing_epsrel = {
    "n_eval": 250,
}

QWE_p_dict_missing_n_eval = {"eps_rel": 1e-9}

DE_Ooura_p_dict = {"n_eval": 150, "eps_rel": 1e-3}

DE_Ogata_p_dict = {"n_eval": 250, "f_max": 1e-9}


def sf_poch(a, x):
    return math.gamma(a + x) / math.gamma(a)


def dab(q, params):
    """
    Example Debye-Anderson-Brumberger (DAB) model.
    """
    XI = params[0]
    H = params[1]
    ETA = params[2]

    factor1 = (2 * XI) ** 3
    factor2 = sf_poch(H, 1.5)

    numer = (factor1 * factor2 * ETA) ** 2 * (np.pi**3)
    denom = (1 + (q * XI) ** 2) ** (1.5 + H)

    return numer / denom


@pytest.mark.parametrize(
    "form_factor, x_arr, strategy_name, strategy_p_dict",
    [
        (dab, INPUT_X_ARR, "QWE_Chave", QWE_p_dict),
        (dab, INPUT_X_ARR, "DHT_10", {}),
        ("gdab", INPUT_X_ARR, "QWE_Chave", QWE_p_dict),
        ("gdab", INPUT_X_ARR, "DHT_10", {}),
    ],
)
def test_hankel_transform_returns_error_when_nu_wrong(
    form_factor,
    x_arr,
    strategy_name,
    strategy_p_dict,
):
    """
    Test that the Hankel transform function raises a ValueError
    when nu is wrong, i.e. neither 0 nor 1.
    """
    nu = 2
    params_gdab = np.array([10.0, 0.5, 1e-4])
    expected_error = "nu needs to be 0 or 1 in order to use the selected strategy"

    with pytest.raises(ValueError, match=expected_error):
        libhankel.hankel_transform(
            nu,
            form_factor,
            x_arr,
            params_gdab,
            strategy_name,
            strategy_p_dict,
        )


@pytest.mark.parametrize(
    "form_factor, x_arr, strategy_name, strategy_p_dict",
    [
        (dab, INPUT_X_ARR, "QWE_Chave", QWE_p_dict_missing_epsrel),
        (dab, INPUT_X_ARR, "QWE_Key", QWE_p_dict_missing_epsrel),
        ("gdab", INPUT_X_ARR, "QWE_Chave", QWE_p_dict_missing_epsrel),
        ("gdab", INPUT_X_ARR, "QWE_Key", QWE_p_dict_missing_epsrel),
    ],
)
def test_hankel_transform_QWE_raises_error_when_eps_rel_missing(
    form_factor,
    x_arr,
    strategy_name,
    strategy_p_dict,
):
    """
    Test that the Hankel transform function raises a ValueError
    when eps_rel is missing from the strategy_params dict, but
    it is required by the strategy.
    """
    nu = 0
    params_gdab = np.array([10.0, 0.5, 1e-4])
    expected_error = "Error: eps_rel must be provided and cannot be zero"

    with pytest.raises(ValueError, match=expected_error):
        libhankel.hankel_transform(
            nu,
            form_factor,
            x_arr,
            params_gdab,
            strategy_name,
            strategy_p_dict,
        )


@pytest.mark.parametrize(
    "form_factor, x_arr, strategy_name, strategy_p_dict",
    [
        (dab, INPUT_X_ARR, "QWE_Chave", QWE_p_dict_missing_n_eval),
        (dab, INPUT_X_ARR, "QWE_Key", QWE_p_dict_missing_n_eval),
        ("gdab", INPUT_X_ARR, "QWE_Chave", QWE_p_dict_missing_n_eval),
        ("gdab", INPUT_X_ARR, "QWE_Key", QWE_p_dict_missing_n_eval),
    ],
)
def test_hankel_transform_QWE_raises_error_when_n_eval_missing(
    form_factor,
    x_arr,
    strategy_name,
    strategy_p_dict,
):
    """
    Test that the Hankel transform function raises a ValueError
    when n_eval is missing from the strategy_params dict, but
    it is required by the strategy.
    """
    nu = 0
    params_gdab = np.array([10.0, 0.5, 1e-4])
    expected_error = "Error: n_eval must be provided and cannot be zero"

    with pytest.raises(ValueError, match=expected_error):
        libhankel.hankel_transform(
            nu,
            form_factor,
            x_arr,
            params_gdab,
            strategy_name,
            strategy_p_dict,
        )


@pytest.mark.parametrize(
    "form_factor, x_arr",
    [
        (dab, INPUT_X_ARR),
        ("gdab", INPUT_X_ARR),
    ],
)
def test_hankel_transform_DE_Ogata_raises_error_when_fmax_missing(
    form_factor,
    x_arr,
):
    """
    Test that the Hankel transform function raises a ValueError
    when f_max is missing from the strategy_params dict, but
    it is required by the strategy.
    """
    nu = 0
    params_gdab = np.array([10.0, 0.5, 1e-4])
    expected_error = "Error: f_max must be provided and cannot be zero"

    with pytest.raises(ValueError, match=expected_error):
        libhankel.hankel_transform(
            nu,
            form_factor,
            x_arr,
            params_gdab,
            "DE_Ogata",
            {"n_eval": 250},
        )


@pytest.mark.parametrize(
    "form_factor, x_arr",
    [
        (dab, INPUT_X_ARR),
        ("gdab", INPUT_X_ARR),
    ],
)
def test_hankel_transform_raises_correct_error_when_not_converged(
    form_factor,
    x_arr,
):
    nu = 0
    params_gdab = np.array([10.0, 0.5, 1e-4])
    expected_error = "Failed to converge"

    with pytest.raises(RuntimeError, match=expected_error):
        libhankel.hankel_transform(
            nu,
            form_factor,
            x_arr,
            params_gdab,
            "QWE_Chave",
            {"n_eval": 4, "eps_rel": 1e-9},
        )
