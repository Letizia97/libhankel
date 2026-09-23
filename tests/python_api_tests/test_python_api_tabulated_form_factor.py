import libhankel
import numpy as np
import pytest


# Simple tabulated form factor data (sphere-like)
Q_DATA = np.array([0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0])
F_DATA = np.array([1.0, 0.98, 0.92, 0.72, 0.48, 0.28, 0.15, 0.05, 0.01, 0.002])

X_ARR = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0])


def test_tabulated_form_factor_basic():
    """Test basic tabulated form factor with linear interpolation and power_law tail."""
    nu = 0
    result = libhankel.hankel_transform(
        nu,
        {
            "q": Q_DATA,
            "f": F_DATA,
            "interp_type": "linear",
            "tail": "power_law",
        },
        X_ARR,
        [],
        "DHT_Key_201",
        {},
    )

    assert result is not None
    assert len(result) == len(X_ARR)
    assert all(np.isfinite(result))


def test_tabulated_form_factor_cubic_interp():
    """Test tabulated form factor with cubic interpolation."""
    nu = 0
    result = libhankel.hankel_transform(
        nu,
        {
            "q": Q_DATA,
            "f": F_DATA,
            "interp_type": "cubic",
            "tail": "zero",
        },
        X_ARR,
        [],
        "DHT_Key_201",
        {},
    )

    assert result is not None
    assert len(result) == len(X_ARR)
    assert all(np.isfinite(result))


def test_tabulated_form_factor_loglinear_interp():
    """Test tabulated form factor with log-linear interpolation."""
    nu = 0
    result = libhankel.hankel_transform(
        nu,
        {
            "q": Q_DATA,
            "f": F_DATA,
            "interp_type": "loglinear",
            "tail": "power_law",
        },
        X_ARR,
        [],
        "DHT_Key_201",
        {},
    )

    assert result is not None
    assert len(result) == len(X_ARR)
    assert all(np.isfinite(result))


def test_tabulated_missing_q_key():
    """Test that missing 'q' key raises ValueError."""
    nu = 0
    with pytest.raises(ValueError, match="must contain"):
        libhankel.hankel_transform(
            nu,
            {
                "f": F_DATA,
                "interp_type": "linear",
                "tail": "power_law",
            },
            X_ARR,
            [],
            "DHT_Key_201",
            {},
        )


def test_tabulated_missing_f_key():
    """Test that missing 'f' key raises ValueError."""
    nu = 0
    with pytest.raises(ValueError, match="must contain"):
        libhankel.hankel_transform(
            nu,
            {
                "q": Q_DATA,
                "interp_type": "linear",
                "tail": "power_law",
            },
            X_ARR,
            [],
            "DHT_Key_201",
            {},
        )


def test_tabulated_missing_interp_type_key():
    """Test that missing 'interp_type' key raises ValueError."""
    nu = 0
    with pytest.raises(ValueError, match="must contain"):
        libhankel.hankel_transform(
            nu,
            {
                "q": Q_DATA,
                "f": F_DATA,
                "tail": "power_law",
            },
            X_ARR,
            [],
            "DHT_Key_201",
            {},
        )


def test_tabulated_missing_tail_key():
    """Test that missing 'tail' key raises ValueError."""
    nu = 0
    with pytest.raises(ValueError, match="must contain"):
        libhankel.hankel_transform(
            nu,
            {
                "q": Q_DATA,
                "f": F_DATA,
                "interp_type": "linear",
            },
            X_ARR,
            [],
            "DHT_Key_201",
            {},
        )


def test_tabulated_invalid_interp_type():
    """Test that invalid interp_type raises ValueError."""
    nu = 0
    with pytest.raises(ValueError, match="interp_type must be"):
        libhankel.hankel_transform(
            nu,
            {
                "q": Q_DATA,
                "f": F_DATA,
                "interp_type": "spline",
                "tail": "power_law",
            },
            X_ARR,
            [],
            "DHT_Key_201",
            {},
        )


def test_tabulated_invalid_tail():
    """Test that invalid tail raises ValueError."""
    nu = 0
    with pytest.raises(ValueError, match="tail must be"):
        libhankel.hankel_transform(
            nu,
            {
                "q": Q_DATA,
                "f": F_DATA,
                "interp_type": "linear",
                "tail": "exponential",
            },
            X_ARR,
            [],
            "DHT_Key_201",
            {},
        )


def test_tabulated_with_explicit_exponent():
    """Test tabulated form factor with explicit exponent for power_law tail."""
    nu = 0
    result = libhankel.hankel_transform(
        nu,
        {
            "q": Q_DATA,
            "f": F_DATA,
            "interp_type": "cubic",
            "tail": "power_law",
            "exponent": 4.0,
        },
        X_ARR,
        [],
        "DHT_Key_201",
        {},
    )

    assert result is not None
    assert len(result) == len(X_ARR)
    assert all(np.isfinite(result))
