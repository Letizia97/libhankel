import libhankel
import numpy as np
import pytest


# Simple tabulated form factor data (sphere-like)
Q_DATA = np.array([0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0])
F_DATA = np.array([1.0, 0.98, 0.92, 0.72, 0.48, 0.28, 0.15, 0.05, 0.01, 0.002])

X_ARR = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0])

# Dense tabulated sphere FF for numerical accuracy tests (R=10, eta=1).
# Matches form_factor_sphere in C: F(q) = (eta * 4pi * (sin(qR) - qR*cos(qR)) / q^3)^2
# with a Taylor expansion for small qR.
_R, _ETA = 10.0, 1.0
_Q_DENSE = np.logspace(-3, 2, 1000)
_qR = _Q_DENSE * _R
_large = _ETA * 4 * np.pi * (np.sin(_qR) - _qR * np.cos(_qR)) / _Q_DENSE**3
_small = _ETA * (4 / 3) * np.pi * _R**3 * (
    1 - _qR**2 / 10 + _qR**4 / 280 - _qR**6 / 15120
)
_F_DENSE = np.where(_qR < 1e-4, _small, _large) ** 2

# x values in a range where the sphere transform is smooth and well-conditioned
_X_ACCURACY = np.array([1.0, 2.0, 5.0, 10.0, 20.0])


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


def _builtin_sphere_transform(x_arr):
    """Reference: hankel_transform using the built-in sphere form factor."""
    return np.array(
        libhankel.hankel_transform(0, "sphere", x_arr, [_R, _ETA], "DHT_Key_201", {})
    )


def _tabulated_transform(interp_type, tail, x_arr, exponent=None):
    ff = {"q": _Q_DENSE, "f": _F_DENSE, "interp_type": interp_type, "tail": tail}
    if exponent is not None:
        ff["exponent"] = exponent
    return np.array(libhankel.hankel_transform(0, ff, x_arr, [], "DHT_Key_201", {}))


def test_tabulated_loglinear_rejects_non_positive_f():
    """Log-linear interpolation requires all f values to be strictly positive."""
    f_with_zero = F_DATA.copy()
    f_with_zero[3] = 0.0
    with pytest.raises(ValueError, match="Failed to create tabulated form factor"):
        libhankel.hankel_transform(
            0,
            {"q": Q_DATA, "f": f_with_zero, "interp_type": "loglinear", "tail": "power_law"},
            X_ARR,
            [],
            "DHT_Key_201",
            {},
        )


def test_tabulated_loglinear_matches_builtin_sphere():
    """Tabulated sphere with loglinear interp matches built-in to within 1%.
    atol=5 handles near-zero crossings where rtol is unreliable (peak ~4e5)."""
    ref = _builtin_sphere_transform(_X_ACCURACY)
    tab = _tabulated_transform("loglinear", "power_law", _X_ACCURACY)
    np.testing.assert_allclose(tab, ref, rtol=0.01, atol=5.0)


def test_tabulated_cubic_matches_builtin_sphere():
    """Tabulated sphere with cubic interp matches built-in to within 1%."""
    ref = _builtin_sphere_transform(_X_ACCURACY)
    tab = _tabulated_transform("cubic", "power_law", _X_ACCURACY)
    np.testing.assert_allclose(tab, ref, rtol=0.01, atol=5.0)


def test_tabulated_linear_matches_builtin_sphere():
    """Tabulated sphere with linear interp matches built-in to within 2%.
    Linear interpolation of a rapidly oscillating function is less accurate
    near zero crossings, so a slightly looser atol applies there."""
    ref = _builtin_sphere_transform(_X_ACCURACY)
    tab = _tabulated_transform("linear", "power_law", _X_ACCURACY)
    np.testing.assert_allclose(tab, ref, rtol=0.02, atol=5.0)


def test_tabulated_explicit_exponent_matches_builtin_sphere():
    """Tabulated sphere with explicit power-law exponent matches built-in to within 1%."""
    ref = _builtin_sphere_transform(_X_ACCURACY)
    # Sphere FF ~ q^-4 in the Porod regime; exponent=4 in the library convention
    tab = _tabulated_transform("loglinear", "power_law", _X_ACCURACY, exponent=4.0)
    np.testing.assert_allclose(tab, ref, rtol=0.01, atol=5.0)
