import math

import libhankel
import numpy as np
import pytest


def exponential_ff(q, _params):
    return math.exp(-q)


def test_hankel_G0_builtin_sphere_regression():
    result = libhankel.hankel_G0("sphere", [10.0, 1.0], 1e-9)
    np.testing.assert_allclose(result, 394781.721568501, rtol=1e-6)


def test_hankel_G0_builtin_gdab_regression():
    result = libhankel.hankel_G0("gdab", [10.0, 0.5, 1e-4], 1e-9)
    np.testing.assert_allclose(result, 0.0315827340834857, rtol=1e-6)


def test_hankel_G0_custom_callable_analytical():
    # int_0^inf q * exp(-q) dq = Gamma(2) = 1
    result = libhankel.hankel_G0(exponential_ff, [], 1e-9)
    assert abs(result - 1.0) < 1e-9


def test_hankel_G0_raises_for_invalid_builtin_name():
    with pytest.raises(ValueError):
        libhankel.hankel_G0("not_a_form_factor", [], 1e-6)
