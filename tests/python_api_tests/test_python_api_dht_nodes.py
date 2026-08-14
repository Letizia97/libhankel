import re

import libhankel
import numpy as np
import pytest

INPUT_X_ARR = np.array([30.0, 147.0, 402.0, 480.0])

FILTER_LENGTHS = [
    ("DHT_Key_51", 51),
    ("DHT_Key_101", 101),
    ("DHT_Key_201", 201),
    ("DHT_Anderson_801", 801),
]

ADAPTIVE_STRATEGIES = [
    "DHT_Guptasarma",
    "DHT_Guptasarma_Fast",
    "Fixed_DE_Ogata",
    "Adaptive_DE_Ooura",
    "QWE_Chave",
    "QWE_Key",
]


def decaying_exponential(x, params):
    return np.exp(-params[0] * x)


@pytest.mark.parametrize("strategy_name, length", FILTER_LENGTHS)
@pytest.mark.parametrize("nu", [0, 1])
def test_dht_nodes_returns_one_abscissa_and_one_weight_per_filter_point(strategy_name, length, nu):
    abscissae, weights = libhankel.dht_nodes(strategy_name, nu)

    assert len(abscissae) == length
    assert len(weights) == length


@pytest.mark.parametrize("strategy_name, _length", FILTER_LENGTHS)
def test_dht_nodes_abscissae_are_positive_and_increasing(strategy_name, _length):
    abscissae, _weights = libhankel.dht_nodes(strategy_name, 0)
    abscissae = np.asarray(abscissae)

    assert np.all(abscissae > 0)
    assert np.all(np.diff(abscissae) > 0)


@pytest.mark.parametrize("strategy_name, _length", FILTER_LENGTHS)
def test_dht_nodes_returns_a_different_weight_column_for_each_order(strategy_name, _length):
    abscissae_j0, weights_j0 = libhankel.dht_nodes(strategy_name, 0)
    abscissae_j1, weights_j1 = libhankel.dht_nodes(strategy_name, 1)

    assert abscissae_j0 == abscissae_j1
    assert weights_j0 != weights_j1


@pytest.mark.parametrize("strategy_name, _length", FILTER_LENGTHS)
@pytest.mark.parametrize("nu", [0, 1])
def test_dht_nodes_reproduce_hankel_transform(strategy_name, _length, nu):
    """
    The point of the accessor: a caller who evaluates f at the abscissae itself
    and applies the weights must get what hankel_transform would have returned.
    """
    params = np.array([0.01])
    expected = np.asarray(
        libhankel.hankel_transform(nu, decaying_exponential, INPUT_X_ARR, params, strategy_name, {})
    )

    abscissae, weights = libhankel.dht_nodes(strategy_name, nu)
    abscissae = np.asarray(abscissae)
    weights = np.asarray(weights)

    q = abscissae[:, None] / INPUT_X_ARR[None, :]
    f = decaying_exponential(q, params)
    result = np.sum(f * q * weights[:, None] / INPUT_X_ARR[None, :], axis=0)

    np.testing.assert_allclose(result, expected, rtol=1e-12)


@pytest.mark.parametrize("strategy_name", ADAPTIVE_STRATEGIES)
def test_dht_nodes_raises_correct_error_for_strategies_without_a_node_table(strategy_name):
    expected_error = (
        "Error: nodes are only tabulated for the fixed-abscissa filters, one of: "
        "'DHT_Key_51', 'DHT_Key_101', 'DHT_Key_201', 'DHT_Anderson_801'."
    )

    with pytest.raises(ValueError, match=re.escape(expected_error)):
        libhankel.dht_nodes(strategy_name, 0)


def test_dht_nodes_raises_correct_error_for_unknown_strategy_name():
    with pytest.raises(ValueError):
        libhankel.dht_nodes("DHT_1", 0)


@pytest.mark.parametrize("nu", [-1, 2])
def test_dht_nodes_raises_correct_error_when_nu_is_out_of_range(nu):
    with pytest.raises(ValueError, match=re.escape("nu must be 0 or 1")):
        libhankel.dht_nodes("DHT_Key_201", nu)
