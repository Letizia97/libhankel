import libhankel
import numpy as np
import pytest

INPUT_X_ARR = np.array(
    [
        30.0,
        49.0,
        69.0,
        88.0,
        108.0,
        127.0,
        147.0,
        167.0,
        186.0,
        206.0,
        225.0,
        225.0,
        245.0,
        265.0,
        284.0,
        304.0,
        323.0,
        343.0,
        362.0,
        382.0,
        402.0,
        421.0,
        441.0,
        460.0,
        480.0,
    ]
)

QWE_p_dict = {"n_eval": 250, "eps_rel": 1e-9}

QWE_p_dict_missing_epsrel = {
    "n_eval": 250,
}

QWE_p_dict_missing_n_eval = {"eps_rel": 1e-9}

DE_Ooura_p_dict = {"n_eval": 150, "eps_rel": 1e-3}

DE_Ogata_p_dict = {"n_eval": 250, "f_max": 1e-9}


EXPECTED_QWE_CHAVE = np.array(
    [
        15.347799569863255,
        14.754529001964462,
        13.846264717613977,
        12.739784556790136,
        11.353852485895926,
        9.86578297557279,
        8.165195138959756,
        6.377754678115403,
        4.647758740014918,
        2.8440418142266046,
        1.1928792455066153,
        1.1928792455066153,
        -0.4337308219482863,
        -1.9039958830278936,
        -3.122688416817851,
        -4.19016323871282,
        -4.980787075862623,
        -5.5657139669079445,
        -5.88293908666179,
        -5.970271973755689,
        -5.8172944576660015,
        -5.467798696558961,
        -4.9106780365473695,
        -4.230423782590558,
        -3.390155516018234,
    ]
)

EXPECTED_QWE_KEY = np.array(
    [
        15.347799569869006,
        14.754529001968024,
        13.846264717615314,
        12.739784556785377,
        11.353852485861177,
        9.86578297557316,
        8.16519513895994,
        6.377754678115468,
        4.647758748702756,
        2.8440418158982874,
        1.192879245506616,
        1.192879245506616,
        -0.433730821894878,
        -1.9039958830278967,
        -3.12268841681787,
        -4.190163238712818,
        -4.980787075862545,
        -5.565713966907964,
        -5.8829390866623905,
        -5.970271973754746,
        -5.817294457665575,
        -5.467798696557315,
        -4.9106780365336915,
        -4.2304237824456115,
        -3.390155513781909,
    ]
)

EXPECTED_DE_OOURA = np.array(
    [
        15.347814688723965,
        14.754529226195814,
        13.846266014908949,
        12.739785224181167,
        11.353852723194382,
        9.865783015355285,
        8.165195206905562,
        6.376933532411812,
        4.647758765440497,
        2.8440418212818255,
        1.1928792447535888,
        1.1928792447535888,
        -0.43373082212322694,
        -1.9039958859283146,
        -3.1226884198974774,
        -4.190163248317258,
        -4.980787089609459,
        -5.565713976846757,
        -5.882939103327825,
        -5.970271989399333,
        -5.817294475597886,
        -5.467798708798445,
        -4.910678046707662,
        -4.230423786986428,
        -3.3901555158478267,
    ]
)

EXPECTED_DE_OGATA = np.array(
    [
        5.215540452856403e-09,
        1.952950546187973e-09,
        9.844087961421144e-10,
        6.050585360585201e-10,
        4.0164339516349116e-10,
        2.9042339458983e-10,
        2.1675363938056697e-10,
        1.6793423865941656e-10,
        1.3537068422304434e-10,
        1.1035641200492283e-10,
        9.250228228830715e-11,
        9.250228228830715e-11,
        7.801398405385191e-11,
        6.668097796146051e-11,
        5.8056129830674474e-11,
        5.066748929214705e-11,
        4.4881219548991537e-11,
        3.9799265421055466e-11,
        3.5730634827358127e-11,
        3.208677153972563e-11,
        2.8973170771239905e-11,
        2.6416785034483806e-11,
        2.4074826019893963e-11,
        2.2126939769657564e-11,
        2.0321292746578363e-11,
    ]
)

EXPECTED_DHT_6 = np.array(
    [
        21.98842036499114,
        10.816528401842662,
        19.857681627890095,
        14.76541564845125,
        13.710227009818745,
        13.705536807727487,
        5.817441659395433,
        6.485201112071521,
        4.269668058176977,
        3.1817339811677123,
        0.48286782320950133,
        0.48286782320950133,
        -0.21072123677388235,
        -1.1024963500005902,
        -3.9997998700602806,
        -4.940576282963947,
        -3.0523050153848525,
        -5.1011871743792785,
        -8.67046291402422,
        -6.072364369611737,
        -3.5076777944765234,
        -4.254808741350566,
        -6.660015186684702,
        -6.011895842184716,
        -3.138755904417303,
    ]
)

EXPECTED_DHT_7 = np.array(
    [
        23.184913768050198,
        24.751740035803923,
        14.72245755962198,
        19.092660914196305,
        16.394330074972324,
        4.391754030395695,
        14.536248779621479,
        2.776592465195398,
        6.758447662443891,
        2.929113432076164,
        0.192537084988817,
        0.192537084988817,
        -1.3152832567132864,
        -1.0539578898267064,
        -1.234081934544303,
        -4.064027344098499,
        -9.039585519101509,
        -6.19351154835753,
        -2.8202913537305534,
        -2.676770725561007,
        -5.549142391323917,
        -9.606403099542248,
        -8.185871820544921,
        -4.094751526321018,
        -1.6772283185838566,
    ]
)

EXPECTED_DHT_8 = np.array(
    [
        5.159657102916072,
        12.60711748738437,
        10.91039875121369,
        8.874382496969877,
        7.601780506491344,
        13.171462377076196,
        9.019075952990963,
        4.075820926995248,
        6.867869562086511,
        2.0873845984315027,
        0.7088050776732974,
        0.7088050776732974,
        0.029811420342042826,
        -1.6169960894530955,
        -4.561584624825414,
        -3.2687524301332807,
        -3.157795766313723,
        -7.09516798441115,
        -8.223983807538522,
        -4.414317359320974,
        -3.5713373368184693,
        -5.735794260186137,
        -7.370960138036142,
        -4.933510499866119,
        -2.317867654528054,
    ]
)

EXPECTED_DHT_9 = np.array(
    [
        89.14690771745127,
        -4.216519147678679,
        12.26454195127962,
        6.542138836701038,
        9.658528247746444,
        9.423980584346149,
        6.6630303854593285,
        4.927966998947857,
        6.575940002710602,
        2.006248428949638,
        0.7944017997368219,
        0.7944017997368219,
        0.22282374234455019,
        -1.5668984265692611,
        -4.112061242833336,
        -3.6054074260143913,
        -5.440183380970249,
        -5.254616601357,
        -5.105440630050645,
        -7.092114176228941,
        -5.539490825462417,
        -4.704747656146645,
        -5.685423678879062,
        -4.686915558040081,
        -2.764170653976346,
    ]
)

EXPECTED_DHT_10 = np.array(
    [
        14.771625088102969,
        14.489437040851723,
        14.895597160982659,
        13.748865558951504,
        10.833464382681628,
        10.619741249205504,
        8.620407501629384,
        6.689190513979386,
        5.040338909609998,
        2.9059049117713918,
        1.0512202921036806,
        1.0512202921036806,
        -0.30188943980349575,
        -2.1080829684706863,
        -2.8459733470501596,
        -4.537679828331834,
        -4.564342955637244,
        -6.0499161953790255,
        -5.437429557424883,
        -6.26551984122112,
        -5.809710118610776,
        -5.098023929516242,
        -5.359677992239709,
        -4.113025824589891,
        -3.1340433478541456,
    ]
)


@pytest.mark.parametrize(
    "x_arr, strategy_name, strategy_p_dict, expected",
    [
        (INPUT_X_ARR, "QWE_Chave", QWE_p_dict, EXPECTED_QWE_CHAVE),
        (INPUT_X_ARR, "QWE_Key", QWE_p_dict, EXPECTED_QWE_KEY),
        (INPUT_X_ARR, "DE_Ogata", DE_Ogata_p_dict, EXPECTED_DE_OGATA),
        (INPUT_X_ARR, "DE_Ooura", DE_Ooura_p_dict, EXPECTED_DE_OOURA),
        (INPUT_X_ARR, "DHT_6", {}, EXPECTED_DHT_6),
        (INPUT_X_ARR, "DHT_7", {}, EXPECTED_DHT_7),
        (INPUT_X_ARR, "DHT_8", {}, EXPECTED_DHT_8),
        (INPUT_X_ARR, "DHT_9", {}, EXPECTED_DHT_9),
        (INPUT_X_ARR, "DHT_10", {}, EXPECTED_DHT_10),
    ],
)
def test_hankel_transform_strategies_with_builtin_form_factor(
    x_arr, strategy_name, strategy_p_dict, expected
):
    """
    Test that hankel_transform strategies returns expected result,
    when called with a built-in form factor.
    """
    nu = 0
    params_broad_peak = np.array([10e5, 1000, 0.01, 2, 2])

    result = libhankel.hankel_transform(
        nu,
        "broad_peak",
        x_arr,
        params_broad_peak,
        strategy_name,
        strategy_p_dict,
    )
    np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize(
    "x_arr, strategy_name, strategy_p_dict",
    [
        (INPUT_X_ARR, "QWE_Chave", QWE_p_dict),
        (INPUT_X_ARR, "DHT_10", {}),
    ],
)
def test_hankel_transform_returns_error_when_nu_wrong(
    x_arr,
    strategy_name,
    strategy_p_dict,
):
    """
    Test that the Hankel transform function raises a ValueError
    when nu is wrong, i.e. neither 0 nor 1.
    """
    nu = 2
    params_broad_peak = np.array([10e5, 1000, 0.01, 2, 2])
    expected_error = "nu needs to be 0 or 1 in order to use the selected strategy"

    with pytest.raises(ValueError, match=expected_error):
        libhankel.hankel_transform(
            nu,
            "broad_peak",
            x_arr,
            params_broad_peak,
            strategy_name,
            strategy_p_dict,
        )


@pytest.mark.parametrize(
    "x_arr, strategy_name, strategy_p_dict",
    [
        (INPUT_X_ARR, "QWE_Chave", QWE_p_dict_missing_epsrel),
        (INPUT_X_ARR, "QWE_Key", QWE_p_dict_missing_epsrel),
    ],
)
def test_hankel_transform_QWE_raises_error_when_eps_rel_missing(
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
    params_broad_peak = np.array([10e5, 1000, 0.01, 2, 2])
    expected_error = "Error: eps_rel must be provided and cannot be zero"

    with pytest.raises(ValueError, match=expected_error):
        libhankel.hankel_transform(
            nu,
            "broad_peak",
            x_arr,
            params_broad_peak,
            strategy_name,
            strategy_p_dict,
        )


@pytest.mark.parametrize(
    "x_arr, strategy_name, strategy_p_dict",
    [
        (INPUT_X_ARR, "QWE_Chave", QWE_p_dict_missing_n_eval),
        (INPUT_X_ARR, "QWE_Key", QWE_p_dict_missing_n_eval),
    ],
)
def test_hankel_transform_QWE_raises_error_when_n_eval_missing(
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
    params_broad_peak = np.array([10e5, 1000, 0.01, 2, 2])
    expected_error = "Error: n_eval must be provided and cannot be zero"

    with pytest.raises(ValueError, match=expected_error):
        libhankel.hankel_transform(
            nu,
            "broad_peak",
            x_arr,
            params_broad_peak,
            strategy_name,
            strategy_p_dict,
        )


def test_hankel_transform_DE_Ogata_raises_error_when_fmax_missing():
    """
    Test that the Hankel transform function raises a ValueError
    when f_max is missing from the strategy_params dict, but
    it is required by the strategy.
    """
    x_arr = np.array([1, 2, 3])
    nu = 0
    params_broad_peak = np.array([10e5, 1000, 0.01, 2, 2])
    expected_error = "Error: f_max must be provided and cannot be zero"

    with pytest.raises(ValueError, match=expected_error):
        libhankel.hankel_transform(
            nu,
            "broad_peak",
            x_arr,
            params_broad_peak,
            "DE_Ogata",
            {"n_eval": 250},
        )
