.. _strategy-parameters:

Strategy Parameters
================================

The ``hankel_transform`` function is designed to receive a both a ``strategy_name``, 
which allows the user the select one of the strategies available, and ``strategy_params``, 
which will contain the parameters required by the particular strategy being used. 
The latter will be a dict when using the Python API and a struct when using the C API.

The table below reports the parameters required by each strategy.


.. table:: Parameters needed by each strategy.

    +----------------------+--------------------------------------------+----------------------+
    | Strategy Name        | Field names needed in ``strategy_params``  | Index (as in SASfit) |
    +======================+============================================+======================+
    | Adaptive_DE_Ooura    | ``n_eval``, ``eps_rel``                    | 0                    |
    +----------------------+--------------------------------------------+----------------------+
    | Fixed_DE_Ogata       | ``n_eval``, ``f_max``                      | 1                    |
    +----------------------+--------------------------------------------+----------------------+
    | all ``DHT_*``        | N/A                                        | 6 to 11              |
    | filters              |                                            |                      |
    +----------------------+--------------------------------------------+----------------------+
    | QWE_Key              | ``n_eval``, ``eps_rel``                    | 12                   |
    +----------------------+--------------------------------------------+----------------------+
    | QWE_Chave            | ``n_eval``, ``eps_rel``                    | 13                   |
    +----------------------+--------------------------------------------+----------------------+

See :ref:`strategy-selection` for the full list of ``DHT_*`` filter names.


Here is an explanation of each parameter:

.. doxygenstruct:: strategy_params
   :project: libhankel
   :members:


Please note, for the digital-filter strategies (any name starting with **DHT_**) feel free to supply an
empty struct in C (or an empty dict in Python) for ``strategy_params`` . Any field within it will be
ignored in these strategies as they do not require such additional parameters.
