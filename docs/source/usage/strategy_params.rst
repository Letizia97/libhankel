.. _strategy-parameters:

Strategy Parameters
================================

The ``hankel_transform`` function is designed to receive a both a ``strategy_name``, 
which allows the user the select one of the strategies available, and ``strategy_params``, 
which will contain the parameters required by the particular strategy being used. 
The latter will be a dict when using the Python API and a struct when using the C API.

The table below reports the parameters required by each strategy.


.. table:: Parameters needed by each strategy.

    +---------------------+--------------------------------------------+--------------------------+
    | Strategy Name       | Field names needed in ``strategy_params``  | Index (as in SASfit)     |
    +=====================+============================================+==========================+
    | DE_Ooura            | ``n_eval``, ``eps_rel``                    | 0                        |
    +---------------------+--------------------------------------------+--------------------------+
    | DE_Ogata            | ``n_eval``, ``f_max``                      | 1                        |
    +---------------------+--------------------------------------------+--------------------------+
    | DHT_6 to DHT_11     | N/A                                        | 6 to 11                  |
    +---------------------+--------------------------------------------+--------------------------+
    | QWE_Key             | ``n_eval``, ``eps_rel``                    | 12                       |
    +---------------------+--------------------------------------------+--------------------------+
    | QWE_Chave           | ``n_eval``, ``eps_rel``                    | 13                       |
    +---------------------+--------------------------------------------+--------------------------+


Here is an explanation of each parameter:

.. doxygenstruct:: strategy_params
   :project: libhankel
   :members:


Please note, for DHT strategies (that is, **DHT_6** to **DHT_11**) feel free to supply an empty struct in C 
(or an empty dict in Python) for ``strategy_params`` . Any field within it will be ignored in these strategies 
as they do not require such additional parameters.  
