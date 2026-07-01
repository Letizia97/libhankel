.. _strategy-parameters:

Strategy Parameters
================================

Independently of the interface being used, the ``hankel_transform`` function is designed to receive a 
parameter named ``strategy_name``, which allows the user the select one of the strategies available.
This could be any of the strings in the "Strategy Name" column in the table below.

The user is also required to supply a ``strategy_params`` input, which will contain the 
parameters required by the particular strategy being used. 
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

- ``n_eval``: integer indicating number of function evaluations

- ``eps_rel``: relative error allowed e.g. 1e-9 

- ``f_max``: 	float indicating starting guess for max in form factor

Please note, for DHT strategies (that is, **DHT_6** to **DHT_11**) feel free to supply an empty struct 
for ``strategy_params`` . Any field within it will be ignored in these strategies as they do not 
require such additional parameters.  
