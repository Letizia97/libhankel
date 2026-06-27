.. _strategy-parameters:

Strategy Parameters
-------------------


The ``hankel_transform`` function is designed to receive ...

The parameter ``strategy_name`` allows the user the select one of the strategies available.
This could be any of the strings in the "Strategy Name" column in the table below.

The user is also required to supply ``strategy_params``, which will have different
form depending on whether the Python or C API is used, but in both cases will need to contain 
a number of parameters which vary based on ``strategy_name``. These are reported in the table 
below in the last column on the left. It's important that ``strategy_params`` contains the right 
parameters required by the selected strategy. If it doesn't, an error will be thrown. 


+--------------------------+---------------------+--------------------------------------------+
| index (as in SASfit)     | Strategy Name       | Field names needed in ``strategy_params``  |
+==========================+=====================+============================================+
| 0                        | DE_Ooura            | ``n_eval``, ``eps_rel``                    |
+--------------------------+---------------------+--------------------------------------------+
| 1                        | DE_Ogata            | ``n_eval``, ``f_max``                      |
+--------------------------+---------------------+--------------------------------------------+
| 6 to 11                  | DHT_6 to DHT_11     | N/A                                        |
+--------------------------+---------------------+--------------------------------------------+
| 12                       | QWE_Key             | ``n_eval``, ``eps_rel``                    |
+--------------------------+---------------------+--------------------------------------------+
| 13                       | QWE_Chave           | ``n_eval``, ``eps_rel``                    |
+--------------------------+---------------------+--------------------------------------------+


Here is an explanation of each parameter:
- ``n_eval``: integer indicating number of function evaluations
- ``eps_rel``: relative error allowed e.g. 1e-9 
- ``f_max``: 	float indicating starting guess for max in form factor

Please note, for DHT strategies (that is, DHT_6 to DHT_11) feel free to supply an empty struct 
for ``strategy_params`` . Any field within it will be ignored in these strategies as they do not 
require such additional parameters.  


.. _status-codes:

Status codes
------------

When successful, the ``hankel_transform`` function returns 0 . When unsuccessful, it returns an error code. 
To find out more about the error, the table below lists all possible error codes one might receive plus 
corresponding explanations.

+--------------------+-------------------------------------------------------------------------+
| Status codes       | Explanation                                                             |      
+====================+=========================================================================+
| 0                  | Success                                                                 |
+--------------------+-------------------------------------------------------------------------+
| -1                 | Wrong nu (order of bessel function)                                     |                                              
+--------------------+-------------------------------------------------------------------------+
| -2                 | Wrong number for DHT strategy, must be int between 6 and 11             |
+--------------------+-------------------------------------------------------------------------+
| -3                 | Failed to allocate variables                                            |
+--------------------+-------------------------------------------------------------------------+
| -4                 | Did not converge                                                        |
+--------------------+-------------------------------------------------------------------------+
| -5                 | Internal error: division by zero                                        |
+--------------------+-------------------------------------------------------------------------+
| -6                 | Internal error: wrong nzeros in function bessel_j_zero (must be >= 1)   |
+--------------------+-------------------------------------------------------------------------+
| -7                 | Internal error: wrong n of iterations in pade sum (must be >= 1)        |
+--------------------+-------------------------------------------------------------------------+
| -8                 | Wrong strategy parameters supplied to hankel transform function         |
+--------------------+-------------------------------------------------------------------------+
| -9                 | Wrong parameters supplied to form factor                                |
+--------------------+-------------------------------------------------------------------------+