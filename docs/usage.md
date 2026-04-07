# Usage Guide {#usage}

## Main function

The main function a user will be using is `hankel_transform()`. Please click on the name to find out what parameters the user receives as well as their types. The parameter `strategy_name` essentially allows the user the select one of the 14 Hankel transform strategies available. This is further explained under "Note" at `hankel_transform()`.

Among the input parameters to the function, there is a C struct called `strategy_params`. This struct needs to contain different parameters depending on `strategy_name`. Clicking on `strategy_params` will give an idea of all possible fields this struct can contain. Each strategy only uses some of these, and it's important that `strategy_params` contains the right parameters required by the selected strategy. If it doesn't, an error will be thrown. 

Below is a summary of the parameters required by each strategy. 


| index (as in SASfit)      | Strategy Name       | Field names needed in `strategy_params`    |
|:--------------------------|:-------------------:|--------------------------------------------|
| 0                         | DE_Ooura            | `n_eval`, `eps_rel`                        |
| 1                         | DE_Ogata            | `n_eval`, `f_max`                          |
| 6 to 11                   | DHT_6 to DHT_11     |  N/A                                       |
| 12                        | QWE_Key             | `n_eval`, `eps_rel`                        |
| 13                        | QWE_Chave           | `n_eval`, `eps_rel`                        |
| ...                       | ...                 | ...                                        |

Here is an explanation of each parameter:
- `n_eval`: integer indicating number of function evaluations
- `eps_rel`: relative error allowed e.g. 1e-9 
- `f_max`: 	float indicating starting guess for max in form factor

Please note, for DHT strategies (that is, DHT_6 to DHT_11) feel free to supply an empty struct for `strategy_params` . Any field within it will be ignored in these strategies as they do not require such additional parameters.  

## Status codes

When successful, the `hankel_transform()` function returns 0 . When unsuccessful, it returns an error code. To find out more about the error, the table below lists all possible error codes one might receive plus corresponding explanations.

| Status codes   | Explanation                                                                   |
|:-------------|:--------------------------------------------------------------------------------|
| 0            | Success                                                                         |
| -1           | Wrong nu (order of bessel function)                                             |
| -2           | Wrong number for DHT strategy, must be int between 6 and 11                     |
| -3           | Failed to allocate variables                                                    |
| -4           | Did not converge                                                                |
| -5           | Internal error: division by zero                                                |
| -6           | Internal error: wrong nzeros in function bessel_j_zero (must be >= 1)           |
| -7           | Internal error: wrong n of iterations in pade sum (must be >= 1)                |
| -8           | Wrong strategy parameters supplied to hankel transform function                 |
| -9           | Wrong parameters supplied to form factor                                        |