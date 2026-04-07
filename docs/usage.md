# Usage Guide {#usage}

## Main function

The main function a user will be using is `hankel_transform()`. Please click on the name to find out what parameters the user receives as well as their types. The parameter `strategy_name` essentially allows the user the select one of the 14 Hankel transform strategies available. This is further explained under "Note" at `hankel_transform()`.


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