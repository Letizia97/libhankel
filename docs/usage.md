# Usage Guide {#usage}

## Status codes


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