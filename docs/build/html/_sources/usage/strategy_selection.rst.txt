.. _strategy-selection:

Strategy Selection
================================

LibHankel includes various implementations of the Hankel transform, corresponding to 
different ways of approximating the Hankel transform integral. 

Independently of the interface being used, the ``hankel_transform`` function is designed 
to receive a parameter named ``strategy_name``, which allows the user the select one of 
the strategies available. 

The following table lists the strategies that are available in LibHankel. 
The column named *Strategy Name* in the table essentially contains the possible strings
the user could pass in as ``strategy_name`` when calling the ``hankel_transform`` function.
The column named *Index* in the table essentially refers to the numbering used in SASfit.


.. table:: Hankel transform strategies available in LibHankel.

    +---------------------+---------------------------------------+--------------------------+
    | Strategy Name       | Type / Method                         | Index (as in SASfit)     |
    +=====================+=======================================+==========================+
    | DE_Ooura            | Double-exponential quadrature         | 0                        |
    +---------------------+---------------------------------------+--------------------------+
    | DE_Ogata            | Double-exponential quadrature         | 1                        |
    +---------------------+---------------------------------------+--------------------------+
    | DHT_6 to DHT_11     | Digital filters                       | 6 to 11                  |
    +---------------------+---------------------------------------+--------------------------+
    | QWE_Key             | QWE with continued fraction expansion | 12                       |
    +---------------------+---------------------------------------+--------------------------+
    | QWE_Chave           | QWE with Shanks transformation        | 13                       |
    +---------------------+---------------------------------------+--------------------------+


Strategies 6 to 11 in SASfit, which in this library are referred to as *DHT_6*, *DHT_7* ... 
*DHT_11*,  are **fixed-grid digital filter** strategies. They are fast and perform well on simple form
factors but their accuracy decreases on more complex or oscillatory ones. They are not iterative, 
rely on pre-computed constants, and do not allow the user to control the accuracy of the solution 
in any way. 

All other strategies instead are **iterative**, and some of them require the user to supply a tolerance 
parameter (called ``eps_rel``) which is used to determine when to stop interating. Specifically, 
iterations are stopped when the difference in the solution of two consecutive iterations is less than 
the tolerance. For this reasons, strategies that require eps_rel generally perform better than the DHT 
ones on complex form factors, though they are more slow. 
