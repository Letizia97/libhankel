
#ifndef FORM_FACTORS_H
#define FORM_FACTORS_H


/*
Compute the g_dab form factor. 

Receives: 
    - double q          value at which to compute the Hankel t
    - double (*params)  pointer to an array of params for the function

    Params must contain :
        -   XI
        -   H
        -   ETA
*/
double form_factor_g_dab(double q, double (*params)[50]); 



/*
Compute the sphere form factor. 

Receives: 
    - double q          value at which to compute the Hankel t
    - double (*params)  pointer to an array of params for the function
    
    Params must contain :
        -   R
        -   ETA
*/
double form_factor_sphere(double q, double (*params)[50]);



/*
Compute the broad peak form factor. 

Receives: 
    - double q          value at which to compute the Hankel t
    - double (*params)  pointer to an array of params for the function
    Params must contain :
        -   I0
        -   XI
        -   Q0
        -   M
        -   P
*/
double form_factor_broad_peak(double q, double (*params)[50]);


#endif // FORM_FACTORS_H
