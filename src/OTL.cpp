
#include "OTL.h"

#define SINKHORN_ITERATIONS 250
#define SIGMA 1

// Finds the Gaussian filter of vector v
VectorXd Wasserstein::Gauss_filter(VectorXd v)
{
    // return v.array().exp().matrix();
    return v;
}


// The Wasserstein-2 distance using Sinkhorn's algorithm
double Wasserstein::W2(VectorXd p, VectorXd q)
{
    VectorXd v(p.size()), w(q.size());
    v.setOnes();
    w.setOnes();
    
    for (int i = 0; i < SINKHORN_ITERATIONS; i++)
    {
            v << p.cwiseQuotient(Gauss_filter(w));
            w << q.cwiseQuotient(Gauss_filter(v));
    }
    
    // // Return the actual Wasserstein distance
    return SIGMA * (
            p.dot(v.array().log10().matrix()) + 
            q.dot(w.array().log10().matrix())
        );
}



