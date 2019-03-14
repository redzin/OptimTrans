#include "OTL.h"
#include "CImg.h"
#include "Utils.h"
#include <iostream>
#include <cassert>

#define SINKHORN_ITERATIONS 5
#define SINKHORN_VERBOSE false
#define SIGMA 1.0
#define SMALLEST_NON_ZERO_VALUE pow(10.0, -300.0)

using namespace cimg_library;

// L2-norm (Euclidean distance)
double OTL::D2(CImg<> p, CImg<> q)
{
    return (p-q).magnitude();
}


// Entropically Smoothed Wasserstein Distance
// Computed using Sinkhorn's algorithm
double OTL::W2_Sinkhorn(CImg<> p, CImg<> q, bool useDeriche)
{

    assert(p.width() == q.width());
    assert(p.height() == q.height());
    assert(p.depth() == q.depth());
    assert(p.spectrum() == q.spectrum());
    
    CImg<double> v(p.width(), p.height(), p.depth(), p.spectrum(),1.0);
    CImg<double> w(q.width(), q.height(), q.depth(), q.spectrum(),1.0);
    
    // std::cout << p.sum() << std::endl;
    // std::cout << q.sum() << std::endl;
    std::cout << "----------------------" << std::endl;
    
    if (!useDeriche){
        for (int i = 0; i < SINKHORN_ITERATIONS; i++)
        {
            // No-boundary (Dirichlet) condition Van Vliet Gaussian filter
            v = p.get_div(w.get_blur(SIGMA, false, true).get_max(SMALLEST_NON_ZERO_VALUE));
            w = q.get_div(v.get_blur(SIGMA, false, true).get_max(SMALLEST_NON_ZERO_VALUE));
            //std::cout << v.max() << std::endl;
            //Utils::printImageAsMatrix(p.get_mul(v.max(SMALLEST_NON_ZERO_VALUE).log()));
            // CImgDisplay disp(w);
            // while (!disp.is_closed()) {
            //     disp.wait();
            // }
        }
    } else {
        for (int i = 0; i < SINKHORN_ITERATIONS; i++)
        {
            // No-boundary (Dirichlet) condition Deriche approximate Gaussian filter
            v = p.get_div(w.get_blur(SIGMA, SIGMA, SIGMA, false, false).get_max(SMALLEST_NON_ZERO_VALUE));
            w = q.get_div(v.get_blur(SIGMA, SIGMA, SIGMA, false, false).get_max(SMALLEST_NON_ZERO_VALUE));
        }
    }
    return (p.get_mul(v.max(SMALLEST_NON_ZERO_VALUE).log()) + q.get_mul(w.max(SMALLEST_NON_ZERO_VALUE).log())).sum() * SIGMA;

}

CImg<double> OTL::SinkhornEvalR(CImg<> v, CImg<> w, CImg<> a)
{
    return v.get_mul((w.get_mul(a)).blur(SIGMA, false, true));
}

CImg<double> OTL::SinkhornEvalL(CImg<> v, CImg<> w, CImg<> a)
{
    return w.get_mul((v.get_mul(a)).blur(SIGMA, false, true));
}

// CImg<double> OTL::QuadrantTest(CImg<> im_1, CImg<> im_2)
// {
//     double W2 = OTL::W2_Sinkhorn(im_1, im_2);

// }



// Wasserstein Distance Computed Directly
// double OTL::W2_direct(CImg<> p, CImg<> q)
// {
    // CImg<double> d2(6,6,1,1, 0.0);
    // CImg<double> A_eq(6*6,2*6,1,1, 0.0);
    // CImg<double> b_eq(p+p);

    // for(int i = 0; i < 6; i++)
    // {
    //     for (int j = 0; j < 6; j++)
    //     {
    //         d2(j,i) = (i-j)*(i-j);
    //         A_eq(j,j*6+i) = 1.0;
    //         A_eq(6+j,i*6+j) = 1.0;
    //     }
    // }
    
    // CImg<double> c = d2.get_resize(6*6);
    
    // printf("A shape: %d, %d", A_eq.height(), A_eq.width());

    // CImg<double> x = b_eq.get_solve(A_eq);

    // return c.dot(x);

//     return 0.0;
// }
