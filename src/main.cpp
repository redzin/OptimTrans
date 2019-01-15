
#include <iostream>
#include "OTL.h"
#include "Eigen/Dense"

using namespace Eigen;

int main()
{
    VectorXd p(6), q(6), r(6);
    p << 0.1,0.5,0.1,0.1,0.1,0.2;
    q << 0.1,0.1,0.1,0.5,0.1,0.1;
    r << 0.1,0.1,0.1,0.1,0.1,0.5;

    MatrixXd m(6,3);
    m << p,q,r;

    double w = OTL::W2(p,q);

    std::cout << w << std::endl;
    
}


