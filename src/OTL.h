
#include "Eigen/Dense"

using namespace Eigen;

class OTL {
    public:
        OTL();
        static double W2(VectorXd, VectorXd);
    private:
        static VectorXd Gauss_filter(VectorXd);
};
