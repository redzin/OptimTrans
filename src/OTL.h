
#include "Eigen/Dense"

using namespace Eigen;

class Wasserstein {
    public:
        Wasserstein();
        static double W2(VectorXd, VectorXd);
    private:
        static VectorXd Gauss_filter(VectorXd);
};
