
#include "CImg.h"

using namespace cimg_library;

//template <class T>
class OTL
{
    public:
        OTL();
        static double D2(CImg<> p, CImg<> q);
        // static double W2_direct(CImg<> p, CImg<> q);
        static double W2_Sinkhorn(CImg<> p, CImg<> q, bool useDeriche = false);
        static CImg<double> SinkhornEvalR(CImg<> v, CImg<> w, CImg<> a);
        static CImg<double> SinkhornEvalL(CImg<> v, CImg<> w, CImg<> a);
        static CImg<double> QuadrantTest(CImg<> im_1, CImg<> im_2);
};

/*
class OTLEigen
{
    public:
        OTLEigen();
        static double W2(MatrixXd p, MatrixXd q);
        static MatrixXd Gauss_filter(MatrixXd v);
        static MatrixXd stdNestedArrayToEigen(array2d<double> v);
};
*/

