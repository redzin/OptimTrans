
#include <iostream>
#include "CImg.h"
#include "OTL.h"
#include "Utils.h"
#include <cassert>
//#include "Eigen/Dense"

using namespace std;
using namespace cimg_library;
//using namespace Eigen;


int main(int argc, char *argv[])
{
    

    if (argc == 3)
    {
        CImg<double> image_01(argv[1]), image_02(argv[2]);
        image_01 = Utils::makeDistribution(image_01);
        image_02 = Utils::makeDistribution(image_02);
        int combined_width = image_01.width() + image_02.width();
        int combined_height = max(image_02.height(), image_01.height());
        CImg<double> image_combiner(combined_width, combined_height, 1, 3, 0);
        image_combiner.draw_image(0, 0, image_01);
        image_combiner.draw_image(image_01.width(), 0, image_02);
        double d2 = OTL::D2(image_01 / image_01.sum(), image_02 / image_02.sum());
        double w2 = OTL::W2_Sinkhorn(image_01 / image_01.sum(), image_02 / image_02.sum());
        char title[200];
        sprintf(title, "W2 = %5.4f, D2 = %5.4f", w2, d2);
        CImgDisplay main_disp(image_combiner, title);
        while (!main_disp.is_closed()) {
            main_disp.wait();
        }
    } else if (std::string(argv[1]) == "blobs") {
        CImg<double>    image_01("images/one-blob.png"),
                        image_02("images/one-blob-moved.png"),
                        image_03("images/one-blob-moved-even-more.png"),
                        image_04("images/one-blob-moved-even-more-again.png");
        
        image_01 = Utils::makeDistribution(image_01);
        image_02 = Utils::makeDistribution(image_02);
        image_03 = Utils::makeDistribution(image_03);
        image_04 = Utils::makeDistribution(image_04);

        int combined_width = image_01.width() + image_02.width() + image_03.width() + image_04.width();
        int combined_height = max(max(max(image_01.height(), image_02.height()), image_03.height()), image_04.height());
        CImg<double> image_combiner(combined_width, combined_height, 1, 3, 0);
        image_combiner.draw_image(0, 0, image_01);
        image_combiner.draw_image(image_01.width(), 0, image_02);
        image_combiner.draw_image(image_01.width()+image_02.width(), 0, image_03);
        image_combiner.draw_image(image_01.width()+image_02.width()+image_03.width(), 0, image_04);
        double d2_2 = OTL::D2(image_01, image_02);
        double d2_3 = OTL::D2(image_01, image_03);
        double d2_4 = OTL::D2(image_01, image_04);
        double w2_2 = OTL::W2_Sinkhorn(image_01, image_02);
        double w2_3 = OTL::W2_Sinkhorn(image_01, image_03);
        double w2_4 = OTL::W2_Sinkhorn(image_01, image_04);
        cout << "1 and 2: W2 = " << w2_2 << ", D2 = " << d2_2 << endl;
        cout << "1 and 3: W2 = " << w2_3 << ", D2 = " << d2_3 << endl;
        cout << "1 and 4: W2 = " << w2_4 << ", D2 = " << d2_4 << endl;
        char title[200];
        sprintf(title, "Blobs");
        //sprintf(title, "W2 = %5.4f, D2 = %5.4f", w2, d2);
        CImgDisplay main_disp(image_combiner, title);
        while (!main_disp.is_closed()) {
            main_disp.wait();
        }
    

    } else if (std::string(argv[1]) == "smile"){
        
        CImg<double>    image_01("images/smile.png"),
                        image_02("images/dots.png");
        
        image_01 = Utils::makeDistribution(image_01);
        image_02 = Utils::makeDistribution(image_02);
        
        int combined_width = image_01.width() + image_02.width();
        int combined_height = max(image_01.height(), image_02.height());
        CImg<double> image_combiner(combined_width, combined_height, 1, 3, 0);
        image_combiner.draw_image(0, 0, image_02);
        image_combiner.draw_image(image_02.width(), 0, image_01);
        double d2_2 = OTL::D2(image_01, image_02);
        double w2_2 = OTL::W2_Sinkhorn(image_01, image_02);
        cout << "1 and 2: W2 = " << w2_2 << ", D2 = " << d2_2 << endl;
        char title[200];
        sprintf(title, "Smile");
        //sprintf(title, "W2 = %5.4f, D2 = %5.4f", w2, d2);
        CImgDisplay main_disp(image_combiner, title);
        while (!main_disp.is_closed()) {
            main_disp.wait();
        }
    
        

    } else {
        CImg<double>    p(6, 1, 1, 1, 0.1, 0.5, 0.1, 0.1, 0.1, 0.1),
                        q(6, 1, 1, 1, 0.1, 0.1, 0.1, 0.5, 0.1, 0.1),
                        r(6, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.5);
        
        printf("D2(p,p) = %5.2f, W2(p,p) = %5.2f\n", OTL::D2(p,p), OTL::W2_Sinkhorn(p,p));
        printf("D2(p,q) = %5.2f, W2(p,q) = %5.2f\n", OTL::D2(p,q), OTL::W2_Sinkhorn(p,q));
        printf("D2(p,r) = %5.2f, W2(p,r) = %5.2f\n", OTL::D2(p,r), OTL::W2_Sinkhorn(p,r));
        

        // Test of the Gaussian
        // printf("------------------\n");
        // Utils::printImageAsMatrix(p);
        // Utils::printImageAsMatrix(p.get_vanvliet(1.0, 0, 'x', false));
        // Utils::printImageAsMatrix(p.get_vanvliet(1.0, 0, 'x', true));
        // Utils::printImageAsMatrix(p.get_deriche (1.0, 0, 'x', false));
        // Utils::printImageAsMatrix(p.get_deriche (1.0, 0, 'x', true));
        // Utils::printImageAsMatrix(p);

    }
    
    return 0;
    
}



