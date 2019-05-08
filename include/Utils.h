
#include "CImg.h"
#include <cassert>


class Utils{
    public:
        static void printImageAsMatrix(CImg<double> x)
        {
            assert(x.depth() == 1);
            //assert(x.spectrum() == 1);
            for (int i = 0; i < x.height(); i++)
            {
                for (int j = 0; j < x.width(); j++)
                {
                    printf("%f ", x(j,i));
                }
                printf("\n");
            }
        }

        static CImg<double> makeDistribution(CImg<double> x)
        {
            x.blur(1.0, false, true);
            return x / x.sum();
        }

        /*static CImgList<double> boxes(int width, int height)
        {
            CImgList<double> out_list;

            int w2 = width/2;
            int h2 = height/2;

            int rngs [][] = {{0,0}, {w2,0}, {0,h2}, {w2,h2}};
            for (int i = 0; i < 4; i++) {
                
            }

            return out_list;
        }*/
        
};
