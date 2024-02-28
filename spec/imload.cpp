#include <opencv2/opencv.hpp>
#include "image.h"
#include "imload.h"
using namespace cv;
using namespace std;

namespace surf {

Image *ImLoad::readImage(const char *fn) {
    // Učitavanje slike koristeći OpenCV

    string path = "../data/";
    path += fn;
    Mat colorImg = imread(path, IMREAD_COLOR);
    if (colorImg.empty()) {
        cerr << "Sorry, could not open or find the image: " << path << endl;
        exit(0);
    }

    // Konverzija slike u grayscale
    Mat grayImg;
    cvtColor(colorImg, grayImg, COLOR_BGR2GRAY);

    // Kreiranje Image objekta kompatibilnog sa SURF
    Image *im = new Image(grayImg.cols, grayImg.rows);
    for (int y = 0; y < grayImg.rows; y++) {
        for (int x = 0; x < grayImg.cols; x++) {
            im->setPix(x, y, (double)grayImg.at<uchar>(y, x) / 255.0);
        }
    }

    return im;
}

}
