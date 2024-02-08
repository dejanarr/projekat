/*
 * Speeded-Up Robust Features (SURF)
 * https://github.com/herbertbay/SURF
 *
 * Authors: Herbert Bay, Andreas Ess, Geert Willems
 *
 * Copyright (2006): ETH Zurich, Switzerland
 * Katholieke Universiteit Leuven, Belgium
 * All rights reserved.
 *
 * For details, see the paper:
 * Herbert Bay,  Tinne Tuytelaars,  Luc Van Gool,
 *  "SURF: Speeded Up Robust Features"
 * Proceedings of the ninth European Conference on Computer Vision, May 2006
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for educational, research, and not-for-profit
 * purposes, without fee and without a signed licensing agreement, is
 * hereby granted, provided that the above copyright notice and this
 * paragraph appear in all copies modifications, and distributions.
 *
 * Any commercial use or any redistribution of this software
 * requires a license from one of the above mentioned establishments.
 *
 * For further details, contact Herbert Bay (bay@vision.ee.ethz.ch).
 */

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
