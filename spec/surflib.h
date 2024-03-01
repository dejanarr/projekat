#ifndef __SURFLIB_H
#define __SURFLIB_H

#include "ipoint.h"
#include "fasthessian.h"
#include "surf.h"
#include "image.h"

namespace surf {

/**
 * Identify interest points and calculate their descriptor
 *
 * @param im pointer to double image
 * @param ipts (return) vector of interest points
 * @param thres blob response threshold
 * @param doubleImageSize double image size
 * @param initLobe custom lobe size
 * @param samplingStep initial sampling step
 * @param octaves number of octaves
 * @param indexSize descriptor size
 **/
inline void surfDetDes(Image *im, std::vector< Ipoint >& ipts,
        double thres = 4.0, bool doubleImageSize = false,
        int initLobe = 3, int samplingStep = 2, int octaves = 4, int indexSize = 4) {
  // Create the integral image
  Image iimage(im, doubleImageSize);

  // Extract interest points with Fast-Hessian
  FastHessian fh(&iimage, /* pointer to integral image */
         ipts, /* interest point vector to be filled */
                 thres, /* blob response threshold */
                 doubleImageSize, /* double image size flag */
                 initLobe * 3 /* 3 times lobe size equals the mask size */,
                 samplingStep, /* subsample the blob response map */
                 octaves /* number of octaves to be analysed */);

  // Extract them and get their pointer
  fh.getInterestPoints();

  // Initialise the SURF descriptor
  Surf des(&iimage, /* pointer to integral image */
           doubleImageSize, /* double image size flag */
           indexSize /* square size of the descriptor window (default 4x4)*/);

  // Compute the orientation and the descriptor for every interest point
  for (unsigned n=0; n<ipts.size(); n++){
    // set the current interest point
    des.setIpoint(&ipts[n]);
    // assign reproducible orientation
    des.assignOrientation();
    // make the SURF descriptor
    des.makeDescriptor();
  }
}

/**
 * Calculate descriptor for given interest points
 *
 * @param im pointer to double image
 * @param ipts (return) vector of interest points
 * @param doubleImageSize double image size
 * @param indexSize descriptor size
 **/
inline void surfDes(Image *im, std::vector< Ipoint >& ipts,
       bool doubleImageSize = false, int indexSize = 4) {
  // Create the integral image
  Image iimage(im, doubleImageSize);

  // Initialise the SURF descriptor
  Surf des(&iimage, /* pointer to integral image */
           doubleImageSize, /* double image size flag */
           indexSize /* square size of the descriptor window (default 4x4)*/);

  // Compute the orientation and the descriptor for every interest point
  for (unsigned n=0; n<ipts.size(); n++){
    //for (Ipoint *k = ipts; k != NULL; k = k->next){
    // set the current interest point
    des.setIpoint(&ipts[n]);
    // assign reproducible orientation
    des.assignOrientation();
    // make the SURF descriptor
    des.makeDescriptor();
  }
}

}

#endif
