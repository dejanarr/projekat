#ifndef IPOINT_H
#define IPOINT_H

#include <stdlib.h>

namespace surf {

class Ipoint {
  public:
  // Constructor
  Ipoint(){
    ivec = NULL;
    ori = 0.0;
  };

  // Destructor
  ~Ipoint(){
    if (ivec)
      delete [] ivec;
  };

  // Allocate space
  void allocIvec(const int si){
    ivec = new double[si];
  };

    // x, y value of the interest point
    double x, y;
    // detected scale
    double scale;
    // strength of the interest point
    double strength;
    // orientation
    double ori;
    // sign of Laplacian
    int laplace;
    // descriptor
    double *ivec;
};

}

#endif // IPOINT_H
