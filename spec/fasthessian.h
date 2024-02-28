#ifndef FASTHESSIAN_H
#define FASTHESSIAN_H

#include <vector>

namespace surf {

class Ipoint;
class Image;

class FastHessian {
  //private:
    // Constructor
  //  FastHessian();

  public:
    //! Destructor
    ~FastHessian();

    //! Constructor with parameters
    FastHessian(Image *im, std::vector< Ipoint >& ip, double thres = 0.2, bool doub = false,
                short int initMasksize = 9, short int samplingStep = 2,
                short int octaves=4);

    //! Pass the integral image
    void setIimage( Image *iim );

    //! Detect the interest Points, write into ipts
    void getInterestPoints();

    //! Create a new ipoint at location (x, y),  at a certain scale
    //! and corner response strength
    void makeIpoint(double x, double y, double scale, double strength=0);

  protected:
    //! Allocate scale layers for one octave
    void allocateOctave();

    //! Fast non-maximum-suppression
    void findMaximum(int *borders, int o, int octave);
    void interpFeature(int s, int row, int col, Image *map,
                       int o, int octave, int movesRemain,
                       int *borders);
    double fitQuadrat(int s, int r, int c);

  private:
    //! integral image
    Image *_Iimage;
    //! octave
    Image **_scaleLevel;
    //! vector of variables
    int _vas[9];
    //! threshold for interest point detection
    double _threshold;
    //! boolean indicating whether the image size was doubled or not
    //! default is false
    bool _doubled;
    //! reference to vector of interest points passed from outside
    std::vector< Ipoint >& _ipts;
    //! initial lobe size for the second derivative in one direction
    //! default is 3
    short int _initLobe;
    //! number scales
    short int _maxScales;
    //! number octaves
    short int _maxOctaves;
    //! the sampling step
    short int _sampling;
    //! Integral image dimensions
    int _width;
    int _height;
    //! Result of fitting quadratic
    double _offset[3];
};

}

#endif // FASTHESSIAN_H
