#ifndef SURF_H
#define SURF_H

#define SC_INCLUDE_FX
#include <systemc>
#include <vector>

typedef sc_dt::sc_int<11> num_i;
typedef sc_dt::sc_fixed<48, 30, sc_dt::SC_TRN, sc_dt::SC_SAT> num_f;


namespace surf {

class Image;

class Surf {  
  public:
    //! Constructor
    Surf();

    //! Constructor with parameters
    Surf(Image *im, bool dbl=false, int insi=4);

    //! Destructor
    ~Surf();

    //! Get length of the descriptor vector
    int getVectLength();

    //! set Ipoint for which a descriptor has to be computed
    void setIpoint(Ipoint *ipt);

    //! Assign reproducible orienation
    void assignOrientation();

    //! Compute the robust features
    void makeDescriptor();
    
    std::vector<double> getRValues() const {return rValues; }
    
    //void ispisiIndex() const; 
    
    //void ispisiPixels() const;

  protected:
    //! Create the vector
    void createVector(double scale,
                      double row, double col);

    //! Add sample to the vector
    void AddSample(num_i r, num_i c, num_f rpos,
                   num_f cpos, num_f rx, num_f cx, num_i step);

    //! Place sample to index in vector
    void PlaceInIndex(num_f mag1, num_i ori1,
                      num_f mag2, num_i ori2, num_f rx, num_f cx);

    //! Normalise descriptor vector for illumination invariance for
    //! Lambertian surfaces
    void normalise();

    //! Create Lookup tables
    void createLookups(); 


  private:
    Image *_iimage;
    Ipoint *_current;
    //num_f ***_index;
    std::vector<std::vector<std::vector<num_f>>> _index;
    bool _doubleImage;
    num_i _VecLength;
    num_i _IndexSize;
    num_i _MagFactor;
    num_i _OriSize;
    num_i _width, _height;

    num_f _sine, _cose;
    std::vector<std::vector<num_f>> _Pixels;

    num_f _lookup1[83], _lookup2[40];
    
    std::vector<double> rValues;
};

}

#endif // SURF_H
