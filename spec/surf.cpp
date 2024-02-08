#define SC_INCLUDE_FX
#include <math.h>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <systemc>
#include <fstream>

#include "ipoint.h"
#include "surf.h"
#include "image.h"

namespace surf {

#define OriHistTh 0.8
#define window M_PI/3
#define IndexSigma 1.0

#define get_sum(I, x1, y1, x2, y2) (I[y1+1][x1+1] + I[y2][x2] - I[y2][x1+1] - I[y1+1][x2])
#define get_wavelet1(IPatch, x, y, size) (get_sum(IPatch, x + size, y, x - size, y - size) - get_sum(IPatch, x + size, y + size, x - size, y))
#define get_wavelet2(IPatch, x, y, size) (get_sum(IPatch, x + size, y + size, x, y - size) - get_sum(IPatch, x, y + size, x - size, y - size))

#define MAX(x,y)  (((x) > (y)) ? (x) : int(y))

using namespace std;

// Constructor
Surf::Surf(){
  _iimage = NULL;
  _Pixels.clear();
  _doubleImage = false;
  _extended = false;
  _IndexSize = 4;
  _MagFactor = 3;
  _OriSize = 4;
  // calculate length of the descriptor vector
  _VecLength = _IndexSize * _IndexSize * _OriSize;
  // rotation invariance
  _upright = false;
  // allocate _index
  _index.resize(_IndexSize, std::vector<std::vector<num_f>>(
       _IndexSize, std::vector<num_f>(_OriSize, 0.0f)));
  // create _lookup tables
  createLookups();
  _sine = 0.0;
  _cose = 1.0;
}

// Constructor with parameters
Surf::Surf(Image *im, bool dbl, bool usurf,
           bool ext, int insi){
  // set image
  _iimage = im;
  //_Pixels = _iimage->getPixels();
  _doubleImage = dbl;
  
  //std::cout << "_iimage: " << im << ", _doubleImage: " << dbl << endl;

  _IndexSize = insi;
  _MagFactor = 12/insi;
  _extended = ext;
  _OriSize = 4 + _extended*4;
  // calculate length of the descriptor vector
  _VecLength = _IndexSize * _IndexSize * _OriSize;
  _upright = usurf;
  _width = _iimage->getWidth();
  _height = _iimage->getHeight();
 
  

  // create _lookup tables
  createLookups();

  double** tempPixels = _iimage->getPixels();
  _Pixels.resize(_height);
  for (int i = 0; i < _height; ++i) {
      _Pixels[i].resize(_width);
      for (int j = 0; j < _width; ++j) {
          _Pixels[i][j] = static_cast<num_f>(tempPixels[i][j]);
      }
  }

 // allocate _index
_index.resize(_IndexSize, std::vector<std::vector<num_f>>(
    _IndexSize, std::vector<num_f>(_OriSize, 0.0f)));

  // initial sine and cosine
  _sine = 0.0;
  _cose = 1.0;
}

// Destructor
Surf::~Surf() {
}

// Get length of the descriptor vector
int Surf::getVectLength(){
  return _VecLength;
}

// set Ipoint for which a descriptor has to be computed
void Surf::setIpoint(Ipoint* ipt){
  _current = ipt;
  //std::cout << "_current: " << ipt << endl;
}

// Assign orienationt
void Surf::assignOrientation() {
  

  double scale = (1.0+_doubleImage) * _current->scale;
  int x = (int)((1.0+_doubleImage) * _current->x + 0.5);
  int y = (int)((1.0+_doubleImage) * _current->y + 0.5);
  
  //std::cout << "doubleImg: " << _doubleImage << ", current: " << _current << endl;
  
  //std::cout << "y: " << y << ", x: " << x << endl;
  
  int pixSi = (int)(2*scale + 1.6);
  //std::cout << "pixSi: " << pixSi << endl;
  const int pixSi_2 = (int)(scale + 0.8);
  //std::cout << "pixSi_2: " << pixSi_2 << endl;
  double weight;
  const int radius=9;
  double dx=0, dy=0, magnitude, angle, distsq;
  const double radiussq = 81.5;
  int y1, x1;
  int yy, xx;

  vector< pair< double, double > > values;
  for (yy = y - pixSi_2*radius, y1= -radius; y1 <= radius; y1++,
       yy+=pixSi_2){
    for (xx = x - pixSi_2*radius, x1 = -radius; x1 <= radius; x1++,
         xx+=pixSi_2) {
      // Do not use last row or column, which are not valid
      if (yy + pixSi + 2 < _height &&
          xx + pixSi + 2 < _width &&
          yy - pixSi > -1 &&
          xx - pixSi > -1) {
        distsq = (y1 * y1 + x1 * x1);
        //std::cout << "x1: " << x1 << ", y1:" << y1 << endl;
        if (distsq < radiussq) {
          weight = _lookup1[(int)distsq];
          dx = get_wavelet2(_iimage->getPixels(), xx, yy, pixSi);
          dy = get_wavelet1(_iimage->getPixels(), xx, yy, pixSi);

          magnitude = sqrt(dx * dx + dy * dy);
          if (magnitude > 0.0){
            angle = atan2(dy, dx);
            values.push_back( make_pair( angle, weight*magnitude ) );
          }
        }
      }
    }
  }
  //std::cout << "angle: " << angle << ", magnitude: " << magnitude << ", distsq: " << distsq << endl;
  //std::cout << "xx: " << xx << ", yy:" << yy << ", x1: " << x1 << ", y1: " << y1 << endl;
  double best_angle = 0;

  if (values.size()) {
    sort( values.begin(), values.end() );
    int N = values.size();

    float d2Pi = 2.0*M_PI;
    //std::cout << "N: " << N << ", d2Pi: " << d2Pi <<endl;
    for( int i = 0; i < N; i++ ) {
      values.push_back( values[i] );
      values.back().first += d2Pi;
    }
    //std::cout << "values: " << values.back().first << endl;

    double part_sum = values[0].second;
    double best_sum = 0;
    double part_angle_sum = values[0].first * values[0].second;
    //part_angle_sum je negativno
    //std::cout << "part_sum: " << part_sum << ", part_angle_sum:" << part_angle_sum<< endl;
    //std::cout << "values[0].first: " << values[0].first << ", values[0].second:" << values[0].second<< endl;

    for( int i = 0, j = 0; i < N && j<2*N; ) {
      if( values[j].first - values[i].first < window ) {
        if( part_sum > best_sum ) {
          best_angle  = part_angle_sum / part_sum;
          best_sum = part_sum;
          //std::cout << "best_angle: " << best_angle << ", best_sum:" << best_sum<< endl;
        }
        j++;
        part_sum += values[j].second;
        part_angle_sum += values[j].second * values[j].first;
        //std::cout << "part_sum: " << part_sum << ", part_angle_sum:" << part_angle_sum<< endl;
      }
      else {
        part_sum -= values[i].second;
        part_angle_sum -= values[i].second * values[i].first;
        i++;
        //std::cout << "part_sum: " << part_sum << ", part_angle_sum:" << part_angle_sum<< endl;
      }
    }
  }
  _current->ori = best_angle;
  //std::cout << "best_angle: " << best_angle << endl;
}

// Compute the robust features
void Surf::makeDescriptor() {
  _current->allocIvec(_VecLength);
  //std::cout << "_VecLength: " << _VecLength << endl;
  // Initialize _index array
  for (int i = 0; i < _IndexSize; i++) {
    for (int j = 0; j < _IndexSize; j++) {
      for (int k = 0; k < _OriSize; k++)
        _index[i][j][k] = 0.0;
    }
  }

    // calculate _sine and co_sine once
    _sine = sin(_current->ori);
    _cose = cos(_current->ori);
    //std::cout << "_sine: " << _sine << ",_cose: " << _cose << endl;
    
    //std::cout << "_current " << _current << endl;
    
    // Produce _upright sample vector
    createVector(1.65*(1+_doubleImage)*_current->scale,
                 (1+_doubleImage)*_current->y,
                 (1+_doubleImage)*_current->x);
 
  int v = 0;
  for (int i = 0; i < _IndexSize; i++){
    for (int j = 0; j < _IndexSize; j++){
      for (int k = 0; k < _OriSize; k++)
        _current->ivec[v++] = _index[i][j][k];
    }
  }
  normalise();
}

// -------------------------------------------------------------------------------
// protected:
// -------------------------------------------------------------------------------

// Create Descriptor vector
void Surf::createVector(double scale, double y, double x) {

  //std::cout << "y: " << y << ", x: " << x << endl;

  int i, j, iradius, iy, ix;
  double spacing, radius, rpos, cpos, rx, cx;
  int step = MAX((int)(scale/2 + 0.5),1);
  //std::cout << "step: " << step << endl;
  //iy = static_cast<int>(y);
  //ix = static_cast<int>(x);
  iy = (int) (y + 0.5);
  ix = (int) (x + 0.5);
  
  //std::cout << "current2: " << _current << endl;
  //std::cout << "DblImg: " << _doubleImage << endl;
  
  //std::cout << "iy: " << iy << ", ix: " << ix << endl;

  double fracy = y-iy;
  double fracx = x-ix;
  double fracr =   _cose * fracy + _sine * fracx;
  double fracc = - _sine * fracy + _cose * fracx;
  //Odstupaju resenja za ova dva
  //std::cout << "fracy: " << fracy << ", fracx: " << fracx << endl;
  //std::cout << "fracc: " << fracc << ", fracr: " << fracr << endl;
  // The spacing of _index samples in terms of pixels at this scale
  spacing = scale * _MagFactor;
  //std::cout << "spacing: " << spacing << ", scale: " << scale << endl;

  // Radius of _index sample region must extend to diagonal corner of
  // _index patch plus half sample for interpolation.
  radius = 1.4 * spacing * (_IndexSize + 1) / 2.0;
  iradius = (int) (radius/step + 0.5);
  //std::cout << "radius: " << radius << ", iradius: " << iradius << endl;

  // Examine all points from the gradient image that could lie within the
  // _index square.
  for (i = -iradius; i <= iradius; i++)
    for (j = -iradius; j <= iradius; j++) {
      // Rotate sample offset to make it relative to key orientation.
      // Uses (x,y) coords.  Also, make subpixel correction as later image
      // offset must be an integer.  Divide by spacing to put in _index units.
      rpos = (step*(_cose * i + _sine * j) - fracr) / spacing;
      cpos = (step*(- _sine * i + _cose * j) - fracc) / spacing;
      //std::cout << "rpos: " << rpos << ", cpos: " << cpos << endl;
      
      // Compute location of sample in terms of real-valued _index array
      // coordinates.  Subtract 0.5 so that rx of 1.0 means to put full
      // weight on _index[1] (e.g., when rpos is 0 and _IndexSize is 3.
      rx = rpos + _IndexSize / 2.0 - 0.5;
      cx = cpos + _IndexSize / 2.0 - 0.5;
      //std::cout << "rx: " << rx << ", cx: " << cx << endl;

      // Test whether this sample falls within boundary of _index patch
      if (rx > -1.0 && rx < (double) _IndexSize  &&
          cx > -1.0 && cx < (double) _IndexSize)
        AddSample(iy + i*step, ix + j*step, rpos, cpos,
                  rx, cx, (int)(scale));
    }
}

/*void Surf::ispisiIndex() const {
    for (const auto& matrica : _index) {
        for (const auto& red : matrica) {
            for (const auto& element : red) {
                std::cout << element << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}*/

/*void Surf::ispisiPixels() const {
    for (const auto &red: _Pixels) {
        for (const auto &element: red) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}*/

// Add Sample in the descriptor vector
void Surf::AddSample(num_i r, num_i c, num_f rpos,
                     num_f cpos, num_f rx, num_f cx, num_i step) {
  num_f weight;
  num_f dx, dy;

  //this->ispisiIndex();
  //this->ispisiPixels();
  
  // Clip at image boundaries.
  if (r < 1+step  ||  r >= _height - 1-step  ||
      c < 1+step  ||  c >= _width - 1-step)
     return;

  weight = _lookup2[num_i(rpos * rpos + cpos * cpos)];
  //std::cout << "weight: " << weight << endl;
  num_f dxx, dyy;

  dxx = weight*get_wavelet2(_Pixels, c, r, step);
  dyy = weight*get_wavelet1(_Pixels, c, r, step);
  dx = _cose*dxx + _sine*dyy;
  dy = _sine*dxx - _cose*dyy;
  //std::cout << "dyy: " << dyy << ", dxx: " << dxx << endl << endl;
  //std::cout << "dy: " << dy << ", dx: " << dx << endl << endl;
  
  PlaceInIndex(dx, (dx<0?0:1), dy, (dy<0?2:3), rx, cx);

}

void Surf::PlaceInIndex(num_f mag1, num_i ori1, num_f mag2, num_i ori2, num_f rx, num_f cx) {
  // Uverite se da ste inicijalizovali _index sa pravim dimenzijama pre korišćenja.

  // Konverzija rx i cx u indekse ri i ci unutar validnih granica
  num_i ri = std::max(0, std::min(static_cast<int>(_IndexSize - 1), static_cast<int>(rx)));
  num_i ci = std::max(0, std::min(static_cast<int>(_IndexSize - 1), static_cast<int>(cx)));
  
  //std::cout<< "ri: " << ri << ", ci: " << ci << ", rx: " << rx << ", cx: " << cx << endl;

  // Izračunavanje frakcionih delova i težina
  num_f rfrac = rx - ri;
  num_f cfrac = cx - ci;
  
  rfrac = std::max(0.0f, std::min(float(rfrac), 1.0f));
  cfrac = std::max(0.0f, std::min(float(cfrac), 1.0f));
  
  num_f rweight1 = mag1 * (1.0 - rfrac);
  num_f rweight2 = mag2 * (1.0 - rfrac);
  num_f cweight1 = rweight1 * (1.0 - cfrac);
  num_f cweight2 = rweight2 * (1.0 - cfrac);
  
  //std::cout << "rfrac: " << rfrac << ", cfrac: " << cfrac << endl;
  //std::cout << "rweight1: " << rweight1 << ", rweight2: " << rweight2 << endl;
  //std::cout << "mag1: " << mag1 << ", mag2: " << mag2 << endl << endl; 
  //std::cout << "ori1: " << ori1 << ", ori2: " << ori2 << endl; 

  // Pre nego što pristupamo _index, proveravamo da li su ri i ci unutar granica
  if (ri >= 0 && ri < _IndexSize && ci >= 0 && ci < _IndexSize) {
    _index[ri][ci][ori1] += cweight1;
    _index[ri][ci][ori2] += cweight2;
  }

  // Proverite da li je ci + 1 unutar granica pre pristupa
  if (ci + 1 < _IndexSize) {
    _index[ri][ci + 1][ori1] += rweight1 * cfrac;
    _index[ri][ci + 1][ori2] += rweight2 * cfrac;
  }

  // Proverite da li je ri + 1 unutar granica pre pristupa
  if (ri + 1 < _IndexSize) {
    _index[ri + 1][ci][ori1] += mag1 * rfrac * (1.0 - cfrac);
    _index[ri + 1][ci][ori2] += mag2 * rfrac * (1.0 - cfrac);
  }
}


// Normalise descriptor vector for illumination invariance for
// Lambertian surfaces
void Surf::normalise() {
  num_f val, sqlen = 0.0, fac;
  for (num_i i = 0; i < _VecLength; i++){
    val = _current->ivec[i];
    sqlen += val * val;
  }
  fac = 1.0/sqrt(sqlen);
  for (num_i i = 0; i < _VecLength; i++)
    _current->ivec[i] *= fac;
}

// Create _lookup tables
void Surf::createLookups(){
  for (int n=0;n<83;n++)
    _lookup1[n]=exp(-((double)(n+0.5))/12.5);

  for (int n=0;n<40;n++)
    _lookup2[n]=exp(-((double)(n+0.5))/8.0);
}

}
