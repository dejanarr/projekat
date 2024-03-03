#ifndef IMLOAD_H
#define IMLOAD_H

#include <opencv2/opencv.hpp>
#include "image.h"

namespace surf {

class ImLoad {
public:
    static Image *readImage(const char *filename);
};

} // namespace surf

#endif // IMLOAD_H

