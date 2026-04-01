#ifndef GEOMETRIC_OPTICS_H
#define GEOMETRIC_OPTICS_H

#include "mirror.h"

enum class LensType {
    Converging,
    Diverging
};

struct OpticalImage {
    bool valid = false;
    bool atInfinity = false;
    bool real = false;
    float distance = 0.0f;
    float height = 0.0f;
    float magnification = 0.0f;
};

OpticalImage calculateMirrorImage(MirrorType type, float objectDistance, float objectHeight, float curvatureRadius);
OpticalImage calculateLensImage(LensType type, float objectDistance, float objectHeight, float focalLength);

#endif // GEOMETRIC_OPTICS_H
