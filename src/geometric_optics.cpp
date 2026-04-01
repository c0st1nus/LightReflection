#include "geometric_optics.h"

#include <cmath>

namespace {

constexpr float kEpsilon = 1.0e-4f;

} // namespace

OpticalImage calculateMirrorImage(MirrorType type, float objectDistance, float objectHeight, float curvatureRadius) {
    OpticalImage image;
    if (objectDistance <= kEpsilon) {
        return image;
    }

    image.valid = true;
    if (type == MirrorType::Flat) {
        image.distance = -objectDistance;
        image.height = objectHeight;
        image.magnification = 1.0f;
        image.real = false;
        return image;
    }

    float focalLength = 0.5f * curvatureRadius;
    if (type == MirrorType::Convex) {
        focalLength *= -1.0f;
    }

    const float denominator = objectDistance - focalLength;
    if (std::fabs(denominator) <= kEpsilon) {
        image.atInfinity = true;
        image.real = focalLength > 0.0f;
        return image;
    }

    image.distance = (focalLength * objectDistance) / denominator;
    image.magnification = -image.distance / objectDistance;
    image.height = image.magnification * objectHeight;
    image.real = image.distance > 0.0f;
    return image;
}

OpticalImage calculateLensImage(LensType type, float objectDistance, float objectHeight, float focalLength) {
    OpticalImage image;
    if (objectDistance <= kEpsilon || focalLength <= kEpsilon) {
        return image;
    }

    image.valid = true;
    const float signedFocalLength = (type == LensType::Converging) ? focalLength : -focalLength;
    const float denominator = objectDistance - signedFocalLength;
    if (std::fabs(denominator) <= kEpsilon) {
        image.atInfinity = true;
        image.real = signedFocalLength > 0.0f;
        return image;
    }

    image.distance = (signedFocalLength * objectDistance) / denominator;
    image.magnification = -image.distance / objectDistance;
    image.height = image.magnification * objectHeight;
    image.real = image.distance > 0.0f;
    return image;
}
