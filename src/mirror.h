#ifndef MIRROR_H
#define MIRROR_H

#include "optics.h"

enum class MirrorType {
    Flat,
    Concave,
    Convex
};

struct MirrorSurface {
    MirrorType type = MirrorType::Flat;
    float vertexY = -0.8f;
    float halfWidth = 0.8f;
    float radius = 1.2f;
};

MirrorSurface makeMirrorSurface(MirrorType type, float vertexY, float halfWidth, float radius);
float clampMirrorPosition(const MirrorSurface& mirror, float x);
Vector2 mirrorPointAtX(const MirrorSurface& mirror, float x);
SurfaceHit intersectMirror(const Vector2& origin, const Vector2& direction, const MirrorSurface& mirror);
void drawMirrorSurface(const MirrorSurface& mirror);

#endif // MIRROR_H
