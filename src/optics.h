#ifndef OPTICS_H
#define OPTICS_H

#include "vector2.h"

enum class MaterialType {
    Vacuum,
    Mirror,
    Dielectric
};

enum class SurfaceType {
    FlatMirror,
    FlatInterface,
    LensSurface,
    PrismFace
};

struct OpticalMaterial {
    MaterialType type = MaterialType::Vacuum;
    float refractiveIndex = 1.0f;
};

struct SurfaceHit {
    bool hit = false;
    Vector2 point;
    Vector2 normal;
    OpticalMaterial incomingMaterial;
    OpticalMaterial outgoingMaterial;
    SurfaceType surfaceType = SurfaceType::FlatInterface;
    float distance = 0.0f;
};

struct OpticalSurface {
    SurfaceType type = SurfaceType::FlatInterface;
    Vector2 point;
    Vector2 tangent = Vector2(1.0f, 0.0f);
    Vector2 normal = Vector2(0.0f, 1.0f);
    float minProjection = -1.0f;
    float maxProjection = 1.0f;
    OpticalMaterial positiveSideMaterial;
    OpticalMaterial negativeSideMaterial;
};

struct RefractionResult {
    bool totalInternalReflection = false;
    Vector2 direction;
};

OpticalMaterial makeVacuumMaterial();
OpticalMaterial makeMirrorMaterial();
OpticalMaterial makeDielectricMaterial(float refractiveIndex);

OpticalSurface makeFlatSurface(
    SurfaceType type,
    const Vector2& point,
    const Vector2& tangent,
    const Vector2& normal,
    float minProjection,
    float maxProjection,
    const OpticalMaterial& positiveSideMaterial,
    const OpticalMaterial& negativeSideMaterial
);

SurfaceHit intersectRayWithSurface(const Vector2& origin, const Vector2& direction, const OpticalSurface& surface);
Vector2 orientNormalAgainstRay(const Vector2& direction, const Vector2& normal);
Vector2 reflectDirection(const Vector2& incident, const Vector2& normal);
RefractionResult refractDirection(const Vector2& incident, const Vector2& normal, float n1, float n2);
float incidentAngleDegrees(const Vector2& incident, const Vector2& normal);

#endif // OPTICS_H
