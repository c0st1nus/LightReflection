#include "optics.h"

#include <algorithm>
#include <cmath>

namespace {

constexpr float kParallelEpsilon = 1.0e-5f;
constexpr float kIntersectionEpsilon = 1.0e-4f;

} // namespace

OpticalMaterial makeVacuumMaterial() {
    return { MaterialType::Vacuum, 1.0f };
}

OpticalMaterial makeMirrorMaterial() {
    return { MaterialType::Mirror, 1.0f };
}

OpticalMaterial makeDielectricMaterial(float refractiveIndex) {
    return { MaterialType::Dielectric, refractiveIndex };
}

OpticalSurface makeFlatSurface(
    SurfaceType type,
    const Vector2& point,
    const Vector2& tangent,
    const Vector2& normal,
    float minProjection,
    float maxProjection,
    const OpticalMaterial& positiveSideMaterial,
    const OpticalMaterial& negativeSideMaterial
) {
    OpticalSurface surface;
    surface.type = type;
    surface.point = point;
    surface.tangent = tangent.normalized();
    surface.normal = normal.normalized();
    surface.minProjection = std::min(minProjection, maxProjection);
    surface.maxProjection = std::max(minProjection, maxProjection);
    surface.positiveSideMaterial = positiveSideMaterial;
    surface.negativeSideMaterial = negativeSideMaterial;
    return surface;
}

SurfaceHit intersectRayWithSurface(const Vector2& origin, const Vector2& direction, const OpticalSurface& surface) {
    SurfaceHit result;
    const Vector2 normalizedDirection = direction.normalized();
    const float denominator = dot(normalizedDirection, surface.normal);
    if (std::fabs(denominator) < kParallelEpsilon) {
        return result;
    }

    const float distance = dot(surface.point - origin, surface.normal) / denominator;
    if (distance < 0.0f) {
        return result;
    }

    const Vector2 point = origin + normalizedDirection * distance;
    const float projection = dot(point - surface.point, surface.tangent);
    if (projection < surface.minProjection - kIntersectionEpsilon || projection > surface.maxProjection + kIntersectionEpsilon) {
        return result;
    }

    result.hit = true;
    result.point = point;
    result.normal = denominator < 0.0f ? surface.normal : surface.normal * -1.0f;
    result.incomingMaterial = denominator < 0.0f ? surface.positiveSideMaterial : surface.negativeSideMaterial;
    result.outgoingMaterial = denominator < 0.0f ? surface.negativeSideMaterial : surface.positiveSideMaterial;
    result.surfaceType = surface.type;
    result.distance = distance;
    return result;
}

Vector2 orientNormalAgainstRay(const Vector2& direction, const Vector2& normal) {
    const Vector2 normalizedDirection = direction.normalized();
    const Vector2 normalizedNormal = normal.normalized();
    return dot(normalizedDirection, normalizedNormal) <= 0.0f ? normalizedNormal : normalizedNormal * -1.0f;
}

Vector2 reflectDirection(const Vector2& incident, const Vector2& normal) {
    const Vector2 normalizedIncident = incident.normalized();
    const Vector2 orientedNormal = orientNormalAgainstRay(normalizedIncident, normal);
    return (normalizedIncident - 2.0f * dot(normalizedIncident, orientedNormal) * orientedNormal).normalized();
}

RefractionResult refractDirection(const Vector2& incident, const Vector2& normal, float n1, float n2) {
    RefractionResult result;
    const Vector2 normalizedIncident = incident.normalized();
    const Vector2 orientedNormal = orientNormalAgainstRay(normalizedIncident, normal);

    const float eta = n1 / n2;
    const float cosI = -dot(normalizedIncident, orientedNormal);
    const float sinT2 = eta * eta * std::max(0.0f, 1.0f - cosI * cosI);
    if (sinT2 > 1.0f) {
        result.totalInternalReflection = true;
        return result;
    }

    const float cosT = std::sqrt(std::max(0.0f, 1.0f - sinT2));
    result.direction = (eta * normalizedIncident + (eta * cosI - cosT) * orientedNormal).normalized();
    return result;
}

float incidentAngleDegrees(const Vector2& incident, const Vector2& normal) {
    const Vector2 normalizedIncident = incident.normalized();
    const Vector2 orientedNormal = orientNormalAgainstRay(normalizedIncident, normal);
    const float cosine = std::clamp(-dot(normalizedIncident, orientedNormal), -1.0f, 1.0f);
    return std::acos(cosine) * 180.0f / static_cast<float>(M_PI);
}
