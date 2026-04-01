#include "mirror.h"

#include "simulation.h"

#include <algorithm>
#include <cmath>

namespace {

constexpr float kMirrorColor = 0.7f;
constexpr float kMirrorWidth = 3.0f;
constexpr float kIntersectionEpsilon = 1.0e-4f;

Vector2 mirrorCircleCenter(const MirrorSurface& mirror) {
    if (mirror.type == MirrorType::Concave) {
        return Vector2(0.0f, mirror.vertexY + mirror.radius);
    }
    return Vector2(0.0f, mirror.vertexY - mirror.radius);
}

bool isPointOnMirrorArc(const MirrorSurface& mirror, const Vector2& point) {
    if (std::fabs(point.x) > mirror.halfWidth + kIntersectionEpsilon) {
        return false;
    }

    const Vector2 center = mirrorCircleCenter(mirror);
    if (mirror.type == MirrorType::Concave) {
        return point.y <= center.y + kIntersectionEpsilon;
    }
    return point.y >= center.y - kIntersectionEpsilon;
}

Vector2 mirrorNormalAtPoint(const MirrorSurface& mirror, const Vector2& point) {
    if (mirror.type == MirrorType::Flat) {
        return Vector2(0.0f, 1.0f);
    }

    const Vector2 center = mirrorCircleCenter(mirror);
    if (mirror.type == MirrorType::Concave) {
        return (center - point).normalized();
    }
    return (point - center).normalized();
}

} // namespace

MirrorSurface makeMirrorSurface(MirrorType type, float vertexY, float halfWidth, float radius) {
    MirrorSurface mirror;
    mirror.type = type;
    mirror.vertexY = vertexY;
    mirror.halfWidth = std::max(0.1f, halfWidth);
    mirror.radius = std::max(radius, mirror.halfWidth + 0.05f);
    return mirror;
}

float clampMirrorPosition(const MirrorSurface& mirror, float x) {
    return std::clamp(x, -mirror.halfWidth, mirror.halfWidth);
}

Vector2 mirrorPointAtX(const MirrorSurface& mirror, float x) {
    const float clampedX = clampMirrorPosition(mirror, x);
    if (mirror.type == MirrorType::Flat) {
        return Vector2(clampedX, mirror.vertexY);
    }

    const Vector2 center = mirrorCircleCenter(mirror);
    const float dx = clampedX - center.x;
    const float dy = std::sqrt(std::max(0.0f, mirror.radius * mirror.radius - dx * dx));
    if (mirror.type == MirrorType::Concave) {
        return Vector2(clampedX, center.y - dy);
    }
    return Vector2(clampedX, center.y + dy);
}

SurfaceHit intersectMirror(const Vector2& origin, const Vector2& direction, const MirrorSurface& mirror) {
    if (mirror.type == MirrorType::Flat) {
        OpticalSurface flatMirror = makeFlatSurface(
            SurfaceType::FlatMirror,
            Vector2(0.0f, mirror.vertexY),
            Vector2(1.0f, 0.0f),
            Vector2(0.0f, 1.0f),
            -mirror.halfWidth,
            mirror.halfWidth,
            makeVacuumMaterial(),
            makeMirrorMaterial()
        );
        return intersectRayWithSurface(origin, direction, flatMirror);
    }

    SurfaceHit result;
    const Vector2 normalizedDirection = direction.normalized();
    const Vector2 center = mirrorCircleCenter(mirror);
    const Vector2 toOrigin = origin - center;

    const float a = dot(normalizedDirection, normalizedDirection);
    const float b = 2.0f * dot(normalizedDirection, toOrigin);
    const float c = dot(toOrigin, toOrigin) - mirror.radius * mirror.radius;
    const float discriminant = b * b - 4.0f * a * c;
    if (discriminant < 0.0f) {
        return result;
    }

    const float sqrtDiscriminant = std::sqrt(discriminant);
    const float t1 = (-b - sqrtDiscriminant) / (2.0f * a);
    const float t2 = (-b + sqrtDiscriminant) / (2.0f * a);

    const float candidates[] = { t1, t2 };
    for (float distance : candidates) {
        if (distance < 0.0f) {
            continue;
        }

        const Vector2 point = origin + normalizedDirection * distance;
        if (!isPointOnMirrorArc(mirror, point)) {
            continue;
        }

        result.hit = true;
        result.point = point;
        result.normal = orientNormalAgainstRay(normalizedDirection, mirrorNormalAtPoint(mirror, point));
        result.incomingMaterial = makeVacuumMaterial();
        result.outgoingMaterial = makeMirrorMaterial();
        result.surfaceType = SurfaceType::FlatMirror;
        result.distance = distance;
        return result;
    }

    return result;
}

void drawMirrorSurface(const MirrorSurface& mirror) {
    if (mirror.type == MirrorType::Flat) {
        draw_line(-mirror.halfWidth, mirror.vertexY, mirror.halfWidth, mirror.vertexY, kMirrorColor, kMirrorColor, kMirrorColor, kMirrorWidth);
        return;
    }

    constexpr int segmentCount = 96;
    glLineWidth(kMirrorWidth);
    glBegin(GL_LINE_STRIP);
    glColor3f(kMirrorColor, kMirrorColor, kMirrorColor);
    for (int i = 0; i <= segmentCount; ++i) {
        const float t = static_cast<float>(i) / static_cast<float>(segmentCount);
        const float x = -mirror.halfWidth + (2.0f * mirror.halfWidth * t);
        const Vector2 point = mirrorPointAtX(mirror, x);
        glVertex2f(point.x, point.y);
    }
    glEnd();
}
