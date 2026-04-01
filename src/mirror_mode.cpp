#include "simulation.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace {

float mirrorHalfHeight() {
    return 0.62f * scene_half_height();
}

float mirrorSurfaceXAtY(MirrorType type, float radius, float y) {
    const float clampedY = std::clamp(y, -mirrorHalfHeight(), mirrorHalfHeight());
    const float mirrorX = mirror_diagram_x();
    if (type == MirrorType::Flat) {
        return mirrorX;
    }

    const float safeRadius = std::max(radius, mirrorHalfHeight() + 0.05f);
    const float root = std::sqrt(std::max(0.0f, safeRadius * safeRadius - clampedY * clampedY));
    if (type == MirrorType::Concave) {
        const float centerX = mirrorX - safeRadius;
        return centerX + root;
    }

    const float centerX = mirrorX + safeRadius;
    return centerX - root;
}

void drawSegment(const Vector2& a, const Vector2& b, const float* color, float width, bool dashed) {
    if (dashed) {
        draw_dashed_line(a.x, a.y, b.x, b.y, color[0], color[1], color[2], width);
    } else {
        draw_line(a.x, a.y, b.x, b.y, color[0], color[1], color[2], width);
    }
}

void drawArrowDiagram(float baseX, float signedHeight, const float* color, bool dashed) {
    const float absoluteHeight = std::fabs(signedHeight);
    const float shaftHalfWidth = 0.008f * scene_half_width();
    const float headHeight = std::min(
        absoluteHeight * 0.55f,
        std::max(0.018f * scene_half_height(), absoluteHeight * 0.30f)
    );
    const float headWidth = std::min(
        0.032f * scene_half_width(),
        std::max(shaftHalfWidth * 1.8f, headHeight * 0.55f)
    );
    const float direction = (signedHeight >= 0.0f) ? 1.0f : -1.0f;
    const float headBaseY = signedHeight - direction * headHeight;
    const Vector2 tip(baseX, signedHeight);

    if (dashed) {
        draw_dashed_line(baseX, 0.0f, baseX, headBaseY, color[0], color[1], color[2], 2.4f);
        draw_dashed_line(tip.x, tip.y, tip.x - headWidth, headBaseY, color[0], color[1], color[2], 2.0f);
        draw_dashed_line(tip.x, tip.y, tip.x + headWidth, headBaseY, color[0], color[1], color[2], 2.0f);
        return;
    }

    glBegin(GL_QUADS);
    glColor3f(color[0], color[1], color[2]);
    glVertex2f(baseX - shaftHalfWidth, 0.0f);
    glVertex2f(baseX + shaftHalfWidth, 0.0f);
    glVertex2f(baseX + shaftHalfWidth, headBaseY);
    glVertex2f(baseX - shaftHalfWidth, headBaseY);
    glEnd();

    draw_filled_triangle(
        tip,
        Vector2(baseX - headWidth, headBaseY),
        Vector2(baseX + headWidth, headBaseY),
        color[0], color[1], color[2]
    );
}

void drawMirrorShape(MirrorType type, float radius) {
    const float mirrorColor[3] = { 0.12f, 0.62f, 0.92f };
    const float mirrorX = mirror_diagram_x();
    if (type == MirrorType::Flat) {
        draw_line(mirrorX, -mirrorHalfHeight(), mirrorX, mirrorHalfHeight(), mirrorColor[0], mirrorColor[1], mirrorColor[2], 5.0f);
        return;
    }

    glLineWidth(5.0f);
    glBegin(GL_LINE_STRIP);
    glColor3f(mirrorColor[0], mirrorColor[1], mirrorColor[2]);
    for (int i = 0; i <= 128; ++i) {
        const float t = static_cast<float>(i) / 128.0f;
        const float y = -mirrorHalfHeight() + (2.0f * mirrorHalfHeight() * t);
        glVertex2f(mirrorSurfaceXAtY(type, radius, y), y);
    }
    glEnd();
}

std::string focusLabelForIndex(int index) {
    return (index == 1) ? "F" : (std::to_string(index) + "F");
}

void drawAxisLabels(MirrorType type, float radius) {
    const float axisColor[3] = { 0.30f, 0.30f, 0.30f };
    const float zeroX = mirror_diagram_x();
    draw_point(zeroX, 0.0f, axisColor[0], axisColor[1], axisColor[2], 6.0f);
    draw_scene_label(Vector2(zeroX, -0.08f * scene_half_height()), "0", ImVec4(0.82f, 0.82f, 0.82f, 1.0f));

    if (type == MirrorType::Flat) {
        return;
    }

    const float step = std::max(0.12f, radius * 0.5f);
    int index = 1;
    while (true) {
        bool drewAny = false;
        const float offset = step * static_cast<float>(index);
        const float leftX = zeroX - offset;
        const float rightX = zeroX + offset;

        if (leftX > diagram_axis_min_x()) {
            draw_point(leftX, 0.0f, axisColor[0], axisColor[1], axisColor[2], 6.0f);
            draw_scene_label(Vector2(leftX, -0.08f * scene_half_height()), focusLabelForIndex(index), ImVec4(0.82f, 0.82f, 0.82f, 1.0f));
            drewAny = true;
        }
        if (rightX < diagram_axis_max_x()) {
            draw_point(rightX, 0.0f, axisColor[0], axisColor[1], axisColor[2], 6.0f);
            draw_scene_label(Vector2(rightX, -0.08f * scene_half_height()), focusLabelForIndex(index), ImVec4(0.82f, 0.82f, 0.82f, 1.0f));
            drewAny = true;
        }

        if (!drewAny) {
            break;
        }
        ++index;
    }

    const float leftRadiusX = zeroX - radius;
    const float rightRadiusX = zeroX + radius;
    if (leftRadiusX > diagram_axis_min_x()) {
        draw_scene_label(Vector2(leftRadiusX, 0.06f * scene_half_height()), "R", ImVec4(0.55f, 0.92f, 0.95f, 1.0f));
    }
    if (rightRadiusX < diagram_axis_max_x()) {
        draw_scene_label(Vector2(rightRadiusX, 0.06f * scene_half_height()), "R", ImVec4(0.55f, 0.92f, 0.95f, 1.0f));
    }
}

void drawRayToImage(const Vector2& objectTip, const Vector2& hitPoint, const Vector2& imageTip, bool realImage, const float* incomingColor, const float* outgoingColor) {
    drawSegment(objectTip, hitPoint, incomingColor, 3.0f, false);
    if (realImage) {
        const Vector2 direction = (imageTip - hitPoint).normalized();
        drawSegment(hitPoint, hitPoint + direction * (1.30f * scene_half_width()), outgoingColor, 3.0f, false);
    } else {
        const Vector2 direction = (hitPoint - imageTip).normalized();
        drawSegment(hitPoint, hitPoint + direction * (1.30f * scene_half_width()), outgoingColor, 3.0f, false);
        drawSegment(hitPoint, imageTip, outgoingColor, 2.0f, true);
    }
}

void drawParallelMirrorRay(const Vector2& objectTip, const Vector2& imageTip, const OpticalImage& image) {
    const Vector2 hitPoint(mirrorSurfaceXAtY(currentMirrorType, mirrorCurvatureRadius, objectTip.y), objectTip.y);
    drawRayToImage(objectTip, hitPoint, imageTip, image.real, rayAColor, rayAColor);
}

void drawVertexMirrorRay(const Vector2& objectTip, const Vector2& imageTip, const OpticalImage& image) {
    const Vector2 vertex(mirror_diagram_x(), 0.0f);
    drawRayToImage(objectTip, vertex, imageTip, image.real, rayCColor, rayCColor);
}

void drawInfinityMirrorRays(const Vector2& objectTip) {
    const Vector2 parallelHit(mirrorSurfaceXAtY(currentMirrorType, mirrorCurvatureRadius, objectTip.y), objectTip.y);
    drawSegment(objectTip, parallelHit, rayAColor, 3.0f, false);
    drawSegment(parallelHit, Vector2(diagram_axis_min_x(), objectTip.y), rayAColor, 3.0f, false);

    const Vector2 vertex(mirror_diagram_x(), 0.0f);
    drawSegment(objectTip, vertex, rayCColor, 3.0f, false);
    const Vector2 direction = (vertex - objectTip).normalized();
    drawSegment(vertex, vertex + direction * (0.90f * scene_half_width()), rayCColor, 3.0f, false);
}

} // namespace

void draw_mirror_mode() {
    objectHeight = std::clamp(objectHeight, 0.10f, diagram_object_height_limit());
    objectDistance = clamp_mirror_object_distance(objectDistance);

    const float signedObjectHeight = objectPointsUp ? objectHeight : -objectHeight;
    const float objectBaseX = mirror_diagram_x() - objectDistance;
    const OpticalImage image = calculateMirrorImage(currentMirrorType, objectDistance, signedObjectHeight, mirrorCurvatureRadius);

    draw_line(diagram_axis_min_x(), 0.0f, diagram_axis_max_x(), 0.0f, 0.28f, 0.28f, 0.28f, 2.0f);
    drawMirrorShape(currentMirrorType, mirrorCurvatureRadius);
    drawAxisLabels(currentMirrorType, mirrorCurvatureRadius);

    const float objectColor[3] = { 0.92f, 0.92f, 0.94f };
    drawArrowDiagram(objectBaseX, signedObjectHeight, objectColor, false);

    if (image.valid && !image.atInfinity) {
        const float imageBaseX = mirror_diagram_x() - image.distance;
        const float imageColor[3] = { rayBColor[0], rayBColor[1], rayBColor[2] };
        drawArrowDiagram(imageBaseX, image.height, imageColor, !image.real);

        if (showConstructionRays) {
            const Vector2 objectTip(objectBaseX, signedObjectHeight);
            const Vector2 imageTip(imageBaseX, image.height);
            drawParallelMirrorRay(objectTip, imageTip, image);
            drawVertexMirrorRay(objectTip, imageTip, image);
        }
    } else if (showConstructionRays) {
        drawInfinityMirrorRays(Vector2(objectBaseX, signedObjectHeight));
    }
}
