#include "simulation.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace {

float lensHalfHeight() {
    return 0.66f * scene_half_height();
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

void drawLensShape() {
    const float lensColor[3] = { 0.12f, 0.62f, 0.92f };
    const float lensX = lens_diagram_x();
    glLineWidth(4.5f);

    glBegin(GL_LINE_STRIP);
    glColor3f(lensColor[0], lensColor[1], lensColor[2]);
    for (int i = 0; i <= 128; ++i) {
        const float t = -1.0f + 2.0f * static_cast<float>(i) / 128.0f;
        const float y = t * lensHalfHeight();
        float x = lensX;
        if (currentLensType == LensType::Converging) {
            x = lensX - 0.030f * scene_half_width() - 0.040f * scene_half_width() * (1.0f - t * t);
        } else {
            x = lensX - 0.075f * scene_half_width() + 0.040f * scene_half_width() * (1.0f - t * t);
        }
        glVertex2f(x, y);
    }
    glEnd();

    glBegin(GL_LINE_STRIP);
    glColor3f(lensColor[0], lensColor[1], lensColor[2]);
    for (int i = 0; i <= 128; ++i) {
        const float t = -1.0f + 2.0f * static_cast<float>(i) / 128.0f;
        const float y = t * lensHalfHeight();
        float x = lensX;
        if (currentLensType == LensType::Converging) {
            x = lensX + 0.030f * scene_half_width() + 0.040f * scene_half_width() * (1.0f - t * t);
        } else {
            x = lensX + 0.075f * scene_half_width() - 0.040f * scene_half_width() * (1.0f - t * t);
        }
        glVertex2f(x, y);
    }
    glEnd();
}

std::string focusLabelForIndex(int index) {
    return (index == 1) ? "F" : (std::to_string(index) + "F");
}

void drawAxisLabels() {
    const float axisColor[3] = { 0.30f, 0.30f, 0.30f };
    const float zeroX = lens_diagram_x();
    draw_point(zeroX, 0.0f, axisColor[0], axisColor[1], axisColor[2], 6.0f);
    draw_scene_label(Vector2(zeroX, -0.08f * scene_half_height()), "0", ImVec4(0.82f, 0.82f, 0.82f, 1.0f));

    int index = 1;
    while (true) {
        bool drewAny = false;
        const float offset = lensFocalLength * static_cast<float>(index);
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
}

void drawRayToImage(const Vector2& objectTip, const Vector2& hitPoint, const Vector2& imageTip, bool realImage, const float* color) {
    drawSegment(objectTip, hitPoint, color, 3.0f, false);
    if (realImage) {
        const Vector2 direction = (imageTip - hitPoint).normalized();
        drawSegment(hitPoint, hitPoint + direction * (1.30f * scene_half_width()), color, 3.0f, false);
    } else {
        const Vector2 direction = (hitPoint - imageTip).normalized();
        drawSegment(hitPoint, hitPoint + direction * (1.30f * scene_half_width()), color, 3.0f, false);
        drawSegment(hitPoint, imageTip, color, 2.0f, true);
    }
}

void drawParallelLensRay(const Vector2& objectTip, const Vector2& imageTip, const OpticalImage& image) {
    const Vector2 hitPoint(lens_diagram_x(), objectTip.y);
    drawRayToImage(objectTip, hitPoint, imageTip, image.real, rayAColor);
}

void drawCenterLensRay(const Vector2& objectTip, const Vector2& imageTip, const OpticalImage& image) {
    const Vector2 center(lens_diagram_x(), 0.0f);
    drawRayToImage(objectTip, center, imageTip, image.real, rayCColor);
}

void drawInfinityLensRays(const Vector2& objectTip) {
    const Vector2 hitPoint(lens_diagram_x(), objectTip.y);
    const Vector2 outgoingDirection(1.0f, -objectTip.y / std::max(lensFocalLength, 0.15f));
    drawSegment(objectTip, hitPoint, rayAColor, 3.0f, false);
    drawSegment(hitPoint, hitPoint + outgoingDirection.normalized() * (0.90f * scene_half_width()), rayAColor, 3.0f, false);

    const Vector2 center(lens_diagram_x(), 0.0f);
    const Vector2 centerDirection = (center - objectTip).normalized();
    drawSegment(objectTip, center, rayCColor, 3.0f, false);
    drawSegment(center, center + centerDirection * (0.90f * scene_half_width()), rayCColor, 3.0f, false);
}

} // namespace

void draw_lens_mode() {
    objectHeight = std::clamp(objectHeight, 0.10f, diagram_object_height_limit());
    objectDistance = clamp_lens_object_distance(objectDistance);
    lensFocalLength = std::clamp(lensFocalLength, 0.18f, 0.45f * scene_half_width());

    const float signedObjectHeight = objectPointsUp ? objectHeight : -objectHeight;
    const float objectBaseX = lens_diagram_x() - objectDistance;
    const OpticalImage image = calculateLensImage(currentLensType, objectDistance, signedObjectHeight, lensFocalLength);

    draw_line(diagram_axis_min_x(), 0.0f, diagram_axis_max_x(), 0.0f, 0.28f, 0.28f, 0.28f, 2.0f);
    drawLensShape();
    drawAxisLabels();

    const float objectColor[3] = { 0.92f, 0.92f, 0.94f };
    drawArrowDiagram(objectBaseX, signedObjectHeight, objectColor, false);

    if (image.valid && !image.atInfinity) {
        const float imageBaseX = lens_diagram_x() + image.distance;
        drawArrowDiagram(imageBaseX, image.height, rayBColor, !image.real);

        if (showConstructionRays) {
            const Vector2 objectTip(objectBaseX, signedObjectHeight);
            const Vector2 imageTip(imageBaseX, image.height);
            drawParallelLensRay(objectTip, imageTip, image);
            drawCenterLensRay(objectTip, imageTip, image);
        }
    } else if (showConstructionRays) {
        drawInfinityLensRays(Vector2(objectBaseX, signedObjectHeight));
    }
}
