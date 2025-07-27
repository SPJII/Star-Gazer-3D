#include "Stars.h"
#include <cmath>
#include <random>
#include <tuple>

// Constructor for Star
Star::Star(float x, float y, float z, float size, float r, float g, float b)
    : x(x), y(y), z(z), size(size), r(r), g(g), b(b)
{
    // Initialize flickerInterval (3-5 seconds) and flickerDuration (0.1-0.5 seconds)
    static std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> distInterval(3.0f, 5.0f);   // Flicker every 3-5 seconds
    std::uniform_real_distribution<float> distDuration(0.1f, 0.5f);   // Flicker lasts 0.1-0.5 seconds
    std::uniform_real_distribution<float> distPhase(0.0f, 1.0f);      // Phase offset

    flickerInterval = distInterval(rng);
    flickerDuration = distDuration(rng);
    flickerPhase = distPhase(rng) * flickerInterval; // Random phase to desynchronize flickers
}

// Render the star with flickering effect
void Star::render(float time) const {
    // Calculate time within the current flicker interval
    float t = fmod(time + flickerPhase, flickerInterval);
    float brightness = 1.0f; // Default brightness

    if (t < flickerDuration) {
        // Star is flickering
        float flickerProgress = t / flickerDuration;
        brightness = 0.8f + 0.2f * sinf(flickerProgress * M_PI); // Smooth flicker between 0.8 and 1.0
    }

    glColor3f(r * brightness, g * brightness, b * brightness);
    glPointSize(size);
    glBegin(GL_POINTS);
    glVertex3f(x, y, z);
    glEnd();
}

// Constructor for Stars
Stars::Stars(int numStars, float radius, GalaxyType type, float posX, float posY, float posZ)
    : posX(posX), posY(posY), posZ(posZ),
    rng(std::random_device{}()),
    distAngle(0.0f, 2.0f * M_PI),
    distSize(1.0f, 3.0f),
    distColor(0.5f, 1.0f)
{
    switch (type) {
    case ELLIPTICAL:
        generateEllipticalGalaxy(numStars, radius);
        break;
    case SPIRAL:
        generateSpiralGalaxy(numStars, radius, false);
        break;
    case BARRED_SPIRAL:
        generateSpiralGalaxy(numStars, radius, true);
        break;
    case LENTICULAR:
        generateLenticularGalaxy(numStars, radius);
        break;
    case IRREGULAR:
        generateIrregularGalaxy(numStars, radius);
        break;
    case DWARF:
        generateDwarfGalaxy(numStars, radius);
        break;
    }
}

// Generate stars for an elliptical galaxy
void Stars::generateEllipticalGalaxy(int numStars, float radius) {
    std::uniform_real_distribution<float> distRadius(0.0f, 1.0f);
    float r_e = radius * 0.5f;  // Effective radius

    for (int i = 0; i < numStars; ++i) {
        float theta = distAngle(rng);
        float phi = acosf(1.0f - 2.0f * distRadius(rng));

        // De Vaucouleurs' profile radius
        float L = distRadius(rng);
        float r = r_e * powf(-logf(L), 4.0f);

        float x = r * sinf(phi) * cosf(theta) + posX;
        float y = r * sinf(phi) * sinf(theta) + posY;
        float z = r * cosf(phi) + posZ;

        float size = distSize(rng);
        float rColor = distColor(rng);
        float gColor = distColor(rng);
        float bColor = distColor(rng);

        stars.emplace_back(x, y, z, size, rColor, gColor, bColor);
    }
}

// Generate stars for a spiral or barred spiral galaxy
void Stars::generateSpiralGalaxy(int numStars, float radius, bool barred) {
    int numArms = 2 + (rng() % 3); // 2 to 4 arms
    float armSeparation = 2.0f * M_PI / numArms;
    float maxRadius = radius;

    std::uniform_real_distribution<float> distRadius(0.0f, maxRadius);
    std::normal_distribution<float> distThickness(0.0f, radius * 0.05f);

    for (int i = 0; i < numStars; ++i) {
        float r = distRadius(rng);
        float armOffset = (rng() % numArms) * armSeparation;
        float theta = r * 0.3f + armOffset + distThickness(rng); // Logarithmic spiral with some randomness
        float z = distThickness(rng) * 2.0f; // Thin disk

        float x = r * cosf(theta) + posX;
        float y = r * sinf(theta) + posY;
        float zPos = z + posZ;

        // Bar component for barred spiral galaxies
        if (barred && (rng() % 10 < 2)) { // 20% stars in the bar
            std::uniform_real_distribution<float> distBar(0.0f, radius * 0.3f);
            float barX = distBar(rng);
            float barY = 0.0f; // Bar along the X-axis
            float barZ = 0.0f;
            x += barX;
            y += barY;
            zPos += barZ;
        }

        float size = distSize(rng);
        float rColor = distColor(rng);
        float gColor = distColor(rng);
        float bColor = distColor(rng);

        stars.emplace_back(x, y, zPos, size, rColor, gColor, bColor);
    }
}

// Generate stars for a lenticular galaxy
void Stars::generateLenticularGalaxy(int numStars, float radius) {
    std::uniform_real_distribution<float> distRadius(0.0f, radius);
    std::normal_distribution<float> distThickness(0.0f, radius * 0.02f);

    for (int i = 0; i < numStars; ++i) {
        float r = distRadius(rng);
        float theta = distAngle(rng);
        float z = distThickness(rng);

        float x = r * cosf(theta) + posX;
        float y = r * sinf(theta) + posY;
        float zPos = z + posZ;

        float size = distSize(rng);
        float rColor = distColor(rng);
        float gColor = distColor(rng);
        float bColor = distColor(rng);

        stars.emplace_back(x, y, zPos, size, rColor, gColor, bColor);
    }
}

// Generate stars for an irregular galaxy
void Stars::generateIrregularGalaxy(int numStars, float radius) {
    // Simple fractal implementation using random clusters
    std::uniform_real_distribution<float> distPosition(-radius, radius);

    int numClusters = std::uniform_int_distribution<int>(5, 15)(rng);

    // Generate cluster centers
    std::vector<std::tuple<float, float, float>> clusters;
    clusters.reserve(numClusters);
    for (int i = 0; i < numClusters; ++i) {
        float cx = distPosition(rng) + posX;
        float cy = distPosition(rng) + posY;
        float cz = distPosition(rng) + posZ;
        clusters.emplace_back(cx, cy, cz);
    }

    // Assign stars to clusters
    std::uniform_int_distribution<int> clusterDist(0, numClusters - 1);
    std::normal_distribution<float> distOffset(0.0f, radius * 0.1f);

    for (int i = 0; i < numStars; ++i) {
        auto [cx, cy, cz] = clusters[clusterDist(rng)];

        float x = cx + distOffset(rng);
        float y = cy + distOffset(rng);
        float z = cz + distOffset(rng);

        float size = distSize(rng);
        float rColor = distColor(rng);
        float gColor = distColor(rng);
        float bColor = distColor(rng);

        stars.emplace_back(x, y, z, size, rColor, gColor, bColor);
    }
}

// Generate stars for a dwarf galaxy
void Stars::generateDwarfGalaxy(int numStars, float radius) {
    std::uniform_real_distribution<float> distPosition(-radius, radius);

    for (int i = 0; i < numStars; ++i) {
        float x = distPosition(rng) + posX;
        float y = distPosition(rng) + posY;
        float z = distPosition(rng) + posZ;

        float size = distSize(rng);
        float rColor = distColor(rng);
        float gColor = distColor(rng);
        float bColor = distColor(rng);

        stars.emplace_back(x, y, z, size, rColor, gColor, bColor);
    }
}

// Render all the stars
void Stars::render(float Time) const {
    for (const auto& star : stars) {
        star.render(Time);
    }
}
