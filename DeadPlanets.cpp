// DeadPlanets.cpp

#include "DeadPlanets.h"
#include <cmath>
#include <vector>
#include <random>
#include <numeric> // For std::iota
#include <map>
#include <algorithm> // For std::min
#include <GL/glut.h> // For renderText or other GLUT functionalities

// TerrestrialPlanet implementation

// TerrestrialPlanet implementation

TerrestrialPlanet::TerrestrialPlanet(float radius, GLuint texture, GLuint atmosphereTexture, float orbitRadius, float orbitSpeed)
    : Planet(radius, radius * 1.1f, texture, atmosphereTexture, nullptr, orbitRadius, orbitSpeed),
    distribution(0.0f, 1.0f) {
    // Seed the random number generator.
    noiseSeed = static_cast<unsigned int>(std::time(0));
    rng.seed(noiseSeed);

    // Initialize permutation vector for Perlin noise.
    p.resize(256);
    for (int i = 0; i < 256; ++i) {
        p[i] = i;
    }

    // Shuffle the permutation vector using std::shuffle and rng.
    std::shuffle(p.begin(), p.end(), rng);

    // Duplicate the permutation vector.
    p.insert(p.end(), p.begin(), p.end());

    generateTerrain(); // Generate the planet's terrain.
}

void TerrestrialPlanet::generateTerrain() {
    // Increased mesh resolution for higher detail.
    int longitudeBands = 64; // Increased from 256 to 512
    int latitudeBands = 32;  // Increased from 128 to 256

    // Generate vertices.
    for (int latNumber = 0; latNumber <= latitudeBands; latNumber++) {
        float theta = latNumber * M_PI / latitudeBands; // Polar angle.
        float sinTheta = sin(theta);
        float cosTheta = cos(theta);

        for (int longNumber = 0; longNumber <= longitudeBands; longNumber++) {
            float phi = longNumber * 2 * M_PI / longitudeBands; // Azimuthal angle.
            float sinPhi = sin(phi);
            float cosPhi = cos(phi);

            // Unit sphere coordinates.
            float x = cosPhi * sinTheta;
            float y = cosTheta;
            float z = sinPhi * sinTheta;

            // Get elevation from noise function.
            float elevation = getElevation(x, y, z);

            // Vertex position adjusted by elevation.
            Vertex vertex;
            vertex.x = x * (getRadius() + elevation);
            vertex.y = y * (getRadius() + elevation);
            vertex.z = z * (getRadius() + elevation);

            // Normal vector for lighting.
            vertex.nx = x;
            vertex.ny = y;
            vertex.nz = z;

            // Texture coordinates.
            vertex.u = static_cast<float>(longNumber) / longitudeBands;
            vertex.v = static_cast<float>(latNumber) / latitudeBands;

            // Initialize other properties.
            vertex.elevation = elevation;

            vertices.push_back(vertex);
        }
    }

    // Generate indices for triangle drawing, ensuring seamless mesh.
    for (int latNumber = 0; latNumber < latitudeBands; latNumber++) {
        for (int longNumber = 0; longNumber < longitudeBands; longNumber++) {
            int first = (latNumber * (longitudeBands + 1)) + longNumber;
            int second = first + longitudeBands + 1;

            indices.push_back(first);
            indices.push_back(second);
            indices.push_back(first + 1);

            indices.push_back(second);
            indices.push_back(second + 1);
            indices.push_back(first + 1);
        }
    }
}

// Get elevation based on Perlin noise.
float TerrestrialPlanet::getElevation(float x, float y, float z) {
    // Fractal Perlin noise parameters.
    float elevation = 0.0f;
    float frequency = 0.6f;
    float amplitude = 0.06f;
    float persistence = 0.8f;
    float lacunarity = 5.0f;
    int octaves = 12;

    // Generate fractal Perlin noise for elevation.
    for (int i = 0; i < octaves; ++i) {
        float noiseValue = perlinNoise(x * frequency, y * frequency, z * frequency);
        elevation += noiseValue * amplitude;

        amplitude *= persistence;
        frequency *= lacunarity;
    }

    // Scale elevation based on planet radius.
    elevation *= getRadius() * 0.05f; // Adjust as needed.

    return elevation;
}

// Perlin noise function implementation.
float TerrestrialPlanet::perlinNoise(float x, float y, float z) {
    // Find unit cube that contains point.
    int X = static_cast<int>(floor(x)) & 255;
    int Y = static_cast<int>(floor(y)) & 255;
    int Z = static_cast<int>(floor(z)) & 255;

    // Find relative x, y, z of point in cube.
    x -= floor(x);
    y -= floor(y);
    z -= floor(z);

    // Compute fade curves for each of x, y, z.
    float u = fade(x);
    float v = fade(y);
    float w = fade(z);

    // Hash coordinates of the 8 cube corners.
    int A = p[X] + Y;
    int AA = p[A] + Z;
    int AB = p[A + 1] + Z;
    int B = p[X + 1] + Y;
    int BA = p[B] + Z;
    int BB = p[B + 1] + Z;

    // Add blended results from 8 corners of the cube.
    float res = lerp(w,
        lerp(v,
            lerp(u, grad(p[AA], x, y, z),
                grad(p[BA], x - 1, y, z)),
            lerp(u, grad(p[AB], x, y - 1, z),
                grad(p[BB], x - 1, y - 1, z))
        ),
        lerp(v,
            lerp(u, grad(p[AA + 1], x, y, z - 1),
                grad(p[BA + 1], x - 1, y, z - 1)),
            lerp(u, grad(p[AB + 1], x, y - 1, z - 1),
                grad(p[BB + 1], x - 1, y - 1, z - 1))
        )
    );

    return res;
}

// Fade function as defined by Ken Perlin.
float TerrestrialPlanet::fade(float t) {
    return t * t * t * (t * (t * 6 - 15) + 10);
}

// Linear interpolation function.
float TerrestrialPlanet::lerp(float a, float b, float t) {
    return a + t * (b - a);
}

// Gradient function as defined by Ken Perlin.
float TerrestrialPlanet::grad(int hash, float x, float y, float z) {
    int h = hash & 15;          // Convert low 4 bits of hash code.
    float u = h < 8 ? x : y;    // Into 12 gradient directions.
    float v = h < 4 ? y : (h == 12 || h == 14 ? x : z);
    return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
}

// Render the terrestrial planet.
void TerrestrialPlanet::render() {
    glPushMatrix();
    glTranslatef(getPositionX(), 0.0f, getPositionZ());
    glRotatef(getRotationY(), 0.0f, 1.0f, 0.0f);

    // Enable lighting
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    // Disable texture mapping to use dynamic color palette.
    glDisable(GL_TEXTURE_2D);

    // Begin rendering triangles
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < indices.size(); ++i) {
        Vertex& vertex = vertices[indices[i]];

        // Set color based on elevation with smooth gradients
        setTerrainColor(vertex.elevation);

        // Set normal for lighting
        glNormal3f(vertex.nx, vertex.ny, vertex.nz);

        // Set vertex position
        glVertex3f(vertex.x, vertex.y, vertex.z);
    }
    glEnd();

    // Re-enable texture mapping if needed
    glEnable(GL_TEXTURE_2D);

    // Render the atmosphere
    glPushMatrix();
    glRotatef(getRotationY() + 5.0f, 0.0f, 1.0f, 0.0f); // Atmosphere rotation

    // Bind the atmosphere texture
    glBindTexture(GL_TEXTURE_2D, atmosphereTextureID);

    // Enable blending and manage depth testing
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE); // Disable depth writing

    // Set translucency and color
    glColor4f(1.0f, 1.0f, 1.0f, 0.3f); // Slightly transparent atmosphere

    // Render the atmosphere sphere with higher resolution
    renderSphere(atmosphereRadius, 512, 256); // Higher resolution sphere

    // Restore depth mask and blending settings
    glDepthMask(GL_TRUE); // Re-enable depth writing
    glDisable(GL_BLEND);
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f); // Reset opacity

    glPopMatrix(); // End of atmosphere rendering

    // Render moons if any
    for (auto& moon : moons) {
        moon->render();
    }

    // Restore OpenGL states to defaults
    glEnable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);

    glPopMatrix(); // End of planet rendering
}

// Set terrain color based on elevation with smooth gradients
void TerrestrialPlanet::setTerrainColor(float elevation) {
    // Define elevation thresholds and corresponding colors
    struct ColorPoint {
        float elevationThreshold;
        float r, g, b;
    };

    std::vector<ColorPoint> colorPoints = {
        {0.00f * getRadius(), 0.3f, 0.3f, 0.3f},   // Deep Rock - Dark Gray
        {0.02f * getRadius(), 0.5f, 0.5f, 0.5f},   // Rocky Plains - Light Gray
        {0.04f * getRadius(), 0.6f, 0.3f, 0.1f},   // Brownish Mountains - Brown
        {0.06f * getRadius(), 0.8f, 0.8f, 0.8f},   // Snow-Capped Peaks - White
        {0.08f * getRadius(), 0.9f, 0.5f, 0.2f},   // Volcanic Areas - Orange
        {0.10f * getRadius(), 0.8f, 0.0f, 0.0f}    // Extreme Peaks - Red
    };

    // Find the two color points the elevation falls between
    ColorPoint lower = colorPoints[0];
    ColorPoint upper = colorPoints.back();

    for (size_t i = 0; i < colorPoints.size() - 1; ++i) {
        if (elevation >= colorPoints[i].elevationThreshold && elevation < colorPoints[i + 1].elevationThreshold) {
            lower = colorPoints[i];
            upper = colorPoints[i + 1];
            break;
        }
    }

    // Calculate interpolation factor
    float t = 0.0f;
    if (upper.elevationThreshold - lower.elevationThreshold > 0.0f) {
        t = (elevation - lower.elevationThreshold) / (upper.elevationThreshold - lower.elevationThreshold);
    }

    // Apply a non-linear interpolation for more natural gradients
    t = t * t * (3 - 2 * t); // Smoothstep function

    // Interpolate between lower and upper colors
    float r = lower.r + t * (upper.r - lower.r);
    float g = lower.g + t * (upper.g - lower.g);
    float b = lower.b + t * (upper.b - lower.b);

    glColor3f(r, g, b);
}

void TerrestrialPlanet::update(float movementSpeed) {
    Planet::update(movementSpeed);
    // Additional updates if necessary
}

// Implement GasGiantPlanet class methods

GasGiantPlanet::GasGiantPlanet(float radius, GLuint texture, float orbitRadius, float orbitSpeed)
    : Planet(radius, radius * 1.0f, texture, 0, nullptr, orbitRadius, orbitSpeed), ringSystem(nullptr) {
    initializePerlin();
    generateCloudPatterns();
}

GasGiantPlanet::~GasGiantPlanet() {
    // Delete moons
    for (auto& moon : moons) {
        delete moon;
    }
    moons.clear();

    // Delete ring system
    if (ringSystem) {
        delete ringSystem;
        ringSystem = nullptr;
    }
}

void GasGiantPlanet::initializePerlin() {
    p.resize(256);
    std::iota(p.begin(), p.end(), 0);
    std::shuffle(p.begin(), p.end(), std::mt19937{ std::random_device{}() });
    p.insert(p.end(), p.begin(), p.end());
}

float GasGiantPlanet::perlinNoise(float x, float y, float z) {
    // Same implementation as in TerrestrialPlanet
    int X = static_cast<int>(floor(x)) & 255;
    int Y = static_cast<int>(floor(y)) & 255;
    int Z = static_cast<int>(floor(z)) & 255;

    x -= floor(x);
    y -= floor(y);
    z -= floor(z);

    float u = fade(x);
    float v = fade(y);
    float w = fade(z);

    int A = p[X] + Y;
    int AA = p[A] + Z;
    int AB = p[A + 1] + Z;
    int B = p[X + 1] + Y;
    int BA = p[B] + Z;
    int BB = p[B + 1] + Z;

    float res = lerp(w,
        lerp(v,
            lerp(u, grad(p[AA], x, y, z),
                grad(p[BA], x - 1, y, z)),
            lerp(u, grad(p[AB], x, y - 1, z),
                grad(p[BB], x - 1, y - 1, z))
        ),
        lerp(v,
            lerp(u, grad(p[AA + 1], x, y, z - 1),
                grad(p[BA + 1], x - 1, y, z - 1)),
            lerp(u, grad(p[AB + 1], x, y - 1, z - 1),
                grad(p[BB + 1], x - 1, y - 1, z - 1))
        )
    );
    return res;
}

float GasGiantPlanet::fade(float t) {
    return t * t * t * (t * (t * 6 - 15) + 10);
}

float GasGiantPlanet::lerp(float a, float b, float t) {
    return a + t * (b - a);
}

float GasGiantPlanet::grad(int hash, float x, float y, float z) {
    int h = hash & 15;
    float u = h < 8 ? x : y;
    float v = h < 4 ? y : (h == 12 || h == 14 ? x : z);
    return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
}

void GasGiantPlanet::generateCloudPatterns() {
    // Generate the planet's cloud patterns

    int longitudeBands = 128;
    int latitudeBands = 64;

    float noiseScale = 5.0f;

    // Generate vertices
    for (int latNumber = 0; latNumber <= latitudeBands; latNumber++) {
        float theta = latNumber * M_PI / latitudeBands;
        float sinTheta = sin(theta);
        float cosTheta = cos(theta);

        for (int longNumber = 0; longNumber <= longitudeBands; longNumber++) {
            float phi = longNumber * 2 * M_PI / longitudeBands;
            float sinPhi = sin(phi);
            float cosPhi = cos(phi);

            float x = cosPhi * sinTheta;
            float y = cosTheta;
            float z = sinPhi * sinTheta;

            // Use noise to generate cloud patterns
            float noiseValue = perlinNoise(x * noiseScale, y * noiseScale, z * noiseScale);

            Vertex vertex;
            vertex.x = radius * x;
            vertex.y = radius * y;
            vertex.z = radius * z;

            // Normals
            vertex.nx = x;
            vertex.ny = y;
            vertex.nz = z;

            // Texture coordinates
            vertex.u = 1.0f - ((float)longNumber / longitudeBands);
            vertex.v = 1.0f - ((float)latNumber / latitudeBands);

            // Store noise value for coloring
            vertex.elevation = noiseValue;

            vertices.push_back(vertex);
        }
    }

    // Generate indices
    for (int latNumber = 0; latNumber < latitudeBands; latNumber++) {
        for (int longNumber = 0; longNumber < longitudeBands; longNumber++) {
            int first = (latNumber * (longitudeBands + 1)) + longNumber;
            int second = first + longitudeBands + 1;

            indices.push_back(first);
            indices.push_back(second);
            indices.push_back(first + 1);

            indices.push_back(second);
            indices.push_back(second + 1);
            indices.push_back(first + 1);
        }
    }
}

void GasGiantPlanet::generateMoons() {
    // Generate moons based on planet size
    std::mt19937 rng(std::random_device{}());
    int minMoons = static_cast<int>(10 * (radius / 4.0f)); // Minimum 10 moons for smallest gas giant
    int maxMoons = static_cast<int>(90 * (radius / 11.0f)); // Up to 90 moons for largest gas giant

    if (minMoons < 10) minMoons = 10;
    if (maxMoons > 90) maxMoons = 90;

    std::uniform_int_distribution<int> moonCountDist(minMoons, maxMoons);
    std::uniform_real_distribution<float> moonSizeDist(radius * 0.01f, radius * 0.1f);
    std::uniform_real_distribution<float> orbitDistDist(radius * 1.5f, radius * 5.0f);
    std::uniform_real_distribution<float> orbitSpeedDist(0.1f, 0.5f);

    int moonCount = moonCountDist(rng);
    float lastOrbitDistance = radius;

    for (int i = 0; i < moonCount; ++i) {
        float moonSize = moonSizeDist(rng);
        float orbitDistance = lastOrbitDistance + moonSize * 5.0f; // Ensure separation
        float orbitSpeed = orbitSpeedDist(rng);

        Moon* moon = new Moon(this,orbitDistance, moonSize, 0, 0); // Texture IDs set to 0 for simplicity
        moon->setOrbitSpeed(orbitSpeed);

        addMoon(moon);
        lastOrbitDistance = orbitDistance;
    }
}

void GasGiantPlanet::generateRings() {
    // Generate rings if applicable (e.g., for Saturn)
    if (radius >= 9.0f) { // Assume large gas giants have rings
        ringSystem = new RingSystem(
            radius * 1.1f,    // Inner radius just outside the planet
            radius * 2.0f,    // Outer radius
            100              // Number of asteroids in the ring
        );
    }
}

void GasGiantPlanet::render() {
    glPushMatrix();
    glTranslatef(getPositionX(), 0.0f, getPositionZ());
    glRotatef(getRotationY(), 0.0f, 1.0f, 0.0f);

    // Bind the planet's texture if any
    glBindTexture(GL_TEXTURE_2D, textureID);

    // Render the planet's surface using the generated mesh
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < indices.size(); ++i) {
        Vertex& vertex = vertices[indices[i]];

        // Use the noise value to affect the color (simulate cloud bands)
        float noiseValue = vertex.elevation;
        float color = 0.5f + 0.5f * noiseValue; // Map noise value to [0,1]

        glColor3f(color, color * 0.8f, color * 0.5f); // Adjust colors accordingly

        glNormal3f(vertex.nx, vertex.ny, vertex.nz);
        glTexCoord2f(vertex.u, vertex.v);
        glVertex3f(vertex.x, vertex.y, vertex.z);
    }
    glEnd();

    // Render the ring system if it exists
    if (ringSystem) {
        glPushMatrix();
        ringSystem->render();
        glPopMatrix();
    }

    // Render moons
    for (auto& moon : moons) {
        moon->render();
    }

    // Restore OpenGL states to defaults
    glEnable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);

    glPopMatrix(); // End of planet rendering
}

void GasGiantPlanet::update(float movementSpeed) {
    // Update rotation
    rotationY += 0.1f * movementSpeed;
    if (rotationY >= 360.0f) rotationY -= 360.0f;

    // Update orbit angle
    orbitAngle += orbitSpeed * movementSpeed;
    if (orbitAngle >= 360.0f) orbitAngle -= 360.0f;

    // Update position
    positionX = orbitRadius * cosf(orbitAngle * M_PI / 180.0f);
    positionZ = orbitRadius * sinf(orbitAngle * M_PI / 180.0f);

    // Update moons
    for (auto& moon : moons) {
        moon->update(movementSpeed);
    }

    // Update ring system if it exists
    if (ringSystem) {
        ringSystem->update(movementSpeed);
    }
}