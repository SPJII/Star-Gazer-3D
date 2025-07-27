// CelestialBodies.cpp
// Implementation file for celestial body classes.

#include "CelestialBodies.h"
#include "UI.h"        // For showResourceNodes if needed
#include <GL/glut.h>   // For GLUT functions like renderText
#include <array>
#include <algorithm>   // For std::shuffle, std::max, std::min
#include <cmath>       // For mathematical functions like floor()
#include <iostream>    // For debugging purposes
#include <numeric>     // For std::iota
#include <random>      // For random number generation

/*
    Utility function to render a textured sphere using GLU quadric.
    We tilt the sphere by 90° to align the "poles" with the Y-axis.
*/
void renderSphere(float radius, int slices, int stacks) {
    glPushMatrix();
    glRotatef(90.0f, 1.0f, 0.0f, 0.0f); // Adjust orientation
    GLUquadric* quadric = gluNewQuadric();
    gluQuadricTexture(quadric, GL_TRUE);
    gluSphere(quadric, radius, slices, stacks);
    gluDeleteQuadric(quadric);
    glPopMatrix();
}

/*
    Helper for ray-sphere intersection test.
    rayOrigin: The starting point of the ray (3 floats).
    rayDir:    Normalized direction of the ray (3 floats).
    spherePos: Center of the sphere (3 floats).
    radius:    Sphere radius.
    t:         Output parameter for the 't' where the intersection occurs.
*/
bool raySphereIntersect(const float* rayOrigin,
    const float* rayDir,
    const float* spherePos,
    float radius,
    float& t)
{
    // Compute vector from ray origin to sphere center
    float Lx = spherePos[0] - rayOrigin[0];
    float Ly = spherePos[1] - rayOrigin[1];
    float Lz = spherePos[2] - rayOrigin[2];

    // Project L onto ray direction
    float tca = Lx * rayDir[0] + Ly * rayDir[1] + Lz * rayDir[2];
    if (tca < 0) return false;

    // Compute squared distance from sphere center to the ray
    float d2 = Lx * Lx + Ly * Ly + Lz * Lz - tca * tca;
    if (d2 > radius * radius) return false;

    float thc = sqrtf(radius * radius - d2);
    t = tca - thc;

    return true;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ASTEROID IMPLEMENTATION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Asteroid::Asteroid(float size, float orbitRadius, float orbitSpeed)
    : size(size), orbitRadius(orbitRadius), orbitSpeed(orbitSpeed), orbitAngle(0.0f)
{
    // By default, we can set up some approximate mass for the asteroid
    // and gravitational parameter if needed. Just an example:
    mass = size * 50.0f; // e.g. scaled by some factor
    gravitationalParameter = 6.67430e-11f * mass; // G * M, but often scaled in our scene

    generateShape();

    // Initial position
    positionX = orbitRadius;
    positionY = 0.0f;
    positionZ = 0.0f;

    // orbitInclination can be set externally or remain 0 => near ecliptic
}

void Asteroid::generateShape()
{
    // Generate a rough, jagged shape for the asteroid
    int longitudeBands = 10;
    int latitudeBands = 10;

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> perturbationDist(-0.3f * size, 0.3f * size);

    // Generate vertices
    for (int latNumber = 0; latNumber <= latitudeBands; latNumber++) {
        float theta = latNumber * M_PI / latitudeBands;
        float sinTheta = sin(theta);
        float cosTheta = cos(theta);

        for (int longNumber = 0; longNumber <= longitudeBands; longNumber++) {
            float phi = longNumber * 2.0f * M_PI / longitudeBands;
            float sinPhi = sin(phi);
            float cosPhi = cos(phi);

            float x = cosPhi * sinTheta;
            float y = cosTheta;
            float z = sinPhi * sinTheta;

            // Perturb the vertex position
            float perturbation = perturbationDist(rng);

            Vertex vertex;
            vertex.x = (size + perturbation) * x;
            vertex.y = (size + perturbation) * y;
            vertex.z = (size + perturbation) * z;

            // Normals (approximate)
            vertex.nx = x;
            vertex.ny = y;
            vertex.nz = z;

            // Texture coordinates (not used here)
            vertex.u = 0.0f;
            vertex.v = 0.0f;
            vertex.elevation = 0.0f; // not used for asteroids

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

void Asteroid::render()
{
    glPushMatrix();

    // We apply orbitInclination by rotating around the X-axis (for instance),
    // so that the orbit is slightly out of the ecliptic plane.
    glRotatef(orbitInclination, 1.0f, 0.0f, 0.0f);

    // Position asteroid
    glTranslatef(positionX, positionY, positionZ);

    // For some spin
    glRotatef(orbitAngle * 10.0f, 0.0f, 1.0f, 0.0f);

    // Render asteroid as a rough mesh
    glBegin(GL_TRIANGLES);
    glColor3f(0.5f, 0.5f, 0.5f); // Gray color
    for (size_t i = 0; i < indices.size(); ++i) {
        const Vertex& vertex = vertices[indices[i]];

        glNormal3f(vertex.nx, vertex.ny, vertex.nz);
        glVertex3f(vertex.x, vertex.y, vertex.z);
    }
    glEnd();

    // Restore defaults
    glEnable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);

    glPopMatrix();
}

void Asteroid::update(float movementSpeed)
{
    // Update orbit angle
    orbitAngle += orbitSpeed * movementSpeed;
    if (orbitAngle >= 360.0f) orbitAngle -= 360.0f;

    // Update position in the XZ-plane, then incorporate orbitInclination in render()
    float radOrbitAngle = orbitAngle * (M_PI / 180.0f);
    float cosAngle = cosf(radOrbitAngle);
    float sinAngle = sinf(radOrbitAngle);

    // By default, we keep y=0. For a random small offset or ecliptic variation,
    // we could do: positionY = orbitRadius * sin( tilt ) or so. 
    positionX = orbitRadius * cosAngle;
    positionZ = orbitRadius * sinAngle;
}

bool Asteroid::intersectsRay(const float* rayOrigin, const float* rayDir, float& t)
{
    float spherePos[3] = { positionX, positionY, positionZ };
    return raySphereIntersect(rayOrigin, rayDir, spherePos, size, t);
}

float Asteroid::getPositionX() const { return positionX; }
float Asteroid::getPositionY() const { return positionY; }
float Asteroid::getPositionZ() const { return positionZ; }
float Asteroid::getRadius()   const { return size; }
CelestialBodyType Asteroid::getType() const { return CelestialBodyType::ASTEROID; }

void Asteroid::setOrbitAngle(float angle)
{
    orbitAngle = angle;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// MOON IMPLEMENTATION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Moon::Moon(CelestialBody* parent,
    float distance,
    float size,
    GLuint texture,
    GLuint atmosphereTexture,
    float ascendingNode)
    : parentBody(parent),
    distance(distance),
    size(size),
    textureID(texture),
    atmosphereTextureID(atmosphereTexture),
    orbitAngle(0.0f),
    orbitSpeed(0.5f),
    ascendingNodeAngle(ascendingNode),
    positionX(0.0f),
    positionY(0.0f),
    positionZ(0.0f)
{
    // Example mass or GP:
    mass = size * 500.f;
    gravitationalParameter = 6.67430e-11f * mass;
}

void Moon::render()
{
    if (!parentBody) return;  // safety

    // Compute parent's world position so we can transform relative to it
    float px = parentBody->getPositionX();
    float py = parentBody->getPositionY();
    float pz = parentBody->getPositionZ();

    glPushMatrix();
    // Move to parent's position
    glTranslatef(px, py, pz);

    // 1) Rotate around parent's Y by ascendingNodeAngle -> line of nodes
    glRotatef(ascendingNodeAngle, 0.f, 1.f, 0.f);

    // 2) Tilt orbit by orbitInclination around local X
    glRotatef(orbitInclination, 1.f, 0.f, 0.f);

    // 3) Orbit revolve around Y
    glRotatef(orbitAngle, 0.f, 1.f, 0.f);

    // 4) Translate out to orbit distance
    glTranslatef(distance, 0.f, 0.f);

    // Now the moon is at local (distance,0,0). We can store that final
    // position in world space for picking, etc.:
    // Because we've done glTranslatef(px, py, pz) + these rotates, let's
    // get the final modelview matrix and extract the translation to store
    // in positionX/Y/Z. Another approach is to do the trig ourselves in update().
    // We'll do the trig approach in update() so that picking always has the correct position:

    // Render the moon surface
    glPushMatrix();
    glColor3f(1.0f, 1.0f, 1.0f);
    glBindTexture(GL_TEXTURE_2D, textureID);
    renderSphere(size, 30, 30);
    glPopMatrix();

    // Render atmosphere
    glPushMatrix();
    glRotatef(orbitAngle * 0.5f, 0.f, 1.f, 0.f);
    glBindTexture(GL_TEXTURE_2D, atmosphereTextureID);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);
    glColor4f(1.f, 1.f, 1.f, 0.5f);
    renderSphere(size + 0.05f, 30, 30);
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    glPopMatrix();

    glPopMatrix();
}

void Moon::update(float movementSpeed)
{
    if (!parentBody) return;

    // Advance orbit angle
    orbitAngle += orbitSpeed * movementSpeed;
    if (orbitAngle >= 360.f) orbitAngle -= 360.f;

    // We do the actual trig to track world position for picking collisions, etc.
    // 1) Get parent's position
    float px = parentBody->getPositionX();
    float py = parentBody->getPositionY();
    float pz = parentBody->getPositionZ();

    // 2) Ascending node
    //   We'll treat ascendingNodeAngle as rotation around parent's Y axis
    float radNode = ascendingNodeAngle * (M_PI / 180.f);
    float cosNode = cosf(radNode);
    float sinNode = sinf(radNode);

    // 3) OrbitInclination is around X axis
    float radInc = orbitInclination * (M_PI / 180.f);
    float cosInc = cosf(radInc);
    float sinInc = sinf(radInc);

    // 4) orbitAngle is revolve around Y
    float radOrbit = orbitAngle * (M_PI / 180.f);
    float cosO = cosf(radOrbit);
    float sinO = sinf(radOrbit);

    // The net transformation from local (distance,0,0):
    //  localX = distance
    //  localY = 0
    //  localZ = 0
    // Then rotate by orbitAngle about Y, rotate by inclination about X, and rotate by ascendingNode about Y again.
    // We'll do the multiplications carefully (there are multiple ways, but here's a straightforward approach).

    // Step A: revolve around local Y by orbitAngle
    float rx = distance * cosO;  // revolve in XZ plane
    float ry = 0.f;
    float rz = distance * sinO;

    // Step B: tilt around X by orbitInclination
    //  y' = ry*cosInc - rz*sinInc
    //  z' = ry*sinInc + rz*cosInc
    float ry2 = rx * 0.f + ry * cosInc - rz * sinInc;
    float rz2 = rx * 0.f + ry * sinInc + rz * cosInc;
    float rx2 = rx; // x unaffected by rotation around X

    // Step C: rotate around parent's Y by ascendingNodeAngle
    //  x'' = rx2*cosNode - rz2*sinNode
    //  z'' = rx2*sinNode + rz2*cosNode
    float rx3 = rx2 * cosNode - rz2 * sinNode;
    float rz3 = rx2 * sinNode + rz2 * cosNode;
    float ry3 = ry2; // unchanged

    // So final world position = parent's position + (rx3, ry3, rz3)
    positionX = px + rx3;
    positionY = py + ry3;
    positionZ = pz + rz3;
}

bool Moon::intersectsRay(const float* rayOrigin, const float* rayDir, float& t)
{
    float spherePos[3] = { positionX, positionY, positionZ };
    return raySphereIntersect(rayOrigin, rayDir, spherePos, size, t);
}

float Moon::getPositionX() const { return positionX; }
float Moon::getPositionY() const { return positionY; }
float Moon::getPositionZ() const { return positionZ; }
float Moon::getRadius()   const { return size; }
CelestialBodyType Moon::getType() const { return CelestialBodyType::MOON; }

void Moon::setOrbitSpeed(float speed)
{
    orbitSpeed = speed;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ASTEROIDBELT IMPLEMENTATION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AsteroidBelt::AsteroidBelt(float innerRadius, float outerRadius, int asteroidCount)
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> radiusDist(innerRadius, outerRadius);
    std::uniform_real_distribution<float> orbitSpeedDist(0.01f, 0.05f);
    std::uniform_real_distribution<float> sizeDist(0.1f, 0.5f);
    std::uniform_real_distribution<float> angleDist(0.0f, 360.0f);

    for (int i = 0; i < asteroidCount; ++i) {
        float r = radiusDist(rng);
        float sp = orbitSpeedDist(rng);
        float sz = sizeDist(rng);

        Asteroid* asteroid = new Asteroid(sz, r, sp);
        // random orbit angle
        asteroid->setOrbitAngle(angleDist(rng));

        // Optionally randomize orbit inclination
        float inc = (float)(rand() % 10 - 5); // e.g. -5..+5 degrees
        asteroid->orbitInclination = inc;

        asteroids.push_back(asteroid);
    }
}

AsteroidBelt::~AsteroidBelt()
{
    for (Asteroid* asteroid : asteroids) {
        delete asteroid;
    }
    asteroids.clear();
}

void AsteroidBelt::render()
{
    for (Asteroid* asteroid : asteroids) {
        asteroid->render();
    }
}

void AsteroidBelt::update(float movementSpeed)
{
    for (Asteroid* asteroid : asteroids) {
        asteroid->update(movementSpeed);
    }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// RINGSYSTEM IMPLEMENTATION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RingSystem::RingSystem(float innerRadius, float outerRadius, int asteroidCount)
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> radiusDist(innerRadius, outerRadius);
    std::uniform_real_distribution<float> orbitSpeedDist(0.01f, 0.05f);
    std::uniform_real_distribution<float> sizeDist(0.1f, 0.5f);
    std::uniform_real_distribution<float> angleDist(0.0f, 360.0f);

    for (int i = 0; i < asteroidCount; ++i) {
        float r = radiusDist(rng);
        float sp = orbitSpeedDist(rng);
        float sz = sizeDist(rng);

        Asteroid* asteroid = new Asteroid(sz, r, sp);
        asteroid->setOrbitAngle(angleDist(rng));

        // Typically, rings might have minimal inclination relative
        // to the planet's equatorial plane. You could set that here as well.
        asteroid->orbitInclination = 0.0f; // or planet’s tilt

        asteroids.push_back(asteroid);
    }
}

RingSystem::~RingSystem()
{
    for (Asteroid* asteroid : asteroids) {
        delete asteroid;
    }
    asteroids.clear();
}

void RingSystem::render()
{
    for (Asteroid* asteroid : asteroids) {
        asteroid->render();
    }
}

void RingSystem::update(float movementSpeed)
{
    for (Asteroid* asteroid : asteroids) {
        asteroid->update(movementSpeed);
    }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// PLANET IMPLEMENTATION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Planet::Planet(float radius,
    float atmosphereRadius,
    GLuint texture,
    GLuint atmosphereTexture,
    Moon* moon,
    float orbitRadius,
    float orbitSpeed)
    : radius(radius), atmosphereRadius(atmosphereRadius),
    textureID(texture), atmosphereTextureID(atmosphereTexture),
    rotationY(0.0f), moon(moon), orbitRadius(orbitRadius),
    orbitAngle(0.0f), orbitSpeed(orbitSpeed)
{
    // If we got a single moon, store it in the vector
    if (moon) {
        moons.push_back(moon);
    }

    // Some arbitrary density approach: e.g. mass ~ volume * density
    // We'll just do something simple:
    float volume = (4.0f / 3.0f) * (float)M_PI * radius * radius * radius;
    float density = 5000.0f; // e.g. Earth-like or something
    mass = volume * density;
    gravitationalParameter = 6.67430e-11f * mass;

    // Initial position
    positionX = orbitRadius;
    positionY = 0.0f;
    positionZ = 0.0f;

    // By default, planet is near ecliptic plane => orbitInclination ~ 0
    // Axial tilt can also be set externally (like Earth ~ 23.5f)
    axialTilt = 0.0f;
}

Planet::~Planet()
{
    // If we allocated a single "moon" ourselves, we can delete them here
    // but watch out if we didn't allocate the pointer. 
    // This example doesn't forcibly delete them, but you might choose to:
    // for (auto* m : moons) {
    //     delete m;
    // }
    // moons.clear();
}

void Planet::render()
{
    glPushMatrix();

    // Apply orbit inclination
    glRotatef(orbitInclination, 1.0f, 0.0f, 0.0f);

    // Move out from star
    glTranslatef(positionX, positionY, positionZ);

    // Now tilt the planet by its axialTilt, around the X-axis or Y-axis
    // For instance, if we want a typical Earth-like tilt around the X-axis:
    glRotatef(axialTilt, 0.0f, 0.0f, 1.0f); // tilt axis for "spin"
    glRotatef(rotationY, 0.0f, 1.0f, 0.0f); // spin

    // Render planet surface
    glColor3f(1.0f, 1.0f, 1.0f);
    glBindTexture(GL_TEXTURE_2D, textureID);
    renderSphere(radius, 50, 50);

    // Render atmosphere
    glPushMatrix();
    // optionally rotate atmosphere more slowly
    glBindTexture(GL_TEXTURE_2D, atmosphereTextureID);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);

    glColor4f(1.0f, 1.0f, 1.0f, 0.3f);
    renderSphere(atmosphereRadius, 50, 50);

    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    glPopMatrix(); // atmosphere

    // If planet has moons, render them
    for (Moon* m : moons) {
        m->render();
    }

    glPopMatrix();
}

void Planet::update(float movementSpeed)
{
    // Orbit
    orbitAngle += orbitSpeed * movementSpeed;
    if (orbitAngle >= 360.0f) orbitAngle -= 360.0f;

    // Convert to radians
    float radOrbitAngle = orbitAngle * (M_PI / 180.0f);
    positionX = orbitRadius * cosf(radOrbitAngle);
    positionZ = orbitRadius * sinf(radOrbitAngle);
    positionY = 0.0f; // ecliptic plane, offset by orbitInclination in render

    // Rotate planet on its axis
    // We can scale rotation by e.g. 10 * movementSpeed
    rotationY += 10.0f * movementSpeed;
    if (rotationY >= 360.0f) rotationY -= 360.0f;

    // Update moons
    for (Moon* m : moons) {
        m->update(movementSpeed);
    }
}

bool Planet::intersectsRay(const float* rayOrigin, const float* rayDir, float& t)
{
    // We do a simple sphere test with planet's center
    float spherePos[3] = { positionX, positionY, positionZ };
    return raySphereIntersect(rayOrigin, rayDir, spherePos, radius, t);
}

float Planet::getPositionX() const { return positionX; }
float Planet::getPositionY() const { return positionY; }
float Planet::getPositionZ() const { return positionZ; }
float Planet::getRotationY() const { return rotationY; }
float Planet::getRadius()   const { return radius; }
CelestialBodyType Planet::getType() const { return CelestialBodyType::PLANET; }

void Planet::setRotation(float rotY)
{
    rotationY = rotY;
}

void Planet::setZoom(float z)
{
    // Not used in this sample, but could be used if we want to 'zoom' planet
}

Moon* Planet::getMoon() const
{
    return moon;
}

void Planet::addMoon(Moon* m)
{
    if (m) {
        moons.push_back(m);
    }
}

const std::vector<Moon*>& Planet::getMoons() const
{
    return moons;
}

// Sun implementation

// Constructor
Sun::Sun(SunType type, float orbitRadius, float orbitSpeed, GLuint texture)
    : type(type), orbitRadius(orbitRadius), orbitSpeed(orbitSpeed),
    textureID(texture), rotation(0.0f), noiseOffsetX(0.0f), noiseOffsetY(0.0f), noiseOffsetZ(0.0f),
    noiseSpeed(0.01f), distribution(0.0f, 1.0f)
{
    // Seed the random number generator
    rng.seed(static_cast<unsigned int>(time(0)));

    // Initialize Perlin noise
    initializePerlin();

    // Assign properties based on spectral type
    switch (type) {
    case SunType::O:
        radius = 10.0f * 109.0f / 1.0f; // Scale by 109 relative to planetRadius=1.0f
        temperature = 30000.0f;
        colorR = 0.5f; colorG = 0.5f; colorB = 1.0f; // Blue-White
        break;
    case SunType::B:
        radius = 6.6f * 109.0f / 1.0f;
        temperature = 20000.0f;
        colorR = 0.6f; colorG = 0.4f; colorB = 1.0f; // Blue-White
        break;
    case SunType::A:
        radius = 5.0f * 109.0f / 1.0f;
        temperature = 9000.0f;
        colorR = 0.8f; colorG = 0.8f; colorB = 0.8f; // White
        break;
    case SunType::F:
        radius = 1.6f * 109.0f / 1.0f;
        temperature = 7000.0f;
        colorR = 1.0f; colorG = 1.0f; colorB = 0.8f; // Yellow-White
        break;
    case SunType::G:
        radius = 1.0f * 109.0f / 1.0f; // Earth radius = 1.0f, Sun = 109.0f
        temperature = 5800.0f;
        colorR = 1.0f; colorG = 1.0f; colorB = 0.6f; // Yellow
        break;
    case SunType::K:
        radius = 0.8f * 109.0f / 1.0f;
        temperature = 4500.0f;
        colorR = 1.0f; colorG = 0.7f; colorB = 0.3f; // Orange
        break;
    case SunType::M:
        radius = 0.6f * 109.0f / 1.0f;
        temperature = 3500.0f;
        colorR = 1.0f; colorG = 0.3f; colorB = 0.0f; // Red
        break;
    default:
        radius = 1.0f * 109.0f / 1.0f; // Default to G-type
        temperature = 5800.0f;
        colorR = 1.0f; colorG = 1.0f; colorB = 0.6f; // Yellow
        break;
    }

    // Generate the mesh for the sun
    generateMesh();
}

// Destructor
Sun::~Sun() {
    // Clean up resources if necessary
}

// Initialize Perlin noise permutation vector
void Sun::initializePerlin() {
    p.resize(256);
    // Initialize permutation vector with values 0 to 255
    for (int i = 0; i < 256; ++i)
        p[i] = i;

    // Shuffle using the random number generator
    std::shuffle(p.begin(), p.end(), rng);

    // Duplicate the permutation vector
    p.insert(p.end(), p.begin(), p.end());
}

// Perlin noise function implementation
float Sun::perlinNoise(float x, float y, float z) {
    // Find unit cube that contains point
    int X = static_cast<int>(floor(x)) & 255;
    int Y = static_cast<int>(floor(y)) & 255;
    int Z = static_cast<int>(floor(z)) & 255;

    // Find relative x, y, z of point in cube
    x -= floor(x);
    y -= floor(y);
    z -= floor(z);

    // Compute fade curves for each of x, y, z
    float u = fade(x);
    float v = fade(y);
    float w = fade(z);

    // Hash coordinates of the 8 cube corners
    int A = p[X] + Y;
    int AA = p[A] + Z;
    int AB = p[A + 1] + Z;
    int B = p[X + 1] + Y;
    int BA = p[B] + Z;
    int BB = p[B + 1] + Z;

    // Add blended results from 8 corners of the cube
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

// Fade function as defined by Ken Perlin
float Sun::fade(float t) {
    return t * t * t * (t * (t * 6 - 15) + 10);
}

// Linear interpolation
float Sun::lerp(float a, float b, float t) {
    return a + t * (b - a);
}

// Gradient function as defined by Ken Perlin
float Sun::grad(int hash, float x, float y, float z) {
    int h = hash & 15;          // Convert low 4 bits of hash code
    float u = h < 8 ? x : y;    // Into 12 gradient directions
    float v = h < 4 ? y : (h == 12 || h == 14 ? x : z);
    return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
}

// Generates the sun's mesh with procedural perturbations
void Sun::generateMesh() {
    // Mesh resolution
    int longitudeBands = 128; // Increased for higher resolution
    int latitudeBands = 64;

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

            // Apply Perlin noise for surface perturbation
            float noiseValue = perlinNoise(x * 5.0f + noiseOffsetX, y * 5.0f + noiseOffsetY, z * 5.0f + noiseOffsetZ);
            float elevation = noiseValue * 0.05f; // Adjust amplitude as needed

            SunVertex vertex;
            vertex.x = (radius + elevation) * x;
            vertex.y = (radius + elevation) * y;
            vertex.z = (radius + elevation) * z;

            // Normalize for normals
            float length = sqrt(vertex.x * vertex.x + vertex.y * vertex.y + vertex.z * vertex.z);
            vertex.nx = vertex.x / length;
            vertex.ny = vertex.y / length;
            vertex.nz = vertex.z / length;

            // Texture coordinates
            vertex.u = static_cast<float>(longNumber) / longitudeBands;
            vertex.v = static_cast<float>(latNumber) / latitudeBands;

            vertex.elevation = elevation;

            vertices.push_back(vertex);
        }
    }

    // Generate indices for triangle strips
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

// Animates the sun's surface by updating noise offsets
void Sun::animateSurface() {
    noiseOffsetX += noiseSpeed;
    noiseOffsetY += noiseSpeed * 0.5f;
    noiseOffsetZ += noiseSpeed * 0.2f;

    // Regenerate mesh with new noise offsets
    vertices.clear();
    indices.clear();
    generateMesh();
}

// Render the sun
void Sun::render() {
    glPushMatrix();
    glTranslatef(getPositionX(), getPositionY(), getPositionZ());
    glRotatef(rotation, 0.0f, 1.0f, 0.0f); // Rotate for animation

    // Enable lighting
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    // Bind texture if available
    if (textureID != 0) {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, textureID);
    }
    else {
        glDisable(GL_TEXTURE_2D);
    }

    // Set material properties with emissive color for glow
    GLfloat materialAmbient[] = { colorR * 0.2f, colorG * 0.2f, colorB * 0.2f, 1.0f };
    GLfloat materialDiffuse[] = { colorR, colorG, colorB, 1.0f };
    GLfloat materialSpecular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat materialEmissive[] = { colorR * 0.5f, colorG * 0.5f, colorB * 0.5f, 1.0f }; // Increased emissive for glow
    GLfloat materialShininess[] = { 50.0f };

    glMaterialfv(GL_FRONT, GL_AMBIENT, materialAmbient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, materialDiffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, materialSpecular);
    glMaterialfv(GL_FRONT, GL_EMISSION, materialEmissive); // Set emissive color
    glMaterialfv(GL_FRONT, GL_SHININESS, materialShininess);

    // Render the mesh
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < indices.size(); i += 3) {
        SunVertex& v1 = vertices[indices[i]];
        SunVertex& v2 = vertices[indices[i + 1]];
        SunVertex& v3 = vertices[indices[i + 2]];

        // Normals
        glNormal3f(v1.nx, v1.ny, v1.nz);
        glNormal3f(v2.nx, v2.ny, v2.nz);
        glNormal3f(v3.nx, v3.ny, v3.nz);

        // Texture coordinates
        if (textureID != 0) {
            glTexCoord2f(v1.u, v1.v);
            glTexCoord2f(v2.u, v2.v);
            glTexCoord2f(v3.u, v3.v);
        }

        // Vertex positions
        glVertex3f(v1.x, v1.y, v1.z);
        glVertex3f(v2.x, v2.y, v2.z);
        glVertex3f(v3.x, v3.y, v3.z);
    }
    glEnd();

    // Disable texture mapping if it was enabled
    if (textureID != 0) {
        glBindTexture(GL_TEXTURE_2D, 0);
        glDisable(GL_TEXTURE_2D);
    }

    // Reset emissive color
    GLfloat noEmissive[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    glMaterialfv(GL_FRONT, GL_EMISSION, noEmissive);

    // Restore OpenGL states to defaults
    glEnable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);

    glPopMatrix(); // End of planet rendering
}

// Update the sun's state (animation)
void Sun::update(float movementSpeed) {
    // Update rotation for animation
    rotation += movementSpeed * 0.1f;
    if (rotation >= 360.0f)
        rotation -= 360.0f;

    // Animate the surface
    animateSurface();
}

// Intersection test (ray-sphere)
bool Sun::intersectsRay(const float* rayOrigin, const float* rayDir, float& t) {
    float spherePos[3] = { getPositionX(), getPositionY(), getPositionZ() };
    float radiusSquared = radius * radius;
    // Ray-sphere intersection formula
    float Lx = spherePos[0] - rayOrigin[0];
    float Ly = spherePos[1] - rayOrigin[1];
    float Lz = spherePos[2] - rayOrigin[2];
    float tca = Lx * rayDir[0] + Ly * rayDir[1] + Lz * rayDir[2];
    if (tca < 0) return false;
    float d2 = Lx * Lx + Ly * Ly + Lz * Lz - tca * tca;
    if (d2 > radiusSquared) return false;
    float thc = sqrt(radiusSquared - d2);
    t = tca - thc;
    return true;
}

// Getters for position
float Sun::getPositionX() const { return 0.0f; } // Sun is at the origin
float Sun::getPositionY() const { return 0.0f; }
float Sun::getPositionZ() const { return 0.0f; }

float Sun::getRadius() const
{
    return radius;
}

// Get celestial body type
CelestialBodyType Sun::getType() const { return CelestialBodyType::SUN; }