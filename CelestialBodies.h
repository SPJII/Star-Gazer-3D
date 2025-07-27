// CelestialBodies.h
// Header file for celestial body classes.

#ifndef CELESTIALBODIES_H
#define CELESTIALBODIES_H

#include <GL/glut.h>
#include <vector>
#include <string>
#include <algorithm>
#include <random>

// Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Forward declarations
class Moon;
class RingSystem;

/*
    We keep vertices for rendering (e.g., in asteroids, planets).
    This struct is used by some derived classes (like Asteroid, Planet).
*/
struct Vertex {
    float x, y, z;        // Position coordinates
    float nx, ny, nz;     // Normal vector
    float u, v;           // Texture coordinates
    float elevation;      // Additional data (e.g., for noise, terrain)
};

/*
    Defines the type of the celestial body for easy identification.
*/
enum class CelestialBodyType {
    SUN,
    PLANET,
    MOON,
    ASTEROID
};

/*
    Base class for all celestial bodies. Includes new fields and methods
    to handle gravitational influences and ecliptic-plane modeling.
*/
class CelestialBody {
public:
    virtual ~CelestialBody() = default;

    // Render the celestial body (pure virtual - must be overridden)
    virtual void render() = 0;

    // Update the celestial body's state (e.g., orbit, rotation)
    // movementSpeed is a global speed factor.
    virtual void update(float movementSpeed) = 0;

    // Returns true if a ray intersects this body, for picking
    virtual bool intersectsRay(const float* rayOrigin, const float* rayDir, float& t) = 0;

    // Accessors for position and radius
    virtual float getPositionX() const = 0;
    virtual float getPositionY() const = 0;
    virtual float getPositionZ() const = 0;
    virtual float getRadius()   const = 0;

    // Type of this body
    virtual CelestialBodyType getType() const = 0;

    // New functions for gravitational calculations:
    //    * getMass() - returns an approximate mass
    //    * getGravitationalParameter() - GM
    // These defaults can be overridden in derived classes if needed.
    virtual float getMass() const { return mass; }
    virtual float getGravitationalParameter() const { return gravitationalParameter; }

    // Basic gravitational data
    float mass = 1.0f;                 // Derived classes can compute from volume & density
    float gravitationalParameter = 1.0f; // G*M (some constant for the body)

    // By default, we incorporate a small ecliptic inclination approach
    // or rotation tilt. Derived classes can override or set it properly.
    float orbitInclination = 0.0f;     // In degrees, relative to the ecliptic plane
    float axialTilt = 0.0f;           // Planet or moon's axial tilt

    /*
       We can do more advanced general relativity handling
       (perihelion precession, etc.) but here we only store
       placeholders. If you want to track e.g.
       argumentOfPerihelion, ascendingNode, etc.,
       you can add them here.
    */
};

// -------------------------------------------------------
// Asteroid class
// -------------------------------------------------------
class Asteroid : public CelestialBody {
public:
    Asteroid(float size, float orbitRadius, float orbitSpeed);

    void render() override;
    void update(float movementSpeed) override;
    bool intersectsRay(const float* rayOrigin, const float* rayDir, float& t) override;

    float getPositionX() const override;
    float getPositionY() const override;
    float getPositionZ() const override;
    float getRadius()    const override;
    CelestialBodyType getType() const override;

    // For external code to set the asteroid's initial orbit angle
    void setOrbitAngle(float angle);

private:
    float size;            // Physical radius of the asteroid
    float orbitRadius;     // Distance from some central point (e.g. star or planet)
    float orbitSpeed;      // Orbital speed factor
    float orbitAngle;      // Current angle in degrees
    float positionX, positionY, positionZ;

    // For rendering the asteroid
    void generateShape();
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
};

// -------------------------------------------------------
// Moon class (satellite of a planet)
// -------------------------------------------------------
class Moon : public CelestialBody {
public:
    // Pass in the parent body (planet, sun, etc.), plus ascending node angle.
    Moon(CelestialBody* parent,
        float distance,
        float size,
        GLuint texture,
        GLuint atmosphereTexture,
        float ascendingNodeAngle = 0.0f);

    void render() override;
    void update(float movementSpeed) override;
    bool intersectsRay(const float* rayOrigin, const float* rayDir, float& t) override;

    float getPositionX() const override;
    float getPositionY() const override;
    float getPositionZ() const override;
    float getRadius()   const override;
    CelestialBodyType getType() const override;

    void setOrbitSpeed(float speed);

    // Adjust orbit inclination or ascending node externally if desired
    void setOrbitInclination(float degrees) { orbitInclination = degrees; }
    void setAscendingNodeAngle(float degrees) { ascendingNodeAngle = degrees; }

private:
    // Link to the body we orbit (planet or sun)
    CelestialBody* parentBody;

    float distance;               // Orbital distance from parent
    float size;                   // Physical radius of the moon
    GLuint textureID;             // Surface texture
    GLuint atmosphereTextureID;   // Atmosphere texture

    float orbitAngle;             // Current orbit angle around parent
    float orbitSpeed;             // Speed of orbit
    float ascendingNodeAngle;     // Rotation around parent's Y-axis setting line of nodes

    float positionX, positionY, positionZ; // Current position in world space
};

// -------------------------------------------------------
// AsteroidBelt class
// -------------------------------------------------------
class AsteroidBelt {
public:
    // Construct an asteroid belt from an inner radius to outer radius
    // containing 'asteroidCount' asteroids
    AsteroidBelt(float innerRadius, float outerRadius, int asteroidCount);
    ~AsteroidBelt();

    void render();
    void update(float movementSpeed);

private:
    std::vector<Asteroid*> asteroids;
};

// -------------------------------------------------------
// RingSystem class
// -------------------------------------------------------
class RingSystem {
public:
    // Similar to an asteroid belt but around a planet
    RingSystem(float innerRadius, float outerRadius, int asteroidCount);
    ~RingSystem();

    void render();
    void update(float movementSpeed);

private:
    std::vector<Asteroid*> asteroids;
};

// -------------------------------------------------------
// Planet class
// -------------------------------------------------------
class Planet : public CelestialBody {
public:
    /*
       radius            = planet radius
       atmosphereRadius  = radius of the atmosphere
       texture           = planet surface texture
       atmosphereTexture = planet atmosphere texture
       moon              = optional primary moon pointer
       orbitRadius       = distance from the star (in realistic scale)
       orbitSpeed        = orbital speed factor
    */
    Planet(float radius,
        float atmosphereRadius,
        GLuint texture,
        GLuint atmosphereTexture,
        Moon* moon,
        float orbitRadius,
        float orbitSpeed);

    virtual ~Planet();

    void render() override;
    void update(float movementSpeed) override;
    bool intersectsRay(const float* rayOrigin, const float* rayDir, float& t) override;

    float getPositionX() const override;
    float getPositionY() const override;
    float getPositionZ() const override;
    float getRotationY() const;  // Planet's spin around its axis
    float getRadius()    const override;
    CelestialBodyType getType()  const override;

    void setRotation(float rotY);
    void setZoom(float z); // Not used in the current sample

    Moon* getMoon() const;
    void  addMoon(Moon* moon); // Add more moons
    const std::vector<Moon*>& getMoons() const;

    // Additional method for planet's axial tilt
    void setAxialTilt(float tiltDegrees) { axialTilt = tiltDegrees; }

protected:
    float radius;            // Radius of the planet
    float atmosphereRadius;  // Radius of its atmosphere
    GLuint textureID;        // Planet surface texture
    GLuint atmosphereTextureID; // Atmosphere texture
    float rotationY;         // Rotation around Y-axis

    // The planet can have multiple moons
    Moon* moon;                   // The "primary" moon (if provided)
    std::vector<Moon*> moons;     // All moons

    float orbitRadius;    // Distance from star
    float orbitAngle;     // Planet's current orbit angle
    float orbitSpeed;     // Speed of orbit
    float positionX, positionY, positionZ;

    // We can use axialTilt from base CelestialBody to tilt the planet's axis
};

// -------------------------------------------------------
// Enum for different spectral classes of the sun
// -------------------------------------------------------
enum class SunType {
    O, B, A, F, G, K, M
};

/*
    For the sun's own specialized mesh, we define a separate vertex
    struct to incorporate any animations specific to the star's surface.
*/
struct SunVertex {
    float x, y, z;   // Position coordinates
    float nx, ny, nz; // Normal vector
    float u, v;      // Texture coordinates
    float elevation; // For simulated "solar noise" or flares
};

// -------------------------------------------------------
// Sun class representing the central star
// -------------------------------------------------------
// Sun class representing the central star
class Sun : public CelestialBody {
public:
    Sun(SunType type,
        float orbitRadius,
        float orbitSpeed,
        GLuint texture = 0);
    ~Sun();

    void render() override;
    void update(float movementSpeed) override;
    bool intersectsRay(const float* rayOrigin, const float* rayDir, float& t) override;

    float getPositionX() const override;
    float getPositionY() const override;
    float getPositionZ() const override;
    float getRadius()   const override;
    CelestialBodyType getType() const override;

protected:
    // Added member variables to track the sun's position
    float positionX, positionY, positionZ;

private:
    SunType type;
    float radius;
    float temperature;
    float colorR, colorG, colorB;
    float rotation;
    float orbitRadius;
    float orbitSpeed;
    GLuint textureID;
    std::vector<SunVertex> vertices;
    std::vector<unsigned int> indices;
    float noiseOffsetX, noiseOffsetY, noiseOffsetZ;
    float noiseSpeed;
    std::vector<int> p;
    std::mt19937 rng;
    std::uniform_real_distribution<float> distribution;

    // Helper methods
    void initializePerlin();
    float perlinNoise(float x, float y, float z);
    float fade(float t);
    float lerp(float a, float b, float t);
    float grad(int hash, float x, float y, float z);

    void generateMesh();
    void animateSurface();
};

// -------------------------------------------------------
// Utility function to render a sphere
// -------------------------------------------------------
void renderSphere(float radius, int slices, int stacks);

#endif // CELESTIALBODIES_H
