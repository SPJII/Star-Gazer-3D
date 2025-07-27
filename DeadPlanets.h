// DeadPlanets.h

#ifndef DEADPLANETS_H
#define DEADPLANETS_H

#include <GL/glut.h>
#include <vector>
#include <random>
#include <cmath>
#include "CelestialBodies.h"
#include "GoldilocksPlanet.h"



class TerrestrialPlanet : public Planet {
public:
    TerrestrialPlanet(float radius, GLuint texture, GLuint atmosphereTexture, float orbitRadius, float orbitSpeed);

    void generateTerrain();
    float getElevation(float x, float y, float z);
    float perlinNoise(float x, float y, float z);

    float fade(float t);

    float lerp(float a, float b, float t);

    float grad(int hash, float x, float y, float z);

    void render() override;

private:
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
    std::vector<int> p; // Permutation vector for Perlin noise

    void setTerrainColor(float elevation);

    void update(float movementSpeed);

    std::mt19937 rng;
    std::uniform_real_distribution<float> distribution;
    unsigned int noiseSeed;
};

class GasGiantPlanet : public Planet {
public:
    GasGiantPlanet(float radius, GLuint texture, float orbitRadius, float orbitSpeed);
    ~GasGiantPlanet();

    void render() override;       // Renders the planet

    void generateMoons();         // Generates moons based on planet size
    void generateRings();         // Generates rings if applicable
    void update(float movementSpeed) override;

private:
    std::vector<Vertex> vertices;               // Vertices for the planet's mesh
    std::vector<unsigned int> indices;          // Indices for rendering


    void generateCloudPatterns();               // Generates the planet's cloud patterns
    void initializePerlin();                    // Initializes the permutation vector for Perlin noise
    float perlinNoise(float x, float y, float z); // Perlin noise function
    float fade(float t);                        // Fade function for Perlin noise
    float lerp(float a, float b, float t);      // Linear interpolation
    float grad(int hash, float x, float y, float z); // Gradient function for Perlin noise

    std::vector<int> p;                         // Permutation vector for Perlin noise

    // Add a member to store the ring system
    RingSystem* ringSystem;
};
#endif // !DEADPLANETS_H

