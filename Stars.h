// Stars.h
#ifndef STARS_H
#define STARS_H

#include <GL/glut.h>
#include <vector>
#include <random>
#include <cmath>
#include "CelestialBodies.h"

class Star {
public:
    float x, y, z;          // Position in 3D space
    float size;             // Size of the star
    float r, g, b;          // Base RGB color values
    float flickerPhase;     // Phase offset for flickering
    float flickerInterval;  // Time between flickers
    float flickerDuration;  // Duration of the flicker

    Star(float x, float y, float z, float size, float r, float g, float b);
    void render(float Time) const;
};

class Stars {
public:
    enum GalaxyType { ELLIPTICAL, SPIRAL, BARRED_SPIRAL, LENTICULAR, IRREGULAR, DWARF };

    Stars(int numStars, float radius, GalaxyType type, float posX, float posY, float posZ);

    void render(float time) const;

private:
    std::vector<Star> stars;      // Collection of stars
    float posX, posY, posZ;       // Position of the galaxy
    std::mt19937 rng;             // Random number generator
    std::uniform_real_distribution<float> distAngle;
    std::uniform_real_distribution<float> distSize;
    std::uniform_real_distribution<float> distColor;

    void generateEllipticalGalaxy(int numStars, float radius);
    void generateSpiralGalaxy(int numStars, float radius, bool barred);
    void generateLenticularGalaxy(int numStars, float radius);
    void generateIrregularGalaxy(int numStars, float radius);
    void generateDwarfGalaxy(int numStars, float radius);
};

#endif // STARS_H