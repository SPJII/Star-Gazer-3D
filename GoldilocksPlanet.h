#pragma once

#include "CelestialBodies.h"
#include "resource.h"
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <cmath>

class GoldilocksPlanet : public Planet
{
public:
    // -- Keep your existing BiomeType and ProvinceType enums --
    enum class BiomeType {
        DeepOcean, Ocean, River,
        TropicalRainforest, TropicalDryForest, Savanna, Desert, TropicalMontaneForest,
        SubtropicalForest, Grassland, Steppe, SubtropicalMontaneForest,
        TemperateForest, TemperateMontaneForest,
        BorealForest, Taiga, Tundra, Ice,

        // --- BEACH BIOME ---
        Beach
    };

    enum class ProvinceType {
        SeaCoastal, LandCoastal, Island, Inland, Mountain, Ocean
    };

    // Structure representing a vertex on the planet's mesh
    struct Vertex {
        float x, y, z;        // Position
        float nx, ny, nz;     // Normal
        float u, v;           // Texture coords
        float elevation;      // Elevation
        bool  isRiver;        // Is part of a river
        bool  isOcean;        // Is part of the ocean
        float temperature;
        float moisture;
        float latitude;
        float longitude;
        int   provinceId;
        float oceanGyre;
        BiomeType biome;
        float weatherSeverity;
    };

    // Structure representing a province
    struct Province {
        int id;
        ProvinceType type;
        std::vector<int> vertexIndices;
        float colorR, colorG, colorB, colorA;
    };

    // Structure representing a resource node
    struct ResourceNode {
        float positionX, positionY, positionZ;
        ResourceType type;
    };

    // Axial tilt in degrees
    float axialTilt;

    // Current orbital position (angle around the sun)
    float orbitalPosition;

    // Seasons
    enum class Season { Spring, Summer, Autumn, Winter };
    Season currentSeason;

    // Constructor
    GoldilocksPlanet(
        float r,
        float atmosphereR,
        GLuint texture,
        GLuint atmosphereTexture,
        Moon* m,
        float orbitR,
        float orbitS
    );

    // Override from Planet
    void update(float movementSpeed) override;
    void render() override;

private:

    // ~~~~~~~~~~~~~ LOD Structures & Storage ~~~~~~~~~~~~~
    struct LODLevel {
        std::vector<Vertex>        vertices;
        std::vector<unsigned int>  indices;
    };
    // We’ll keep multiple LOD meshes. Example: 3 levels
    std::vector<LODLevel> lodLevels;

    // The highest LOD data is the one on which we do the advanced
    // generation (erosion, climate, etc.). We keep them separate to
    // avoid changing too much of the existing code that references
    // "vertices"/"indices".
    std::vector<Vertex>        highVertices;
    std::vector<unsigned int>  highIndices;

    // Keep your existing provinces & resource nodes
    std::vector<Province>      provinces;
    std::vector<ResourceNode>  resourceNodes;

    // Perlin noise
    std::vector<int> p;
    std::mt19937 rng;
    std::uniform_real_distribution<float> distribution;
    unsigned int noiseSeed;

    // ~~~~~~~~~~~~~~ LOD Helpers ~~~~~~~~~~~~~~
    // Generate one LOD's sphere mesh (without advanced processes)
    void generateLODTerrain(int resolution, LODLevel& lod);

    // Generate the set of LOD levels (e.g., low, medium)
    void generateAllLODs();

    // Decide which LOD to use based on camera distance
    int  chooseLOD(float distanceToCamera);

    // Basic method to render a given LOD level
    void renderLODLevel(const LODLevel& lod);

    // ~~~~~~~~~~~~~~ Original Terrain & Planet Gen ~~~~~~~~~~~~~~
    // (We will store the “big” result in highVertices/highIndices)
    void generateHighResTerrain();        // The original “generateTerrain()” logic
    void applyErosion();                 // Erosion on highVertices
    void calculateClimate();             // Climate on highVertices
    void generateRivers();               // Rivers on highVertices

    // Province generation (works on highVertices)
    void generateProvinces();
    void classifyProvinceTypes();
    void assignProvinceColors();
    void renderProvinces();

    // Resource generation
    void generateResourceNodes();

    // Season logic
    void updateSeasons(float movementSpeed);
    void calculateSeasonalEffects();

    // Perlin noise + helpers
    float getElevation(float x, float y, float z);
    float perlinNoise(float x, float y, float z);
    float fade(float t);
    float lerp(float a, float b, float t);
    float grad(int hash, float x, float y, float z);

    // Biome determination
    BiomeType determineBiome(Vertex& vertex);
};
