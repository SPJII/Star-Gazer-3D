#include "GoldilocksPlanet.h"
#include "CelestialBodies.h"
#include "resource.h"
#include "UI.h"       // For showResourceNodes, showProvinces
#include <GL/glut.h>  // For GLUT functions (rendering text, etc.)
#include <algorithm>
#include <cmath>
#include <queue>
#include <cfloat>
#include <ctime>
#include <iostream>
#include <map>

// Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// External camera variables (assuming from UI.cpp)
extern float camX, camY, camZ;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//              Constructor & Initial Setup
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GoldilocksPlanet::GoldilocksPlanet(
    float r,
    float atmosphereR,
    GLuint texture,
    GLuint atmosphereTexture,
    Moon* m,
    float orbitR,
    float orbitS
)
    : Planet(r, atmosphereR, texture, atmosphereTexture, m, orbitR, orbitS),
    distribution(0.0f, 1.0f)
{
    // Randomize axial tilt between 0 and 45 degrees
    axialTilt = distribution(rng) * 45.0f;
    orbitalPosition = 0.0f;
    currentSeason = Season::Spring;

    // Seed the random number generator
    noiseSeed = static_cast<unsigned int>(std::time(nullptr));
    rng.seed(noiseSeed);

    // Prepare Perlin permutation vector
    p.resize(256);
    for (int i = 0; i < 256; ++i) {
        p[i] = i;
    }
    std::shuffle(p.begin(), p.end(), rng);
    p.insert(p.end(), p.begin(), p.end());

    // 1) Generate the highest LOD terrain: 1000×1000
    generateHighResTerrain(); // Creates highVertices/highIndices
    applyErosion();
    calculateClimate();
    generateRivers();

    // 2) Generate provinces, resource nodes, etc., on the high-res data
    generateProvinces();
    classifyProvinceTypes();
    assignProvinceColors();
    generateResourceNodes();

    // 3) Create additional lower LOD levels
    generateAllLODs(); // e.g. 128×128, 512×512, etc.
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                       LOD Generation
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void GoldilocksPlanet::generateAllLODs()
{
    // Example: We create 2 additional LOD levels
    lodLevels.resize(2);

    // LOD 0 → 128×128
    generateLODTerrain(128, lodLevels[0]);

    // LOD 1 → 256×256
    generateLODTerrain(256, lodLevels[1]);

    // We already have the highest LOD (512×512) in highVertices/highIndices
    // We'll treat that as "LOD 2" in chooseLOD().
}

void GoldilocksPlanet::generateLODTerrain(int resolution, LODLevel& lod)
{
    lod.vertices.clear();
    lod.indices.clear();

    // We do a simpler approach—just create a sphere & get elevation (no erosion, etc.)
    int longitudeBands = resolution;
    int latitudeBands = resolution;

    std::vector<float> elevations;
    elevations.reserve((latitudeBands + 1) * (longitudeBands + 1));

    // Generate vertices
    for (int lat = 0; lat <= latitudeBands; ++lat) {
        float theta = lat * M_PI / latitudeBands;
        float sinTheta = sin(theta);
        float cosTheta = cos(theta);

        for (int lon = 0; lon <= longitudeBands; ++lon) {
            float phi = lon * 2.f * M_PI / longitudeBands;
            float sinPhi = sin(phi);
            float cosPhi = cos(phi);

            float x = cosPhi * sinTheta;
            float y = cosTheta;
            float z = sinPhi * sinTheta;

            float elevation = getElevation(x, y, z); // Perlin-based
            elevations.push_back(elevation);

            GoldilocksPlanet::Vertex v;
            float finalR = getRadius() + elevation;
            v.x = x * finalR;
            v.y = y * finalR;
            v.z = z * finalR;
            v.nx = x;
            v.ny = y;
            v.nz = z;
            v.u = (float)lon / longitudeBands;
            v.v = (float)lat / latitudeBands;
            v.elevation = elevation;
            v.isRiver = false;
            v.isOcean = false; // We’ll do a simple sea-level pass
            v.latitude = asin(v.ny);
            v.longitude = atan2(v.nz, v.nx);
            v.provinceId = -1;
            v.oceanGyre = 0.f;
            v.biome = BiomeType::Ocean; // default
            v.weatherSeverity = 0.f;
            lod.vertices.push_back(v);
        }
    }

    // Compute sea level quickly for a ~15-40% land coverage
    // same logic as your existing approach
    int desiredLandPct = static_cast<int>(distribution(rng) * 25) + 15; // 15%-40%
    std::vector<float> sortedElev = elevations;
    std::sort(sortedElev.begin(), sortedElev.end());
    int landIndex = (int)((1.0f - (desiredLandPct / 100.f)) * sortedElev.size());
    float seaLevel = sortedElev[landIndex];

    // Adjust all to set isOcean
    for (size_t i = 0; i < lod.vertices.size(); i++) {
        auto& v = lod.vertices[i];
        v.elevation -= seaLevel;
        float radiusWithElev = getRadius() + v.elevation;
        v.x = v.nx * radiusWithElev;
        v.y = v.ny * radiusWithElev;
        v.z = v.nz * radiusWithElev;
        v.isOcean = (v.elevation < 0.f);
        // We skip advanced climate, rivers, etc.
    }

    // Build indices
    for (int lat = 0; lat < latitudeBands; ++lat) {
        for (int lon = 0; lon < longitudeBands; ++lon) {
            int first = (lat * (longitudeBands + 1)) + lon;
            int second = first + longitudeBands + 1;
            int firstPlus1 = (lon + 1) % (longitudeBands + 1) + lat * (longitudeBands + 1);
            int secondPlus1 = firstPlus1 + longitudeBands + 1;

            lod.indices.push_back(first);
            lod.indices.push_back(second);
            lod.indices.push_back(firstPlus1);

            lod.indices.push_back(second);
            lod.indices.push_back(secondPlus1);
            lod.indices.push_back(firstPlus1);
        }
    }
}

int GoldilocksPlanet::chooseLOD(float distanceToCamera)
{
    // Example thresholds:
    if (distanceToCamera > 1800.0f) {
        return 0; // LOD 0 (lowest)
    }
    else if (distanceToCamera > 600.0f) {
        return 1; // LOD 1 (medium)
    }
    else {
        // LOD 2 is actually the "highVertices/highIndices"
        return 2;
    }
}

void GoldilocksPlanet::renderLODLevel(const LODLevel& lod)
{
    glBegin(GL_TRIANGLES);
    for (auto idx : lod.indices) {
        const Vertex& v = lod.vertices[idx];
        // If you want a simpler color approach for low LOD
        if (v.isOcean)
            glColor3f(0.0f, 0.0f, 0.5f);
        else
            glColor3f(0.4f, 0.8f, 0.4f);

        glNormal3f(v.nx, v.ny, v.nz);
        glVertex3f(v.x, v.y, v.z);
    }
    glEnd();
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//       Full-Resolution Terrain & Erosion/Climate Logic
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void GoldilocksPlanet::generateHighResTerrain()
{
    int longitudeBands = 512;
    int latitudeBands = 512;

    // We store into highVertices/highIndices
    highVertices.clear();
    highIndices.clear();

    std::vector<float> elevations;
    elevations.reserve((latitudeBands + 1) * (longitudeBands + 1));

    for (int latNumber = 0; latNumber <= latitudeBands; latNumber++) {
        float theta = latNumber * M_PI / latitudeBands;
        float sinTheta = sin(theta);
        float cosTheta = cos(theta);

        for (int longNumber = 0; longNumber <= longitudeBands; longNumber++) {
            float phi = longNumber * 2.f * M_PI / longitudeBands;
            float sinPhi = sin(phi);
            float cosPhi = cos(phi);

            float x = cosPhi * sinTheta;
            float y = cosTheta;
            float z = sinPhi * sinTheta;

            float elevation = getElevation(x, y, z);
            elevations.push_back(elevation);

            Vertex v;
            v.x = x * (getRadius() + elevation);
            v.y = y * (getRadius() + elevation);
            v.z = z * (getRadius() + elevation);
            v.nx = x;
            v.ny = y;
            v.nz = z;
            v.u = (float)longNumber / longitudeBands;
            v.v = (float)latNumber / latitudeBands;

            v.elevation = elevation;
            v.isRiver = false;
            v.latitude = asin(v.ny);
            v.longitude = atan2(v.nz, v.nx);
            v.provinceId = -1;
            v.oceanGyre = 0.f;
            v.biome = BiomeType::Ocean;
            v.weatherSeverity = 0.f;

            highVertices.push_back(v);
        }
    }

    // Adjust sea level
    int desiredLandPercentage = static_cast<int>(distribution(rng) * 25) + 15; // 15-40%
    std::vector<float> sortedElevations = elevations;
    std::sort(sortedElevations.begin(), sortedElevations.end());
    int landIndex = (int)((1.0f - desiredLandPercentage / 100.f) * sortedElevations.size());
    float seaLevel = sortedElevations[landIndex];

    for (size_t i = 0; i < highVertices.size(); ++i) {
        highVertices[i].elevation -= seaLevel;
        float h = getRadius() + highVertices[i].elevation;
        highVertices[i].x = highVertices[i].nx * h;
        highVertices[i].y = highVertices[i].ny * h;
        highVertices[i].z = highVertices[i].nz * h;
        highVertices[i].isOcean = (highVertices[i].elevation < 0.f);
    }

    // Build indices
    for (int latNumber = 0; latNumber < latitudeBands; latNumber++) {
        for (int longNumber = 0; longNumber < longitudeBands; longNumber++) {
            int first = (latNumber * (longitudeBands + 1)) + longNumber;
            int second = first + longitudeBands + 1;

            int firstPlusOne = (longNumber + 1) % (longitudeBands + 1) + latNumber * (longitudeBands + 1);
            int secondPlusOne = firstPlusOne + longitudeBands + 1;

            highIndices.push_back(first);
            highIndices.push_back(second);
            highIndices.push_back(firstPlusOne);

            highIndices.push_back(second);
            highIndices.push_back(secondPlusOne);
            highIndices.push_back(firstPlusOne);
        }
    }
}

void GoldilocksPlanet::applyErosion()
{
    int iterations = 15;
    int vertexCount = (int)highVertices.size();
    int longitudeBands = 512;
    int latitudeBands = 512;

    // Precompute neighbors
    std::vector<std::vector<int>> neighborIndices(vertexCount);
    for (int i = 0; i < vertexCount; i++) {
        int lat = i / (longitudeBands + 1);
        int lon = i % (longitudeBands + 1);

        std::vector<int> neighbors;
        int offsets[4][2] = { {-1,0},{0,1},{1,0},{0,-1} };
        for (auto& o : offsets) {
            int nLat = lat + o[0];
            int nLon = lon + o[1];

            if (nLon < 0)           nLon = longitudeBands;
            else if (nLon > longitudeBands) nLon = 0;

            if (nLat >= 0 && nLat <= latitudeBands) {
                int idx = nLat * (longitudeBands + 1) + nLon;
                neighbors.push_back(idx);
            }
        }
        neighborIndices[i] = neighbors;
    }

    for (int iter = 0; iter < iterations; ++iter) {
        std::vector<float> deltaElevation(vertexCount, 0.f);
        for (int i = 0; i < vertexCount; ++i) {
            for (auto nIdx : neighborIndices[i]) {
                float delta = highVertices[i].elevation - highVertices[nIdx].elevation;
                float talusAngle = 0.0003f;

                if (delta > talusAngle) {
                    float sediment = delta / 4.f;
                    deltaElevation[i] -= sediment;
                    deltaElevation[nIdx] += sediment;
                }
            }
        }
        // apply changes
        for (int i = 0; i < vertexCount; i++) {
            highVertices[i].elevation += deltaElevation[i];
            float h = getRadius() + highVertices[i].elevation;
            highVertices[i].x = highVertices[i].nx * h;
            highVertices[i].y = highVertices[i].ny * h;
            highVertices[i].z = highVertices[i].nz * h;
            highVertices[i].isOcean = (highVertices[i].elevation < 0.f);
        }
    }
}

void GoldilocksPlanet::calculateClimate()
{
    int vertexCount = (int)highVertices.size();
    int longitudeBands = 512;
    int latitudeBands = 512;

    // init
    for (auto& v : highVertices) {
        v.temperature = 0.f;
        v.moisture = 0.f;
        v.oceanGyre = 0.f;
    }

    // ocean gyres
    for (auto& v : highVertices) {
        if (v.isOcean) {
            // northern hemisphere
            if (v.latitude > 0) {
                if (v.longitude >= 0 && v.longitude < M_PI)
                    v.oceanGyre = 1.f;
                else
                    v.oceanGyre = -1.f;
            }
            else {
                if (v.longitude >= 0 && v.longitude < M_PI)
                    v.oceanGyre = -1.f;
                else
                    v.oceanGyre = 1.f;
            }
        }
    }

    // temperature
    for (auto& v : highVertices) {
        float latFactor = cos(v.latitude);
        float elevFactor = exp(-v.elevation / (getRadius() * 0.005f));
        v.temperature = latFactor * elevFactor + 0.1f;
        if (v.temperature > 1.f) v.temperature = 1.f;
    }

    // moisture transport
    for (int step = 0; step < 3; step++) {
        for (int i = 0; i < vertexCount; i++) {
            auto& v = highVertices[i];
            if (v.isOcean) {
                v.moisture += v.temperature * 0.5f;
                continue;
            }
            int lat = i / (longitudeBands + 1);
            int lon = i % (longitudeBands + 1);

            int dLat = (v.latitude > 0) ? -1 : +1;
            int nLat = lat + dLat;
            int nLon = lon;

            if (nLon < 0)         nLon = longitudeBands;
            else if (nLon > longitudeBands) nLon = 0;

            if (nLat < 0) nLat = 0;
            if (nLat > latitudeBands) nLat = latitudeBands;

            int nIdx = nLat * (longitudeBands + 1) + nLon;
            auto& nbr = highVertices[nIdx];

            float moistureTransfer = v.moisture;
            if (nbr.elevation > v.elevation) {
                float eDiff = nbr.elevation - v.elevation;
                float mLoss = eDiff * 5.f;
                moistureTransfer -= mLoss;
            }
            if (moistureTransfer < 0.f) moistureTransfer = 0.f;

            nbr.moisture += (moistureTransfer * 0.9f);
            v.moisture *= 0.1f;
        }
    }
    // ocean current proximity
    for (auto& v : highVertices) {
        if (!v.isOcean) {
            int lat = (int)((v.latitude + M_PI / 2) / (M_PI)*latitudeBands);
            int lon = (int)((v.longitude + M_PI) / (2.f * M_PI) * longitudeBands);

            int dLon = (v.oceanGyre > 0) ? 1 : -1;
            int nLon = lon + dLon;
            if (nLon < 0) nLon = longitudeBands;
            else if (nLon > longitudeBands) nLon = 0;

            int neighborIdx = lat * (longitudeBands + 1) + nLon;
            if (neighborIdx >= 0 && neighborIdx < vertexCount) {
                auto& nbr = highVertices[neighborIdx];
                if (nbr.isOcean) {
                    v.moisture += nbr.temperature * 0.2f;
                }
            }
        }
        v.moisture -= 0.05f;
        if (v.moisture < 0.f) v.moisture = 0.f;
    }

    // clamp
    for (auto& v : highVertices) {
        if (v.temperature < 0.f) v.temperature = 0.f;
        if (v.temperature > 1.f) v.temperature = 1.f;
        if (v.moisture < 0.f)    v.moisture = 0.f;
        if (v.moisture > 1.f)    v.moisture = 1.f;
    }

    // determine biome
    for (auto& v : highVertices) {
        v.biome = determineBiome(v);
    }
}

void GoldilocksPlanet::generateRivers()
{
    int vertexCount = (int)highVertices.size();
    std::vector<float> flow(vertexCount, 0.f);
    int longitudeBands = 512;
    int latitudeBands = 512;

    // neighbors
    std::vector<std::vector<int>> neighborIndices(vertexCount);
    for (int i = 0; i < vertexCount; i++) {
        int lat = i / (longitudeBands + 1);
        int lon = i % (longitudeBands + 1);
        std::vector<int> neighbors;
        int offsets[8][2] = {
            {-1,-1},{0,-1},{1,-1},
            {-1, 0},       {1, 0},
            {-1, 1},{0, 1},{1, 1}
        };
        for (auto& off : offsets) {
            int nLat = lat + off[0];
            int nLon = lon + off[1];
            if (nLon < 0) nLon = longitudeBands;
            else if (nLon > longitudeBands) nLon = 0;
            if (nLat >= 0 && nLat <= latitudeBands) {
                int idx = nLat * (longitudeBands + 1) + nLon;
                neighbors.push_back(idx);
            }
        }
        neighborIndices[i] = neighbors;
    }

    // flow directions
    std::vector<int> flowDirections(vertexCount, -1);

    for (int i = 0; i < vertexCount; i++) {
        if (highVertices[i].elevation <= 0.f) continue;
        int currentIdx = i;
        int steps = 0;
        while (steps < 100) {
            int lowestNeighbor = -1;
            float lowestElevation = highVertices[currentIdx].elevation;
            for (int nbrIdx : neighborIndices[currentIdx]) {
                if (highVertices[nbrIdx].elevation < lowestElevation) {
                    lowestElevation = highVertices[nbrIdx].elevation;
                    lowestNeighbor = nbrIdx;
                }
            }
            if (lowestNeighbor == -1 || highVertices[currentIdx].elevation <= 0.f) {
                break;
            }
            flow[lowestNeighbor] += 1.f;
            flowDirections[currentIdx] = lowestNeighbor;
            currentIdx = lowestNeighbor;
            steps++;
        }
    }
    float flowThreshold = 30.f;
    for (int i = 0; i < vertexCount; i++) {
        if (flow[i] > flowThreshold && highVertices[i].elevation > 0.f) {
            highVertices[i].isRiver = true;
        }
    }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//        Province & Resource Generation (High-Res Only)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void GoldilocksPlanet::generateProvinces()
{
    int numberOfProvinces = 500;
    provinces.resize(numberOfProvinces);
    for (int i = 0; i < numberOfProvinces; i++) {
        provinces[i].id = i;
        provinces[i].type = ProvinceType::Inland;
        provinces[i].colorR = 0.f;
        provinces[i].colorG = 0.f;
        provinces[i].colorB = 0.f;
        provinces[i].colorA = 0.5f;
    }

    std::vector<bool> assigned(highVertices.size(), false);
    int currentProvinceId = 0;

    // seed initial province centers
    for (int i = 0; i < numberOfProvinces; i++) {
        int seedIndex = rng() % highVertices.size();
        if (!assigned[seedIndex]) {
            provinces[i].vertexIndices.push_back(seedIndex);
            highVertices[seedIndex].provinceId = i;
            assigned[seedIndex] = true;
            currentProvinceId++;
        }
    }

    std::vector<std::queue<int>> provinceQueues(numberOfProvinces);
    for (const Province& p : provinces) {
        for (int vIdx : p.vertexIndices) {
            provinceQueues[p.id].push(vIdx);
        }
    }

    // Precompute neighbors in smaller resolution for BFS
    // (But we must be consistent with high-res)
    int longitudeBands = 512;
    int latitudeBands = 512;

    // NOTE: If your planet is actually 1000×1000, you might want
    // to unify these or do a different approach. This example is simplified.
    std::vector<std::vector<int>> neighbors(highVertices.size());
    for (int i = 0; i < (int)highVertices.size(); i++) {
        int lat = i / (longitudeBands + 1);
        int lon = i % (longitudeBands + 1);
        int offs[4][2] = { {-1,0},{0,1},{1,0},{0,-1} };
        for (auto& o : offs) {
            int nLat = lat + o[0];
            int nLon = lon + o[1];
            if (nLon < 0) nLon = longitudeBands;
            else if (nLon > longitudeBands) nLon = 0;
            if (nLat >= 0 && nLat <= latitudeBands) {
                int neighborIdx = nLat * (longitudeBands + 1) + nLon;
                neighbors[i].push_back(neighborIdx);
            }
        }
    }

    bool expansionPossible = true;
    while (expansionPossible) {
        expansionPossible = false;
        for (int i = 0; i < numberOfProvinces; i++) {
            if (provinceQueues[i].empty()) continue;
            int currentVertex = provinceQueues[i].front();
            provinceQueues[i].pop();
            for (int nb : neighbors[currentVertex]) {
                if (!assigned[nb]) {
                    assigned[nb] = true;
                    highVertices[nb].provinceId = i;
                    provinces[i].vertexIndices.push_back(nb);
                    provinceQueues[i].push(nb);
                    expansionPossible = true;
                }
            }
        }
    }
}

void GoldilocksPlanet::classifyProvinceTypes()
{
    // Example logic that references highVertices
    for (auto& province : provinces) {
        int oceanCount = 0, seaCoastalCount = 0, landCoastalCount = 0, islandCount = 0, inlandCount = 0, mountainCount = 0;

        for (auto idx : province.vertexIndices) {
            if (highVertices[idx].isOcean) {
                oceanCount++;
            }
            else {
                bool isCoastal = false;
                int longitudeBands = 512;
                int latitudeBands = 512;
                int lat = idx / (longitudeBands + 1);
                int lon = idx % (longitudeBands + 1);
                int offs[4][2] = { {-1,0},{0,1},{1,0},{0,-1} };
                for (auto& o : offs) {
                    int nLat = lat + o[0];
                    int nLon = lon + o[1];
                    if (nLon < 0)nLon = longitudeBands;
                    else if (nLon > longitudeBands)nLon = 0;
                    if (nLat >= 0 && nLat <= latitudeBands) {
                        int nIdx = nLat * (longitudeBands + 1) + nLon;
                        if (highVertices[nIdx].isOcean) {
                            isCoastal = true;
                            break;
                        }
                    }
                }
                if (isCoastal) {
                    // ...
                    if (highVertices[idx].elevation < getRadius() * 0.02f)
                        seaCoastalCount++;
                    else
                        landCoastalCount++;
                }
                else {
                    // ...
                    int oceanNeighbors = 0;
                    for (auto& o : offs) {
                        int nLat = lat + o[0];
                        int nLon = lon + o[1];
                        if (nLon < 0)nLon = longitudeBands;
                        else if (nLon > longitudeBands)nLon = 0;
                        if (nLat >= 0 && nLat <= latitudeBands) {
                            int nIdx = nLat * (longitudeBands + 1) + nLon;
                            if (highVertices[nIdx].isOcean)
                                oceanNeighbors++;
                        }
                    }
                    if (oceanNeighbors == 0) inlandCount++;
                    else islandCount++;
                }
                // mountain check
                if (highVertices[idx].elevation > getRadius() * 0.1f) {
                    mountainCount++;
                }
            }
        }

        // final classification
        int total = (int)province.vertexIndices.size();
        if (oceanCount > total / 2)
            province.type = ProvinceType::Ocean;
        else if (seaCoastalCount > 0 && landCoastalCount > 0)
            province.type = ProvinceType::SeaCoastal;
        else if (landCoastalCount > 0)
            province.type = ProvinceType::LandCoastal;
        else if (islandCount > 0)
            province.type = ProvinceType::Island;
        else if (mountainCount > total / 3)
            province.type = ProvinceType::Mountain;
        else
            province.type = ProvinceType::Inland;
    }
}

void GoldilocksPlanet::assignProvinceColors()
{
    std::map<ProvinceType, std::vector<float>> typeColors = {
        { ProvinceType::Ocean,     {0.f,0.f,0.5f,0.5f} },
        { ProvinceType::SeaCoastal,{0.f,0.5f,0.8f,0.5f} },
        { ProvinceType::LandCoastal,{0.4f,0.8f,0.4f,0.5f}},
        { ProvinceType::Island,    {0.8f,0.8f,0.f,0.5f} },
        { ProvinceType::Inland,    {0.6f,0.3f,0.f,0.5f} },
        { ProvinceType::Mountain,  {0.5f,0.5f,0.5f,0.5f}}
    };
    for (auto& p : provinces) {
        if (typeColors.find(p.type) != typeColors.end()) {
            p.colorR = typeColors[p.type][0];
            p.colorG = typeColors[p.type][1];
            p.colorB = typeColors[p.type][2];
            p.colorA = typeColors[p.type][3];
        }
        else {
            p.colorR = 1.f; p.colorG = 1.f; p.colorB = 1.f; p.colorA = 0.5f;
        }
    }
}

void GoldilocksPlanet::renderProvinces()
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Draw provinces as solid
    for (const auto& prov : provinces) {
        glColor4f(prov.colorR, prov.colorG, prov.colorB, prov.colorA);
        glBegin(GL_TRIANGLES);
        for (auto vIdx : prov.vertexIndices) {
            auto& v = highVertices[vIdx];
            glNormal3f(v.nx, v.ny, v.nz);
            glVertex3f(v.x, v.y, v.z);
        }
        glEnd();
    }

    // Province borders
    glLineWidth(1.5f);
    glColor3f(0.f, 0.f, 0.f);
    glBegin(GL_LINES);
    for (const auto& prov : provinces) {
        for (auto vIdx : prov.vertexIndices) {
            int longitudeBands = 512;
            int latitudeBands = 512;
            int lat = vIdx / (longitudeBands + 1);
            int lon = vIdx % (longitudeBands + 1);
            int offs[4][2] = { {-1,0},{0,1},{1,0},{0,-1} };
            for (auto& o : offs) {
                int nLat = lat + o[0];
                int nLon = lon + o[1];
                if (nLon < 0)nLon = longitudeBands;
                else if (nLon > longitudeBands)nLon = 0;
                if (nLat >= 0 && nLat <= latitudeBands) {
                    int neighborIdx = nLat * (longitudeBands + 1) + nLon;
                    if (highVertices[neighborIdx].provinceId != prov.id) {
                        auto& v1 = highVertices[vIdx];
                        auto& v2 = highVertices[neighborIdx];
                        glVertex3f(v1.x, v1.y, v1.z);
                        glVertex3f(v2.x, v2.y, v2.z);
                    }
                }
            }
        }
    }
    glEnd();

    glDisable(GL_BLEND);
}

void GoldilocksPlanet::generateResourceNodes()
{
    int numResources = 700;
    for (int i = 0; i < numResources; i++) {
        ResourceNode node;
        int vIdx = rng() % highVertices.size();
        if (highVertices[vIdx].isOcean) continue;
        auto& v = highVertices[vIdx];
        auto biome = v.biome;
        // example logic
        if (biome == BiomeType::TropicalRainforest ||
            biome == BiomeType::TropicalDryForest ||
            biome == BiomeType::Savanna ||
            biome == BiomeType::TemperateForest ||
            biome == BiomeType::SubtropicalForest ||
            biome == BiomeType::TemperateMontaneForest ||
            biome == BiomeType::BorealForest ||
            biome == BiomeType::Taiga ||
            biome == BiomeType::Tundra ||
            biome == BiomeType::Ice)
        {
            node.positionX = v.x;
            node.positionY = v.y;
            node.positionZ = v.z;
            float chance = distribution(rng);
            // quick distribution logic
            if (biome == BiomeType::TropicalRainforest ||
                biome == BiomeType::TemperateForest ||
                biome == BiomeType::SubtropicalForest)
            {
                if (chance < 0.2f) node.type = ResourceType::Wood;
                else if (chance < 0.4f) node.type = ResourceType::Iron;
                else if (chance < 0.5f) node.type = ResourceType::Copper;
                else if (chance < 0.6f) node.type = ResourceType::Gold;
                else node.type = ResourceType::Zinc;
            }
            else if (biome == BiomeType::Savanna || biome == BiomeType::Steppe)
            {
                if (chance < 0.3f) node.type = ResourceType::Coal;
                else if (chance < 0.5f) node.type = ResourceType::Copper;
                else if (chance < 0.7f) node.type = ResourceType::Tin;
                else node.type = ResourceType::Lead;
            }
            else if (biome == BiomeType::BorealForest ||
                biome == BiomeType::Taiga ||
                biome == BiomeType::Tundra ||
                biome == BiomeType::Ice)
            {
                if (chance < 0.25f) node.type = ResourceType::Coal;
                else if (chance < 0.45f) node.type = ResourceType::Lead;
                else if (chance < 0.6f) node.type = ResourceType::Zinc;
                else node.type = ResourceType::Tin;
            }
            resourceNodes.push_back(node);
        }
    }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                         Updates & Seasons
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void GoldilocksPlanet::update(float movementSpeed)
{
    Planet::update(movementSpeed);

    orbitalPosition += orbitSpeed * movementSpeed;
    if (orbitalPosition >= 360.f) orbitalPosition -= 360.f;

    updateSeasons(movementSpeed);
}

void GoldilocksPlanet::updateSeasons(float movementSpeed)
{
    if (orbitalPosition < 90.f) currentSeason = Season::Spring;
    else if (orbitalPosition < 180.f) currentSeason = Season::Summer;
    else if (orbitalPosition < 270.f) currentSeason = Season::Autumn;
    else currentSeason = Season::Winter;

    calculateSeasonalEffects();
}

void GoldilocksPlanet::calculateSeasonalEffects()
{
    float seasonAngleRad = orbitalPosition * (float)M_PI / 180.f;
    float axialTiltRad = axialTilt * (float)M_PI / 180.f;

    for (auto& v : highVertices) {
        float declination = asinf(sinf(axialTiltRad) * sinf(seasonAngleRad));
        float insolation = cosf(v.latitude - declination);
        float elevFactor = exp(-v.elevation / (getRadius() * 0.005f));
        v.temperature = insolation * elevFactor;

        switch (currentSeason) {
        case Season::Spring: v.temperature += 0.05f; break;
        case Season::Summer: v.temperature += 0.15f; break;
        case Season::Autumn: v.temperature -= 0.05f; break;
        case Season::Winter: v.temperature -= 0.15f; break;
        }
        if (v.temperature < 0.f) v.temperature = 0.f;
        if (v.temperature > 1.f) v.temperature = 1.f;

        v.biome = determineBiome(v);
    }

    for (auto& v : highVertices) {
        switch (v.biome)
        {
        case BiomeType::TropicalRainforest:
        case BiomeType::TropicalDryForest:
        case BiomeType::Savanna:
        case BiomeType::TemperateForest:
        case BiomeType::SubtropicalForest:
        case BiomeType::TemperateMontaneForest:
            if (currentSeason == Season::Winter) v.weatherSeverity = 0.5f;
            else if (currentSeason == Season::Summer) v.weatherSeverity = 0.7f;
            else v.weatherSeverity = 0.3f;
            break;

        case BiomeType::Desert:
        case BiomeType::Steppe:
            v.weatherSeverity = 0.1f;
            break;

        case BiomeType::BorealForest:
        case BiomeType::Taiga:
        case BiomeType::Tundra:
        case BiomeType::Ice:
            if (currentSeason == Season::Winter) v.weatherSeverity = 0.8f;
            else v.weatherSeverity = 0.4f;
            break;

        case BiomeType::Beach:
            v.weatherSeverity = 0.3f;
            break;

        default:
            v.weatherSeverity = 0.5f;
        }
    }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                  Biome & Perlin Helpers
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GoldilocksPlanet::BiomeType GoldilocksPlanet::determineBiome(Vertex& v)
{
    if (v.isRiver) return BiomeType::River;
    if (v.isOcean) {
        if (v.elevation < -getRadius() * 0.02f) return BiomeType::DeepOcean;
        else return BiomeType::Ocean;
    }

    // Check beach:
    if (v.elevation >= 0.f && v.elevation < getRadius() * 0.005f) {
        return BiomeType::Beach;
    }

    float temperature = v.temperature;
    float moisture = v.moisture;
    float elevation = v.elevation;

    bool isMountain = (elevation > getRadius() * 0.1f);
    // bool isHighMountain= (elevation> getRadius()*0.15f);

    // example logic
    if (temperature > 0.8f) {
        if (moisture > 0.7f) {
            return isMountain ? BiomeType::TropicalMontaneForest : BiomeType::TropicalRainforest;
        }
        else if (moisture > 0.4f) {
            return isMountain ? BiomeType::TropicalDryForest : BiomeType::Savanna;
        }
        else return BiomeType::Desert;
    }
    else if (temperature > 0.6f) {
        if (moisture > 0.6f) {
            return isMountain ? BiomeType::SubtropicalMontaneForest : BiomeType::SubtropicalForest;
        }
        else if (moisture > 0.3f) {
            return BiomeType::Grassland;
        }
        else return BiomeType::Desert;
    }
    else if (temperature > 0.3f) {
        if (moisture > 0.5f) {
            return isMountain ? BiomeType::TemperateMontaneForest : BiomeType::TemperateForest;
        }
        else if (moisture > 0.2f) {
            return BiomeType::Grassland;
        }
        else {
            return BiomeType::Steppe;
        }
    }
    else {
        if (moisture > 0.4f) {
            return isMountain ? BiomeType::Taiga : BiomeType::BorealForest;
        }
        else if (moisture > 0.2f) {
            return BiomeType::Tundra;
        }
        else return BiomeType::Ice;
    }
}

float GoldilocksPlanet::getElevation(float x, float y, float z)
{
    // fractal Perlin
    float elevation = 0.f;
    float frequency = 0.6f;
    float amplitude = 0.06f;
    float persistence = 0.8f;
    float lacunarity = 2.5f;
    int octaves = 5;

    for (int i = 0; i < octaves; i++) {
        float noiseVal = perlinNoise(x * frequency, y * frequency, z * frequency);
        elevation += noiseVal * amplitude;
        amplitude *= persistence;
        frequency *= lacunarity;
    }
    elevation = elevation / (1.f - persistence);
    elevation *= getRadius() * 0.015f;
    return elevation;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                      Perlin Noise
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

float GoldilocksPlanet::perlinNoise(float x, float y, float z)
{
    int X = (int)floor(x) & 255;
    int Y = (int)floor(y) & 255;
    int Z = (int)floor(z) & 255;

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

float GoldilocksPlanet::fade(float t)
{
    return t * t * t * (t * (t * 6.f - 15.f) + 10.f);
}

float GoldilocksPlanet::lerp(float a, float b, float t)
{
    return a + t * (b - a);
}

float GoldilocksPlanet::grad(int hash, float x, float y, float z)
{
    int h = hash & 15;
    float u = (h < 8) ? x : y;
    float v = (h < 4) ? y : ((h == 12 || h == 14) ? x : z);
    float res = ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
    return res;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                           Render
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void GoldilocksPlanet::render()
{
    glPushMatrix();
    glTranslatef(getPositionX(), 0.f, getPositionZ());
    glRotatef(getRotationY(), 0.f, 1.f, 0.f);

    // 1) Compute camera distance
    float dx = camX - getPositionX();
    float dy = camY - 0.f;
    float dz = camZ - getPositionZ();
    float dist = sqrtf(dx * dx + dy * dy + dz * dz);

    // 2) Choose LOD
    int lodIndex = chooseLOD(dist);

    // 3) Render the chosen LOD
    if (lodIndex < 2) {
        // Use the lower LOD
        renderLODLevel(lodLevels[lodIndex]);
    }
    else {
        // Use high resolution
        glDisable(GL_TEXTURE_2D);
        glBegin(GL_TRIANGLES);
        for (unsigned int i = 0; i < highIndices.size(); i++) {
            Vertex& v = highVertices[highIndices[i]];

            // Color by biome
            switch (v.biome) {
            case BiomeType::TropicalRainforest:
                glColor3f(0.f, 0.6f, 0.f); break;
            case BiomeType::TropicalDryForest:
                glColor3f(0.4f, 0.5f, 0.f); break;
            case BiomeType::Savanna:
                if (currentSeason == Season::Summer)
                    glColor3f(0.9f, 0.8f, 0.2f);
                else
                    glColor3f(0.8f, 0.7f, 0.1f);
                break;
            case BiomeType::TemperateForest:
            case BiomeType::SubtropicalForest:
            case BiomeType::TemperateMontaneForest:
                if (currentSeason == Season::Autumn)
                    glColor3f(0.8f, 0.5f, 0.f);
                else if (currentSeason == Season::Winter)
                    glColor3f(0.5f, 0.5f, 0.5f);
                else
                    glColor3f(0.f, 0.5f, 0.f);
                break;
            case BiomeType::Desert:
                glColor3f(0.9f, 0.8f, 0.5f); break;
            case BiomeType::Grassland:
                if (currentSeason == Season::Winter)
                    glColor3f(0.6f, 0.6f, 0.3f);
                else if (currentSeason == Season::Autumn)
                    glColor3f(0.7f, 0.8f, 0.3f);
                else
                    glColor3f(0.4f, 0.8f, 0.4f);
                break;
            case BiomeType::Steppe:
                glColor3f(0.6f, 0.6f, 0.3f); break;
            case BiomeType::BorealForest:
            case BiomeType::Taiga:
                glColor3f(0.f, 0.4f, 0.f); break;
            case BiomeType::Tundra:
                glColor3f(0.6f, 0.6f, 0.6f); break;
            case BiomeType::Ice:
                if (currentSeason == Season::Winter)
                    glColor3f(0.9f, 0.9f, 0.9f);
                else
                    glColor3f(0.7f, 0.7f, 0.7f);
                break;
            case BiomeType::River:
                glColor3f(0.f, 0.5f, 1.f); break;
            case BiomeType::Ocean:
                glColor3f(0.f, 0.f, 0.5f); break;
            case BiomeType::DeepOcean:
                glColor3f(0.f, 0.f, 0.3f); break;
            case BiomeType::Beach:
                glColor3f(0.93f, 0.85f, 0.58f); break;
            default:
                glColor3f(1.f, 1.f, 1.f); break;
            }

            glNormal3f(v.nx, v.ny, v.nz);
            glVertex3f(v.x, v.y, v.z);
        }
        glEnd();
        glEnable(GL_TEXTURE_2D);

        // Optionally show provinces or resource nodes if toggled
        if (showResourceNodes) {
            glDisable(GL_TEXTURE_2D);
            glPointSize(5.f);
            glBegin(GL_POINTS);
            for (auto& node : resourceNodes) {
                switch (node.type) {
                case ResourceType::Coal:  glColor3f(0.1f, 0.1f, 0.1f); break;
                case ResourceType::Wood:  glColor3f(0.4f, 0.2f, 0.f);  break;
                case ResourceType::Iron:  glColor3f(0.5f, 0.5f, 0.5f); break;
                case ResourceType::Copper:glColor3f(0.72f, 0.45f, 0.2f); break;
                case ResourceType::Tin:   glColor3f(0.82f, 0.82f, 0.82f); break;
                case ResourceType::Gold:  glColor3f(1.f, 0.84f, 0.f);   break;
                case ResourceType::Zinc:  glColor3f(0.6f, 0.6f, 0.6f); break;
                case ResourceType::Lead:  glColor3f(0.3f, 0.3f, 0.3f); break;
                default: glColor3f(1.f, 1.f, 1.f); break;
                }
                glVertex3f(node.positionX, node.positionY, node.positionZ);
            }
            glEnd();
            glEnable(GL_TEXTURE_2D);
        }
        else if (showProvinces) {
            renderProvinces();
        }
    }

    // Render the atmosphere
    glPushMatrix();
    glRotatef(rotationY + 5.f, 0.f, 1.f, 0.f);
    glBindTexture(GL_TEXTURE_2D, atmosphereTextureID);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);

    glColor4f(1.f, 1.f, 1.f, 0.5f);
    renderSphere(atmosphereRadius, 40, 40);

    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    glColor4f(1.f, 1.f, 1.f, 1.f);

    glPopMatrix(); // end atmosphere

    // If there's a moon
    if (moon) moon->render();

    glEnable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);

    glPopMatrix();
}
