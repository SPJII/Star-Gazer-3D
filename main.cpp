// main.cpp
#include <iostream>
#include <SDL.h>
#include <SDL_image.h>
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <cstdlib>
#include <cmath>
#include <GL/glut.h>
#include <chrono>
#include <vector>
#include <string>

// Include your custom headers
#include "CelestialBodies.h"    // For Sun, Planet, Moon, etc.
#include "GoldilocksPlanet.h"
#include "DeadPlanets.h"       // For TerrestrialPlanet, GasGiantPlanet
#include "Stars.h"             // For class Stars
#include "UI.h"                // For handleInput, handleCameraMovement, etc.

// Universal container for celestial bodies.
static std::vector<CelestialBody*> allCelestialBodies;
std::vector<CelestialBody*>& getAllCelestialBodies() {
    return allCelestialBodies;
}

// Global variable to store the Sun’s radius (set after sun creation)
float globalSunRadius = 0.0f;

// Constants
const float DISTANCE_SCALE = 8.0f;
extern const int SCREEN_WIDTH;  // defined in UI.cpp
extern const int SCREEN_HEIGHT;

void initSDL(int argc, char* argv[], SDL_Window*& window, SDL_GLContext& context);
void initOpenGL();
GLuint loadTexture(const char* filename);
void cleanup(SDL_Window* window, SDL_GLContext context);

// --------------------- Frustum Culling Helpers ---------------------
struct Plane {
    float a, b, c, d;
};

void normalizePlane(Plane& p) {
    float mag = sqrt(p.a * p.a + p.b * p.b + p.c * p.c);
    p.a /= mag;
    p.b /= mag;
    p.c /= mag;
    p.d /= mag;
}

void extractFrustumPlanes(Plane frustum[6]) {
    float proj[16], modl[16], clip[16];
    glGetFloatv(GL_PROJECTION_MATRIX, proj);
    glGetFloatv(GL_MODELVIEW_MATRIX, modl);

    // Multiply modelview and projection matrices.
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            clip[i * 4 + j] =
                modl[i * 4 + 0] * proj[0 * 4 + j] +
                modl[i * 4 + 1] * proj[1 * 4 + j] +
                modl[i * 4 + 2] * proj[2 * 4 + j] +
                modl[i * 4 + 3] * proj[3 * 4 + j];
        }
    }

    // Extract planes: RIGHT, LEFT, BOTTOM, TOP, FAR, NEAR.
    frustum[0].a = clip[3] - clip[0];
    frustum[0].b = clip[7] - clip[4];
    frustum[0].c = clip[11] - clip[8];
    frustum[0].d = clip[15] - clip[12];
    normalizePlane(frustum[0]);

    frustum[1].a = clip[3] + clip[0];
    frustum[1].b = clip[7] + clip[4];
    frustum[1].c = clip[11] + clip[8];
    frustum[1].d = clip[15] + clip[12];
    normalizePlane(frustum[1]);

    frustum[2].a = clip[3] + clip[1];
    frustum[2].b = clip[7] + clip[5];
    frustum[2].c = clip[11] + clip[9];
    frustum[2].d = clip[15] + clip[13];
    normalizePlane(frustum[2]);

    frustum[3].a = clip[3] - clip[1];
    frustum[3].b = clip[7] - clip[5];
    frustum[3].c = clip[11] - clip[9];
    frustum[3].d = clip[15] - clip[13];
    normalizePlane(frustum[3]);

    frustum[4].a = clip[3] - clip[2];
    frustum[4].b = clip[7] - clip[6];
    frustum[4].c = clip[11] - clip[10];
    frustum[4].d = clip[15] - clip[14];
    normalizePlane(frustum[4]);

    frustum[5].a = clip[3] + clip[2];
    frustum[5].b = clip[7] + clip[6];
    frustum[5].c = clip[11] + clip[10];
    frustum[5].d = clip[15] + clip[14];
    normalizePlane(frustum[5]);
}

bool sphereInFrustum(Plane frustum[6], float x, float y, float z, float radius) {
    for (int i = 0; i < 6; i++) {
        if (frustum[i].a * x + frustum[i].b * y + frustum[i].c * z + frustum[i].d <= -radius)
            return false;
    }
    return true;
}

// -------------------------
// GALAXY FACTORY
// -------------------------
// Spawn a galaxy using the Stars factory.
Stars spawnGalaxy(int starCount, float radius, Stars::GalaxyType type, float x, float y, float z) {
    return Stars(starCount, radius, type, x, y, z);
}

// Wrap a Stars object with its center and a bounding radius.
struct Galaxy {
    Stars stars;
    float centerX;
    float centerY;
    float centerZ;
    float boundingRadius;

    Galaxy(Stars&& s, float cx, float cy, float cz, float br)
        : stars(std::move(s)), centerX(cx), centerY(cy), centerZ(cz), boundingRadius(br) {}

    Galaxy(Galaxy&&) = default;
    Galaxy& operator=(Galaxy&&) = default;
    Galaxy(const Galaxy&) = delete;
    Galaxy& operator=(const Galaxy&) = delete;
};

std::vector<Galaxy> generateBackgroundGalaxies(int count = 700) {
    std::vector<Galaxy> galaxies;
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> distDistance(50000.0f, 200000.0f);
    std::uniform_real_distribution<float> distGalaxyRadius(100000.0f, 5000000.0f);
    std::uniform_real_distribution<float> distTheta(0.0f, 2.0f * M_PI);
    std::uniform_real_distribution<float> distCosPhi(-1.0f, 1.0f);
    std::uniform_int_distribution<int> distType(0, 4);

    for (int i = 0; i < count; i++) {
        int starCount = 3000;
        float galaxyRadius = distGalaxyRadius(rng);
        Stars::GalaxyType type;
        switch (distType(rng)) {
        case 0: type = Stars::ELLIPTICAL; break;
        case 1: type = Stars::SPIRAL; break;
        case 2: type = Stars::BARRED_SPIRAL; break;
        case 3: type = Stars::LENTICULAR; break;
        case 4: type = Stars::IRREGULAR; break;
        default: type = Stars::SPIRAL; break;
        }
        float theta = distTheta(rng);
        float cosPhi = distCosPhi(rng);
        float phi = acos(cosPhi);
        float distance = distDistance(rng);
        float x = distance * sin(phi) * cos(theta);
        float y = distance * sin(phi) * sin(theta);
        float z = distance * cos(phi);

        Stars galaxyStars = spawnGalaxy(starCount, galaxyRadius, type, x, y, z);
        float cullingRadius = std::min(galaxyRadius, 100000.0f);
        galaxies.push_back(Galaxy(std::move(galaxyStars), x, y, z, cullingRadius));
    }
    return galaxies;
}


// -------------------------
// SUN FACTORY
// -------------------------
Sun* spawnSun(SunType sunType, float x, float y, GLuint sunTexture) {
    Sun* sun = new Sun(sunType, x, y, sunTexture);
    allCelestialBodies.push_back(sun);
    globalSunRadius = sun->getRadius(); // store for later use in planet factory and camera clamping
    return sun;
}

// -------------------------
// PLANET FACTORY
// -------------------------
Planet* spawnPlanet(const std::string& planetType,
    float orbitAU,
    float radiusEarthUnits,
    float distScale,
    // textures
    GLuint terrTex,
    GLuint atmosTex,
    GLuint gasTex,
    GLuint moonTex,
    GLuint moonAtmos)
{
    float orbitRadius = orbitAU * 150.0f * distScale;
    // Enforce a minimum orbit distance relative to the sun (if available)
    float minOrbit = (globalSunRadius > 0.0f) ? globalSunRadius * 1.2f : 200.0f;
    if (orbitRadius < minOrbit) {
        orbitRadius = minOrbit;
    }
    float planetRadius = radiusEarthUnits;
    float orbitSpeed = 0.1f / orbitAU; // example approach

    Planet* planet = nullptr;
    if (planetType == "Terrestrial") {
        auto* tp = new TerrestrialPlanet(
            planetRadius,
            terrTex,
            atmosTex,
            orbitRadius,
            orbitSpeed
        );
        // Create 0 to 2 moons.
        int moonCount = rand() % 3;
        for (int i = 0; i < moonCount; i++) {
            float mDist = planetRadius * (5.f + i * 2.f);
            float mSize = planetRadius * 0.1f;
            auto* moon = new Moon(tp, mDist, mSize, moonTex, moonAtmos);
            moon->setOrbitSpeed(0.2f + 0.1f * i);
            tp->addMoon(moon);
            allCelestialBodies.push_back(moon);
        }
        planet = tp;
    }
    else if (planetType == "Goldilocks") {
        float mDist = planetRadius * 60.f;
        float mSize = planetRadius * 0.27f;
        auto* gp = new GoldilocksPlanet(
            planetRadius,
            planetRadius * 1.1f,
            terrTex,
            atmosTex,
            nullptr,
            orbitRadius,
            orbitSpeed
        );
        float ascendAngle = static_cast<float>(rand() % 360); // Random ascending node angle.
        Moon* earthMoon = new Moon(gp, mDist, mSize, terrTex, atmosTex, ascendAngle);
        earthMoon->setOrbitSpeed(0.2f);
        gp->addMoon(earthMoon);
        allCelestialBodies.push_back(earthMoon);
        planet = gp;
    }
    else if (planetType == "GasGiant") {
        auto* gg = new GasGiantPlanet(
            planetRadius,
            gasTex,
            orbitRadius,
            orbitSpeed
        );
        gg->generateMoons();
        gg->generateRings();
        planet = gg;
    }
    else {
        auto* fallback = new TerrestrialPlanet(
            planetRadius,
            terrTex,
            atmosTex,
            orbitRadius,
            orbitSpeed
        );
        planet = fallback;
    }
    allCelestialBodies.push_back(planet);
    return planet;
}

// ----------------------------------------------------------
int main(int argc, char* argv[])
{
    srand(unsigned(time(0)));
    SDL_Window* window = nullptr;
    SDL_GLContext context = nullptr;

    initSDL(argc, argv, window, context);
    initOpenGL();

    if (!(IMG_Init(IMG_INIT_PNG | IMG_INIT_JPG) & (IMG_INIT_PNG | IMG_INIT_JPG))) {
        std::cerr << "SDL_image init failed: " << IMG_GetError() << std::endl;
        cleanup(window, context);
        return 1;
    }

    // Load textures.
    GLuint terrestrialTexture = loadTexture("moon.jpg");
    GLuint planetAtmosphereTexture = loadTexture("clouds.png");
    GLuint gasGiantTexture = loadTexture("HabitablePlanet.png");
    GLuint moonTexture = loadTexture("moon.jpg");
    GLuint moonAtmosphereTexture = loadTexture("clouds.png");
    GLuint sunTexture = loadTexture("sun.jpg");

    // 1) Create the Sun.
    Sun* sun = spawnSun(SunType::G, 0.f, 0.f, sunTexture);

    // 2) Create Planets from data.
    struct PlanetData {
        std::string name;
        std::string type;    // "Terrestrial", "Goldilocks", or "GasGiant"
        float orbitAU;
        float radiusEarth;
    };
    std::vector<PlanetData> planetData = {
        {"Mercury",  "Terrestrial", 0.39f, 0.38f},
        {"Venus",    "Goldilocks",  0.72f, 0.95f},
        {"Earth",    "Goldilocks",  1.00f, 1.00f},
        {"Mars",     "Goldilocks",  1.52f, 0.53f},
        {"Jupiter",  "GasGiant",    5.20f, 11.21f},
        {"Saturn",   "GasGiant",    9.58f, 9.45f},
        {"Uranus",   "GasGiant",    19.20f,4.01f},
        {"Neptune",  "GasGiant",    30.05f,3.88f}
    };

    for (auto& pd : planetData) {
        spawnPlanet(pd.type,
            pd.orbitAU,
            pd.radiusEarth,
            DISTANCE_SCALE,
            terrestrialTexture,
            planetAtmosphereTexture,
            gasGiantTexture,
            moonTexture,
            moonAtmosphereTexture);
    }



    std::vector<Galaxy> galaxies = generateBackgroundGalaxies(20);

    // MAIN LOOP.
    bool running = true;
    bool movementFrozen = false;
    SDL_Event event;
    auto startTime = std::chrono::high_resolution_clock::now();

    while (running) {
        while (SDL_PollEvent(&event)) {
            handleInput(event, running, allCelestialBodies, movementFrozen, *sun);
        }

        auto currentTime = std::chrono::high_resolution_clock::now();
        float elapsedSec = std::chrono::duration<float>(currentTime - startTime).count();

        if (!movementFrozen) {
            for (auto* body : allCelestialBodies) {
                body->update(movementSpeed); // Pass dt or a constant as needed.
            }
            // Update galaxies if needed.
        }

        handleCameraMovement();

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();
        setupCamera();

        // Ensure camera does not enter the sun: clamp distance.
        {
            // Compute distance from camera to sun.
            float dx = camX - sun->getPositionX();
            float dy = camY - sun->getPositionY();
            float dz = camZ - sun->getPositionZ();
            float dist = sqrt(dx * dx + dy * dy + dz * dz);
            float safeDistance = sun->getRadius() * 2.0f;
            if (dist < safeDistance) {
                // Move the camera back along the viewing vector.
                float factor = safeDistance / (dist + 0.0001f);
                camX = sun->getPositionX() + dx * factor;
                camY = sun->getPositionY() + dy * factor;
                camZ = sun->getPositionZ() + dz * factor;
            }
        }


        // Render background galaxies as a static sky.
                // Zero out translation so galaxies remain fixed relative to camera rotation.
        GLfloat modelview[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, modelview);
        modelview[12] = modelview[13] = modelview[14] = 0.0f;
        glPushMatrix();
        glLoadMatrixf(modelview);
        {
            Plane frustum[6];
            extractFrustumPlanes(frustum);
            for (auto& galaxy : galaxies) {
                // Use a fixed culling radius for galaxies.
                if (sphereInFrustum(frustum, galaxy.centerX, galaxy.centerY, galaxy.centerZ, galaxy.boundingRadius))
                    galaxy.stars.render(0.0f);  // Pass constant time (0.0) to disable flickering.
            }
        }
        glPopMatrix();


        glPopMatrix();
        // Render celestial bodies with frustum culling.
        {
            Plane frustum[6];
            extractFrustumPlanes(frustum);
            for (auto* body : getAllCelestialBodies()) {
                if (sphereInFrustum(frustum, body->getPositionX(), body->getPositionY(), body->getPositionZ(), body->getRadius()))
                    body->render();
            }
        }

        // Render UI.
        renderUI();

        SDL_GL_SwapWindow(window);
    }

    // Cleanup.
    for (auto* b : allCelestialBodies) {
        delete b;
    }
    allCelestialBodies.clear();

    IMG_Quit();
    cleanup(window, context);
    return 0;
}

// SDL Initialization.
void initSDL(int argc, char* argv[], SDL_Window*& window, SDL_GLContext& context) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "SDL could not initialize! SDL_Error: " << SDL_GetError() << std::endl;
        exit(1);
    }

    // Initialize GLUT for text rendering.
    glutInit(&argc, argv);

    // Set OpenGL version (using 2.1 for compatibility).
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

    window = SDL_CreateWindow("StarGazer 3D",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        SCREEN_WIDTH, SCREEN_HEIGHT,
        SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    if (!window) {
        std::cerr << "Window could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        exit(1);
    }

    context = SDL_GL_CreateContext(window);
    if (!context) {
        std::cerr << "OpenGL context could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        exit(1);
    }

    // Enable VSync.
    if (SDL_GL_SetSwapInterval(1) < 0) {
        std::cerr << "Warning: Unable to set VSync! SDL_Error: " << SDL_GetError() << std::endl;
    }
}

// OpenGL Initialization.
void initOpenGL() {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, (double)SCREEN_WIDTH / (double)SCREEN_HEIGHT, 0.1, 100000.0); // Adjusted far clipping plane.

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    // Set initial light properties.
    GLfloat lightAmbient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat lightDiffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f }; // Increased diffuse light for brightness
    GLfloat lightSpecular[] = { 1.0f, 1.0f, 1.0f, 1.0f };

    glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);

    glShadeModel(GL_SMOOTH);
    glEnable(GL_TEXTURE_2D);

    glEnable(GL_COLOR_MATERIAL);
}

// Load texture function.
GLuint loadTexture(const char* filename) {
    SDL_Surface* surface = IMG_Load(filename);
    if (!surface) {
        std::cerr << "Failed to load texture (" << filename << "): " << IMG_GetError() << std::endl;
        exit(1);
    }

    GLuint textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);

    // Set texture parameters.
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);    // Horizontal wrapping.
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);    // Vertical wrapping.
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); // Minification filter.
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); // Magnification filter.

    // Determine the format.
    GLenum format;
    if (surface->format->BytesPerPixel == 3) {
        format = GL_RGB;
    }
    else if (surface->format->BytesPerPixel == 4) {
        format = GL_RGBA;
    }
    else {
        std::cerr << "Unsupported image format for texture: " << filename << std::endl;
        SDL_FreeSurface(surface);
        exit(1);
    }

    // Upload the texture to OpenGL.
    glTexImage2D(GL_TEXTURE_2D, 0, format, surface->w, surface->h, 0, format, GL_UNSIGNED_BYTE, surface->pixels);
    SDL_FreeSurface(surface);
    return textureID;
}

// Cleanup resources.
void cleanup(SDL_Window* window, SDL_GLContext context) {
    SDL_GL_DeleteContext(context);
    SDL_DestroyWindow(window);
    SDL_Quit();
}
