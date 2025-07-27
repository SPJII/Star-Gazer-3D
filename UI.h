#ifndef UI_H
#define UI_H

#include <SDL.h>
#include <vector>
#include <string>
#include "CelestialBodies.h"  // We assume you have a base class CelestialBody, etc.
#include "Stars.h"            // If you need it

// Screen dimensions
extern const int SCREEN_WIDTH;
extern const int SCREEN_HEIGHT;

// Zoom limits
extern const float MIN_ZOOM;
extern const float MAX_ZOOM;

// UI Constants
extern const int TAB_WIDTH;
extern const int TAB_HEIGHT;
extern const int SLIDER_WIDTH;
extern const int SLIDER_HEIGHT;
extern const int SLIDER_KNOB_WIDTH;
extern const int SLIDER_KNOB_HEIGHT;

// Global movement speed variable
extern float movementSpeed;

// Mouse controls and variables
extern bool dragging;
extern int lastMouseX, lastMouseY;

// UI variables
extern bool tabOpen;
extern bool showResourceNodes;
extern bool showProvinces;
extern float sliderPosition;
extern bool draggingSlider;

// Camera variables
extern float camX, camY, camZ;       // Camera position
extern float camYaw, camPitch;       // Camera orientation
extern bool cameraLocked;
extern CelestialBody* lockedObject;

// Camera follow relative position variables
extern float relativeDistance;
extern float relativeYaw;
extern float relativePitch;

// Camera transition variables
extern bool isTransitioning;
extern float targetCamX, targetCamY, targetCamZ;
extern float transitionSpeed;

// Define default camera positions and orientation for reset
extern float defaultCamX, defaultCamY, defaultCamZ;
extern float defaultYaw, defaultPitch;

// Movement variables
extern bool keys[SDL_NUM_SCANCODES];

// Map mode, if you need multiple display modes
enum MapMode {
    NORMAL_MODE = 0,
    HEIGHT_MAP_MODE,
    HEAT_MAP_MODE,
    SEASON_MAP_MODE,
    BIOME_MAP_MODE
};
extern int currentMapMode;

// ---------- Functions ----------

// Main input function
void handleInput(SDL_Event& event,
    bool& running,
    std::vector<CelestialBody*>& allCelestialBodies,
    bool& movementFrozen,
    Sun& sun);

void handleCameraMovement();
void setupCamera();
void renderUI();              // Renders any 2D UI on top of the 3D scene
void renderText(float x, float y, const std::string& text, void* font);

// For 3D picking
void unProject(int mouseX, int mouseY, float* rayOrigin, float* rayDir);

// The mini-map
void renderMiniMap(const std::vector<CelestialBody*>& bodies);
extern std::vector<CelestialBody*>& getAllCelestialBodies();

#endif // UI_H
