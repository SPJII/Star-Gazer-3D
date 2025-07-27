// UI.cpp
// Implementation of user interface and input handling functions.

#include "UI.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cfloat> // For FLT_MAX
#include <GL/glut.h> // For text rendering
#include "resource.h"
#include "CelestialBodies.h"
#include "GoldilocksPlanet.h"
#include "DeadPlanets.h"  // For GasGiantPlanet, etc.


// Screen dimensions.
const int SCREEN_WIDTH = 1915;
const int SCREEN_HEIGHT = 1000;

// Zoom limits.
const float MIN_ZOOM = 0.5f;
const float MAX_ZOOM = 2000.0f;

// UI Constants.
const int TAB_WIDTH = 30;
const int TAB_HEIGHT = 30;
const int SLIDER_WIDTH = 200;
const int SLIDER_HEIGHT = 20;
const int SLIDER_KNOB_WIDTH = 10;
const int SLIDER_KNOB_HEIGHT = 30;
// Mini-Map location and scale
static const int MINI_MAP_WIDTH = SCREEN_WIDTH - 10;
static const int MINI_MAP_HEIGHT = SCREEN_HEIGHT - 70;
static const int MINI_MAP_X = 5; // bottom-left X
static const int MINI_MAP_Y = SCREEN_HEIGHT - MINI_MAP_HEIGHT - 35; // from top if 0=top
static const float MINI_MAP_SCALE = 1.0f / 450.0f;
// adjust so your system fits in that 200x200 box
const int MINI_MAP_TAB_WIDTH = 30;
const int MINI_MAP_TAB_HEIGHT = 30;
const int MINI_MAP_TAB_X = 5;
const int MINI_MAP_TAB_Y = MINI_MAP_Y + MINI_MAP_HEIGHT ;

// Global movement speed variable.
float movementSpeed = 1.0f; // Default speed.

// Mouse controls and variables.
bool dragging = false;
int lastMouseX = 0, lastMouseY = 0;

// UI variables.
bool tabOpen = false;
bool showResourceNodes = false; // Initialize resource node display to false.
bool showProvinces = false;     // Initialize province display to false.
float sliderPosition = 0.5f;    // Slider position ranges from 0.0f to 1.0f.
bool miniMapOpen = false; // Controls minimap visibility
bool draggingSlider = false;

// Camera variables.
float camX = 0.0f, camY = 0.0f, camZ = 2000.0f; // Camera position.
float camYaw = 0.0f, camPitch = 0.0f;             // Camera orientation.
bool cameraLocked = false;
CelestialBody* lockedObject = nullptr;

// Camera follow relative position variables.
float relativeDistance = 5.0f; // Default distance for Planet.
float relativeYaw = 0.0f;      // Horizontal angle in degrees.
float relativePitch = 20.0f;   // Vertical angle in degrees.

// Camera transition variables.
bool isTransitioning = false;
float targetCamX = 0.0f, targetCamY = 0.0f, targetCamZ = 0.0f;
float transitionSpeed = 0.05f; // Adjust for smoothness.

// Define default camera positions and orientation for reset.
float defaultCamX = 0.0f, defaultCamY = 0.0f, defaultCamZ = 250.0f;
float defaultYaw = 0.0f, defaultPitch = 0.0f; // Facing towards the origin (Sun).

// Movement variables.
bool keys[SDL_NUM_SCANCODES] = { false };

// Map mode
int currentMapMode = MapMode::NORMAL_MODE;


// ——— Toggle buttons for resources & provinces ———
const int TOGGLE_BUTTON_WIDTH = 150;
const int TOGGLE_BUTTON_HEIGHT = 20;

// right‑aligned, 10px from edge
int resourceToggleX = SCREEN_WIDTH - TOGGLE_BUTTON_WIDTH - 5;
int resourceToggleY = SCREEN_HEIGHT - TOGGLE_BUTTON_HEIGHT - 5;

// placed immediately to the left of the resource toggle
int provinceToggleX = resourceToggleX - TOGGLE_BUTTON_WIDTH - 10;
int provinceToggleY = resourceToggleY;

// state flags
bool resourcesVisible = false;
bool provincesVisible = false;

float computeMaxDistance(const std::vector<CelestialBody*>& bodies, float sunX, float sunZ) {
    float maxDist = 0.0f;
    for (auto* b : bodies) {
        float dx = b->getPositionX() - sunX;
        float dz = b->getPositionZ() - sunZ;
        float d = sqrt(dx * dx + dz * dz);
        if (d > maxDist)
            maxDist = d;
    }
    return maxDist;
}

// Function to handle user input.
void handleInput(SDL_Event& event, bool& running, std::vector<CelestialBody*>& allCelestialBodies, bool& movementFrozen, Sun& sun)
{
    switch (event.type) {
    case SDL_QUIT:
        running = false;
        break;
    case SDL_KEYDOWN:
        keys[event.key.keysym.scancode] = true;
        if (event.key.keysym.sym == SDLK_SPACE) {
            movementFrozen = !movementFrozen; // Toggle movement state.
        }
        break;
    case SDL_KEYUP:
        keys[event.key.keysym.scancode] = false;
        break;
    case SDL_MOUSEMOTION:
        if (draggingSlider) {
            int mouseX = event.motion.x;
            int sliderXPos = SCREEN_WIDTH - TAB_WIDTH - SLIDER_WIDTH;
            sliderPosition = float(mouseX - sliderXPos) / (SLIDER_WIDTH - SLIDER_KNOB_WIDTH);
            if (sliderPosition < 0.f) sliderPosition = 0.f;
            if (sliderPosition > 1.f) sliderPosition = 1.f;
            movementSpeed = 0.1f + sliderPosition * 1.9f;
        }
        else if (dragging) {
            int dx = event.motion.x - lastMouseX;
            int dy = event.motion.y - lastMouseY;
            if (cameraLocked) {
                relativeYaw += dx * 0.5f;
                relativePitch += dy * 0.5f;
                if (relativePitch > 89.0f) relativePitch = 89.0f;
                if (relativePitch < -89.0f) relativePitch = -89.0f;
            }
            else {
                camYaw -= dx * 0.1f;
                camPitch -= dy * 0.1f;
                if (camPitch > 89.0f) camPitch = 89.0f;
                if (camPitch < -89.0f) camPitch = -89.0f;
            }
            lastMouseX = event.motion.x;
            lastMouseY = event.motion.y;
        }
        break;
    case SDL_MOUSEBUTTONDOWN:
        if (event.button.button == SDL_BUTTON_LEFT) {
            int mouseX = event.button.x;
            int mouseY = event.button.y;

            // First, if tab is open, check if slider knob is clicked.
            if (tabOpen) {
                int sliderXPos = SCREEN_WIDTH - TAB_WIDTH - SLIDER_WIDTH;
                int sliderYPos = TAB_HEIGHT + 40;
                int knobX = sliderXPos + int(sliderPosition * (SLIDER_WIDTH - SLIDER_KNOB_WIDTH));
                if (mouseX >= knobX && mouseX <= knobX + SLIDER_KNOB_WIDTH &&
                    mouseY >= sliderYPos - 5 && mouseY <= sliderYPos + SLIDER_KNOB_HEIGHT - 5) {
                    draggingSlider = true;
                    lastMouseX = mouseX;
                    lastMouseY = mouseY;
                    break; // exit event handling so slider works.
                }
            }

            // Check if clicking on reset button.
            int resetBtnX = 10, resetBtnY = 10, resetBtnWidth = 100, resetBtnHeight = 30;
            if (mouseX >= resetBtnX && mouseX <= resetBtnX + resetBtnWidth &&
                mouseY >= resetBtnY && mouseY <= resetBtnY + resetBtnHeight) {
                lockedObject = nullptr;
                cameraLocked = false;
                targetCamX = defaultCamX;
                targetCamY = defaultCamY;
                targetCamZ = defaultCamZ;
                isTransitioning = true;
                camYaw = defaultYaw;
                camPitch = defaultPitch;
                break;
            }
            // Resource toggle
            if (cameraLocked &&
                lockedObject &&
                lockedObject->getType() == CelestialBodyType::PLANET &&
                mouseX >= resourceToggleX && mouseX <= resourceToggleX + TOGGLE_BUTTON_WIDTH &&
                mouseY >= resourceToggleY && mouseY <= resourceToggleY + TOGGLE_BUTTON_HEIGHT) {
                resourcesVisible = !resourcesVisible;
                break;
            }

            // Province toggle
            if (cameraLocked &&
                lockedObject &&
                lockedObject->getType() == CelestialBodyType::PLANET &&
                mouseX >= provinceToggleX && mouseX <= provinceToggleX + TOGGLE_BUTTON_WIDTH &&
                mouseY >= provinceToggleY && mouseY <= provinceToggleY + TOGGLE_BUTTON_HEIGHT) {
                provincesVisible = !provincesVisible;
                break;
            }
            // Check tab.
            if (mouseX >= SCREEN_WIDTH - TAB_WIDTH && mouseY <= TAB_HEIGHT) {
                tabOpen = !tabOpen;
                break;
            }
            // Check minimap tab click
            if (mouseX >= MINI_MAP_TAB_X && mouseX <= MINI_MAP_TAB_X + MINI_MAP_TAB_WIDTH &&
                mouseY >= MINI_MAP_TAB_Y && mouseY <= MINI_MAP_TAB_Y + MINI_MAP_TAB_HEIGHT) {
                miniMapOpen = !miniMapOpen;
                break;
            }
            // Then check if click is in mini-map.
            bool inMiniMap = (miniMapOpen && mouseX >= MINI_MAP_X && mouseX <= MINI_MAP_X + MINI_MAP_WIDTH &&
                mouseY >= MINI_MAP_Y && mouseY <= MINI_MAP_Y + MINI_MAP_HEIGHT);
            if (inMiniMap) {
                // 1) Find the Sun
                float sunX = 0.0f, sunY = 0.0f, sunZ = 0.0f;
                for (auto* b : allCelestialBodies) {
                    if (b->getType() == CelestialBodyType::SUN) {
                        sunX = b->getPositionX();
                        sunY = b->getPositionY();
                        sunZ = b->getPositionZ();
                        break;
                    }
                }

                // 2) Compute maxDist so we know how large the system is
                float maxDist = computeMaxDistance(allCelestialBodies, sunX, sunZ);
                float dynamicScale = (maxDist > 0.f)
                    ? ((MINI_MAP_WIDTH * 0.45f) / maxDist) // 0.45 instead of 0.5 to leave margins
                    : (1.0f / 1000.0f);

                // 3) Convert mini-map click to world offset from the sun
                float centerX = MINI_MAP_X + MINI_MAP_WIDTH * 0.5f;
                float centerY = MINI_MAP_Y + MINI_MAP_HEIGHT * 0.5f;
                float dx = (mouseX - centerX) / dynamicScale;
                float dz = (mouseY - centerY) / dynamicScale;

                // 4) Attempt to pick the nearest body within a 'mapThreshold'
                float mapThreshold = maxDist * 0.05f; // 5% of maxDist
                float closestDistSq = 1e30f;
                CelestialBody* picked = nullptr;
                for (auto* body : allCelestialBodies) {
                    float bx = body->getPositionX() - sunX;
                    float bz = body->getPositionZ() - sunZ;
                    float distSq = (bx - dx) * (bx - dx) + (bz - dz) * (bz - dz);
                    if (distSq < closestDistSq) {
                        closestDistSq = distSq;
                        picked = body;
                    }
                }

                // 5) If we found something within 'mapThreshold', lock on it
                if (picked && miniMapOpen && sqrtf(closestDistSq) < mapThreshold) {
                    cameraLocked = true;
                    lockedObject = picked;

                    // Similar logic from your existing code:
                    float minDist = picked->getRadius() * 1.2f;
                    switch (picked->getType()) {
                    case CelestialBodyType::SUN:
                        relativeDistance = picked->getRadius() * 1.5f;
                        break;
                    case CelestialBodyType::PLANET:
                        relativeDistance = picked->getRadius() * 5.0f;
                        break;
                    case CelestialBodyType::MOON:
                        relativeDistance = picked->getRadius() * 2.0f;
                        break;
                    default:
                        relativeDistance = picked->getRadius() * 5.0f;
                        break;
                    }
                    if (relativeDistance < minDist) relativeDistance = minDist;

                    relativeYaw = 0.f;
                    relativePitch = 20.f;
                    isTransitioning = true;
                }
                if(miniMapOpen) {
                    // If no body was found within threshold, just recenter camera
                    cameraLocked = false;
                    lockedObject = nullptr;

                    // Recenter camera to that world location 
                    // so the user can freely explore that area.
                    targetCamX = dx + sunX;
                    targetCamY = camY;  // keep same height
                    targetCamZ = dz + sunZ;
                    isTransitioning = true;

                    // 2) Compute maxDist so we know how large the system is
                    float maxDist = computeMaxDistance(allCelestialBodies, sunX, sunZ);

                    // 3) Convert mini-map click to world offset from the sun
                    float centerX = MINI_MAP_X + MINI_MAP_WIDTH * 0.5f;
                    float centerY = MINI_MAP_Y + MINI_MAP_HEIGHT * 0.5f;
                    // NOTE: Y axis is inverted in screen coords if 0=top; adapt if needed
                    float localX = mouseX - centerX;
                    float localY = mouseY - centerY;
                    float worldOffsetX = localX / dynamicScale;
                    float worldOffsetZ = localY / dynamicScale;

                    // 4) Move camera there
                    cameraLocked = false;
                    lockedObject = nullptr;
                    targetCamX = sunX + worldOffsetX;
                    targetCamY = camY;             // keep current height
                    targetCamZ = sunZ + worldOffsetZ;
                    isTransitioning = true;

                }

                // End of inMiniMap block
                break;
            }
            else {
                // Otherwise, perform 3D picking (if desired)
                float rayOrigin[3], rayDir[3];
                unProject(mouseX, mouseY, rayOrigin, rayDir);
                float tMin = FLT_MAX;
                CelestialBody* pickedObject = nullptr;
                float t;
                for (const auto& body : allCelestialBodies) {
                    if (body->intersectsRay(rayOrigin, rayDir, t) && t < tMin) {
                        tMin = t;
                        pickedObject = body;
                    }
                }
                if (pickedObject) {
                    cameraLocked = true;
                    lockedObject = pickedObject;
                    switch (lockedObject->getType()) {
                    case CelestialBodyType::SUN:
                        relativeDistance = lockedObject->getRadius() * 1.5f;
                        break;
                    case CelestialBodyType::PLANET:
                        relativeDistance = lockedObject->getRadius() * 5.0f;
                        break;
                    case CelestialBodyType::MOON:
                        relativeDistance = lockedObject->getRadius() * 2.0f;
                        break;
                    default:
                        relativeDistance = lockedObject->getRadius() * 5.0f;
                        break;
                    }
                    if (relativeDistance < lockedObject->getRadius() * 1.2f)
                        relativeDistance = lockedObject->getRadius() * 1.2f;
                    relativeYaw = 0.0f;
                    relativePitch = 20.0f;
                    isTransitioning = true;
                }
                else {
                    cameraLocked = false;
                    lockedObject = nullptr;
                    targetCamX = defaultCamX;
                    targetCamY = defaultCamY;
                    targetCamZ = defaultCamZ;
                    isTransitioning = false;
                }
            }
            dragging = true;
            lastMouseX = event.button.x;
            lastMouseY = event.button.y;
        }
        break;
    case SDL_MOUSEBUTTONUP:
        if (event.button.button == SDL_BUTTON_LEFT) {
            dragging = false;
            draggingSlider = false;
        }
        break;
    case SDL_MOUSEWHEEL:
        if (cameraLocked && lockedObject) {
            // Adjust desired distance based on zoom.
            float zoomChange = event.wheel.y * 5.0f; // Adjust sensitivity as needed.
            relativeDistance -= zoomChange;

            // Clamp relativeDistance to prevent zooming past the object's surface or exceeding max zoom.
            float minDistance = lockedObject->getRadius() * 1.2f;
            if (relativeDistance < minDistance) relativeDistance = minDistance;
            if (relativeDistance > 1000.0f) relativeDistance = 1000.0f; // Arbitrary max zoom.
        }
        else {
            // Adjust camera zoom.
            camZ -= event.wheel.y * 50.0f; // Adjust sensitivity as needed.
            if (camZ < MIN_ZOOM) camZ = MIN_ZOOM;
            if (camZ > MAX_ZOOM) camZ = MAX_ZOOM;

            // Ensure camera doesn't go inside the Sun.
            if (camZ < sun.getRadius() * 2.0f) { // Arbitrary factor to prevent being too close.
                camZ = sun.getRadius() * 2.0f;
            }
        }
        break;
    }
}

// Handle free camera movement based on key states.
void handleCameraMovement() {
    float moveSpeed = 5.0f;

    // Calculate forward and right vectors based on camera orientation.
    float radYaw = camYaw * M_PI / 180.0f;
    float radPitch = camPitch * M_PI / 180.0f;

    float forwardX = sin(radYaw) * cos(radPitch);
    float forwardY = sin(radPitch);
    float forwardZ = -cos(radYaw) * cos(radPitch);

    float rightX = cos(radYaw);
    float rightY = 0;
    float rightZ = sin(radYaw);

    // Move camera based on key inputs.
    if (keys[SDL_SCANCODE_W]) {
        camX += forwardX * moveSpeed;
        camY += forwardY * moveSpeed;
        camZ += forwardZ * moveSpeed;
    }
    if (keys[SDL_SCANCODE_S]) {
        camX -= forwardX * moveSpeed;
        camY -= forwardY * moveSpeed;
        camZ -= forwardZ * moveSpeed;
    }
    if (keys[SDL_SCANCODE_A]) {
        camX -= rightX * moveSpeed;
        camZ -= rightZ * moveSpeed;
    }
    if (keys[SDL_SCANCODE_D]) {
        camX += rightX * moveSpeed;
        camZ += rightZ * moveSpeed;
    }
    if (keys[SDL_SCANCODE_Q]) {
        camY -= moveSpeed;
    }
    if (keys[SDL_SCANCODE_E]) {
        camY += moveSpeed;
    }
}

// Render the User Interface.
void renderUI() {
    // Switch to orthographic projection.
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, SCREEN_WIDTH, SCREEN_HEIGHT, 0, -1, 1);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // Disable depth test and lighting for UI rendering.
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);


    // Render the reset button in the top-left corner.
    int resetBtnX = 5; // 0 pixels from the left.
    int resetBtnY = 5; // 0 pixels from the top.
    int resetBtnWidth = 125;
    int resetBtnHeight = 20;

    // Button background.
    glColor4f(0.2f, 0.2f, 0.2f, 0.8f); // Semi-transparent dark gray.
    glBegin(GL_QUADS);
    glVertex2i(resetBtnX, resetBtnY);
    glVertex2i(resetBtnX + resetBtnWidth, resetBtnY);
    glVertex2i(resetBtnX + resetBtnWidth, resetBtnY + resetBtnHeight);
    glVertex2i(resetBtnX, resetBtnY + resetBtnHeight);
    glEnd();

    // Button border.
    glColor3f(1.0f, 1.0f, 1.0f); // White color.
    glLineWidth(2.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2i(resetBtnX, resetBtnY);
    glVertex2i(resetBtnX + resetBtnWidth, resetBtnY);
    glVertex2i(resetBtnX + resetBtnWidth, resetBtnY + resetBtnHeight);
    glVertex2i(resetBtnX, resetBtnY + resetBtnHeight);
    glEnd();

    // Reset Button label.
    {
        glColor3f(1.0f, 1.0f, 1.0f); // White color.
        std::string resetText = "Reset Position";
        // Center the text.
        float textWidth = resetText.length() * 8.0f; // Approximate width per character.
        float textX = 2 + (resetBtnWidth - textWidth) / 2.0f;
        float textY = 15 + (resetBtnHeight - 18.0f) / 2.0f + 5.0f; // 18 is approx. height of the font.

        renderText(textX, textY, resetText, GLUT_BITMAP_HELVETICA_18);
    }
    // Render the Time tab background.

    const int timeTabX = SCREEN_WIDTH - TAB_WIDTH - 5;
    const int timeTabY = 5;

    glBegin(GL_QUADS);
    glColor4f(0.2f, 0.2f, 0.2f, 0.8f);
    glVertex2i(timeTabX, timeTabY);
    glVertex2i(timeTabX + TAB_WIDTH, timeTabY);
    glVertex2i(timeTabX + TAB_WIDTH, timeTabY + TAB_HEIGHT);
    glVertex2i(timeTabX, timeTabY + TAB_HEIGHT);
    glEnd();

    // outline
    glColor3f(1, 1, 1);
    glLineWidth(2);
    glBegin(GL_LINE_LOOP);
    glVertex2i(timeTabX, timeTabY);
    glVertex2i(timeTabX + TAB_WIDTH, timeTabY);
    glVertex2i(timeTabX + TAB_WIDTH, timeTabY + TAB_HEIGHT);
    glVertex2i(timeTabX, timeTabY + TAB_HEIGHT);
    glEnd();

    // plus/minus lines
    glBegin(GL_LINES);
    // horizontal
    glVertex2i(timeTabX + 8, timeTabY + TAB_HEIGHT / 2);
    glVertex2i(timeTabX + TAB_WIDTH - 8, timeTabY + TAB_HEIGHT / 2);
    // vertical only when closed
    if (!tabOpen) {
        glVertex2i(timeTabX + TAB_WIDTH / 2, timeTabY + 8);
        glVertex2i(timeTabX + TAB_WIDTH / 2, timeTabY + TAB_HEIGHT - 8);
    }
    glEnd();

    // If Time tab is open, render the slider and label.
    if (tabOpen) {
        int sliderX = SCREEN_WIDTH - TAB_WIDTH - SLIDER_WIDTH;
        int sliderY = TAB_HEIGHT + 40; // Adjusted Y position to make room for label.

        // Slider background.
        glColor4f(0.3f, 0.3f, 0.3f, 0.8f); // Semi-transparent darker gray.
        glBegin(GL_QUADS);
        glVertex2i(sliderX, sliderY + SLIDER_HEIGHT / 2 - 5);
        glVertex2i(sliderX + SLIDER_WIDTH, sliderY + SLIDER_HEIGHT / 2 - 5);
        glVertex2i(sliderX + SLIDER_WIDTH, sliderY + SLIDER_HEIGHT / 2 + 5);
        glVertex2i(sliderX, sliderY + SLIDER_HEIGHT / 2 + 5);
        glEnd();

        // Slider knob.
        int knobX = sliderX + static_cast<int>(sliderPosition * (SLIDER_WIDTH - SLIDER_KNOB_WIDTH));
        glColor4f(0.6f, 0.6f, 0.6f, 0.9f); // Semi-transparent light gray.
        glBegin(GL_QUADS);
        glVertex2i(knobX, sliderY - 5);
        glVertex2i(knobX + SLIDER_KNOB_WIDTH, sliderY - 5);
        glVertex2i(knobX + SLIDER_KNOB_WIDTH, sliderY + SLIDER_KNOB_HEIGHT - 5);
        glVertex2i(knobX, sliderY + SLIDER_KNOB_HEIGHT - 5);
        glEnd();

        // Render the speed label.
        std::stringstream ss;
        ss << "Speed: " << std::fixed << std::setprecision(2) << movementSpeed;
        std::string speedText = ss.str();

        glColor3f(1.0f, 1.0f, 1.0f); // White color.
        renderText(sliderX, sliderY - 20, speedText, GLUT_BITMAP_HELVETICA_18);
    }

    const float fontHeight = 12.0f;

    // only draw these toggles if locked on a planet
    if (cameraLocked && lockedObject && lockedObject->getType() == CelestialBodyType::PLANET) {
        // ——— Resources toggle button ———
        glColor4f(0.2f, 0.2f, 0.2f, 0.8f);
        glBegin(GL_QUADS);
        glVertex2i(resourceToggleX, resourceToggleY);
        glVertex2i(resourceToggleX + TOGGLE_BUTTON_WIDTH, resourceToggleY);
        glVertex2i(resourceToggleX + TOGGLE_BUTTON_WIDTH, resourceToggleY + TOGGLE_BUTTON_HEIGHT);
        glVertex2i(resourceToggleX, resourceToggleY + TOGGLE_BUTTON_HEIGHT);
        glEnd();

        // border
        glColor3f(1, 1, 1);
        glLineWidth(2);
        glBegin(GL_LINE_LOOP);
        glVertex2i(resourceToggleX, resourceToggleY);
        glVertex2i(resourceToggleX + TOGGLE_BUTTON_WIDTH, resourceToggleY);
        glVertex2i(resourceToggleX + TOGGLE_BUTTON_WIDTH, resourceToggleY + TOGGLE_BUTTON_HEIGHT);
        glVertex2i(resourceToggleX, resourceToggleY + TOGGLE_BUTTON_HEIGHT);
        glEnd();

        // ——— Resources toggle button label ———
        {
            std::string txt = resourcesVisible ? "Hide Resources" : "Show Resources";
            float textWidth = txt.length() * 8.0f;
            float textX = resourceToggleX + (TOGGLE_BUTTON_WIDTH - textWidth) / 7.0f;

            const float fontHeight = 12.0f;
            float textY = resourceToggleY
                + (TOGGLE_BUTTON_HEIGHT - fontHeight) / 2.0f
                + fontHeight;

            glColor3f(1, 1, 1);
            renderText(textX, textY, txt, GLUT_BITMAP_HELVETICA_18);
        }

        // ——— Provinces toggle button ———
        glColor4f(0.2f, 0.2f, 0.2f, 0.8f);
        glBegin(GL_QUADS);
        glVertex2i(provinceToggleX, provinceToggleY);
        glVertex2i(provinceToggleX + TOGGLE_BUTTON_WIDTH, provinceToggleY);
        glVertex2i(provinceToggleX + TOGGLE_BUTTON_WIDTH, provinceToggleY + TOGGLE_BUTTON_HEIGHT);
        glVertex2i(provinceToggleX, provinceToggleY + TOGGLE_BUTTON_HEIGHT);
        glEnd();

        // border
        glColor3f(1, 1, 1);
        glLineWidth(2);
        glBegin(GL_LINE_LOOP);
        glVertex2i(provinceToggleX, provinceToggleY);
        glVertex2i(provinceToggleX + TOGGLE_BUTTON_WIDTH, provinceToggleY);
        glVertex2i(provinceToggleX + TOGGLE_BUTTON_WIDTH, provinceToggleY + TOGGLE_BUTTON_HEIGHT);
        glVertex2i(provinceToggleX, provinceToggleY + TOGGLE_BUTTON_HEIGHT);
        glEnd();

        // ——— Provinces toggle button label ———
        {
            std::string txt = provincesVisible ? "Hide Provinces" : "Show Provinces";
            float textWidth = txt.length() * 8.0f;
            float textX = provinceToggleX + (TOGGLE_BUTTON_WIDTH - textWidth) / 5.0f;

            const float fontHeight = 12.0f;
            float textY = provinceToggleY
                + (TOGGLE_BUTTON_HEIGHT - fontHeight) / 2.0f
                + fontHeight;

            glColor3f(1, 1, 1);
            renderText(textX, textY, txt, GLUT_BITMAP_HELVETICA_18);
        }
    }

    // Determine tab position based on open/closed
    int tabX, tabY;
    if (miniMapOpen) {
        // bottom-left corner of where the map was
        tabX = MINI_MAP_X;
        tabY = MINI_MAP_Y + MINI_MAP_HEIGHT - MINI_MAP_TAB_HEIGHT;
    }
    else {
        // bottom-left corner of where the map was
        tabX = MINI_MAP_X;
        tabY = MINI_MAP_Y + MINI_MAP_HEIGHT - MINI_MAP_TAB_HEIGHT;
    }

    // draw the tab background (at the corrected position)
    glColor3f(0.2f, 0.2f, 0.2f);
    glRecti(
        MINI_MAP_TAB_X,
        MINI_MAP_TAB_Y,
        MINI_MAP_TAB_X + MINI_MAP_TAB_WIDTH,
        MINI_MAP_TAB_Y + MINI_MAP_TAB_HEIGHT
    );

    // draw the white outline
    glColor3f(1, 1, 1);
    glLineWidth(2);
    glBegin(GL_LINE_LOOP);
    glVertex2i(MINI_MAP_TAB_X, MINI_MAP_TAB_Y);
    glVertex2i(MINI_MAP_TAB_X + MINI_MAP_TAB_WIDTH, MINI_MAP_TAB_Y);
    glVertex2i(MINI_MAP_TAB_X + MINI_MAP_TAB_WIDTH, MINI_MAP_TAB_Y + MINI_MAP_TAB_HEIGHT);
    glVertex2i(MINI_MAP_TAB_X, MINI_MAP_TAB_Y + MINI_MAP_TAB_HEIGHT);
    glEnd();

    // plus/minus icon, also at MINI_MAP_TAB_X/Y
    glColor3f(1, 1, 1);
    glLineWidth(2);
    glBegin(GL_LINES);
    // horizontal
    glVertex2i(MINI_MAP_TAB_X + 8, MINI_MAP_TAB_Y + MINI_MAP_TAB_HEIGHT / 2);
    glVertex2i(MINI_MAP_TAB_X + MINI_MAP_TAB_WIDTH - 8, MINI_MAP_TAB_Y + MINI_MAP_TAB_HEIGHT / 2);
    // vertical (only when closed)
    if (!miniMapOpen) {
        glVertex2i(MINI_MAP_TAB_X + MINI_MAP_TAB_WIDTH / 2, MINI_MAP_TAB_Y + 8);
        glVertex2i(MINI_MAP_TAB_X + MINI_MAP_TAB_WIDTH / 2, MINI_MAP_TAB_Y + MINI_MAP_TAB_HEIGHT - 8);
    }
    glEnd();

    // Only draw minimap if open
    if (miniMapOpen) {
        renderMiniMap(getAllCelestialBodies());
    }

    // Restore previous state.
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);

    glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

// Render text function.
void renderText(float x, float y, const std::string& text, void* font) {
    glRasterPos2f(x, y);
    for (char c : text) {
        glutBitmapCharacter(font, c);
    }
}

// Unproject mouse coordinates to a ray.
void unProject(int mouseX, int mouseY, float* rayOrigin, float* rayDir) {
    int viewport[4];
    double modelview[16];
    double projection[16];
    float winX, winY;
    double posX, posY, posZ;

    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);

    winX = static_cast<float>(mouseX);
    winY = static_cast<float>(viewport[3] - mouseY);

    // Get point on the near plane (winZ = 0).
    gluUnProject(winX, winY, 0.0, modelview, projection, viewport, &posX, &posY, &posZ);
    rayOrigin[0] = static_cast<float>(posX);
    rayOrigin[1] = static_cast<float>(posY);
    rayOrigin[2] = static_cast<float>(posZ);

    // Get point on the far plane (winZ = 1).
    gluUnProject(winX, winY, 1.0, modelview, projection, viewport, &posX, &posY, &posZ);
    rayDir[0] = static_cast<float>(posX) - rayOrigin[0];
    rayDir[1] = static_cast<float>(posY) - rayOrigin[1];
    rayDir[2] = static_cast<float>(posZ) - rayOrigin[2];

    // Normalize the ray direction.
    float length = sqrt(rayDir[0] * rayDir[0] + rayDir[1] * rayDir[1] + rayDir[2] * rayDir[2]);
    rayDir[0] /= length;
    rayDir[1] /= length;
    rayDir[2] /= length;
}

// Ease-in-out interpolation function.
float easeInOut(float t) {
    return t * t * (3.0f - 2.0f * t);
}

// Setup camera based on current state (locked, transitioning, or free).
void setupCamera() {
    // Camera setup.
    if (cameraLocked && lockedObject) {
        // Calculate camera position relative to the locked object.
        float radYaw = relativeYaw * M_PI / 180.0f;
        float radPitch = relativePitch * M_PI / 180.0f;

        // Spherical to Cartesian conversion.
        float offsetX = relativeDistance * cosf(radPitch) * sinf(radYaw);
        float offsetY = relativeDistance * sinf(radPitch);
        float offsetZ = relativeDistance * cosf(radPitch) * cosf(radYaw);

        // Desired camera position.
        float desiredCamX = lockedObject->getPositionX() + offsetX;
        float desiredCamY = lockedObject->getPositionY() + offsetY;
        float desiredCamZ = lockedObject->getPositionZ() + offsetZ;

        // Smooth transition to the desired camera position using LERP.
        camX += (desiredCamX - camX) * transitionSpeed;
        camY += (desiredCamY - camY) * transitionSpeed;
        camZ += (desiredCamZ - camZ) * transitionSpeed;

        // Update the camera orientation to look at the locked object.
        gluLookAt(
            camX, camY, camZ,                                     // Camera position.
            lockedObject->getPositionX(), lockedObject->getPositionY(), lockedObject->getPositionZ(), // Look at.
            0.0f, 1.0f, 0.0f                                      // Up vector.
        );

        // Check if camera has reached the desired position.
        float dx = desiredCamX - camX;
        float dy = desiredCamY - camY;
        float dz = desiredCamZ - camZ;
        float distance = sqrt(dx * dx + dy * dy + dz * dz);
        if (distance < 0.1f) { // Threshold.
            isTransitioning = false;
            camX = desiredCamX;
            camY = desiredCamY;
            camZ = desiredCamZ;
        }
    }
    else if (isTransitioning) {
        // Calculate interpolation factor.
        static float transitionProgress = 0.0f;
        float targetYaw = defaultYaw;
        float targetPitch = defaultPitch;
        transitionProgress += transitionSpeed;
        if (transitionProgress > 1.0f) transitionProgress = 1.0f;

        // Apply easing function.
        float t = easeInOut(transitionProgress);

        // Interpolate camera position.
        camX = camX * (1.0f - t) + targetCamX * t;
        camY = camY * (1.0f - t) + targetCamY * t;
        camZ = camZ * (1.0f - t) + targetCamZ * t;

        // Interpolate camera orientation for yaw and pitch.
        camYaw = camYaw * (1.0f - t) + targetYaw * t;
        camPitch = camPitch * (1.0f - t) + targetPitch * t;

        // Update the camera orientation to look at the default point (Sun).
        glLoadIdentity();
        gluLookAt(
            camX, camY, camZ,                          // Camera position.
            defaultCamX, defaultCamY, defaultCamZ,     // Look at default point (Sun).
            0.0f, 1.0f, 0.0f                           // Up vector.
        );

        // Check if transition is complete.
        if (transitionProgress >= 1.0f) {
            isTransitioning = false;
            transitionProgress = 0.0f;
        }
    }
    else {
        // Free camera setup.
        glLoadIdentity();
        gluLookAt(
            camX, camY, camZ,                                  // Camera position.
            camX + sin(camYaw * M_PI / 180.0f),                // Look at point (forward vector X).
            camY + sin(camPitch * M_PI / 180.0f),              // Look at point (forward vector Y).
            camZ - cos(camYaw * M_PI / 180.0f),                // Look at point (forward vector Z).
            0.0f, 1.0f, 0.0f                                   // Up vector.
        );
    }
}

void renderMiniMap(const std::vector<CelestialBody*>& bodies)
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, SCREEN_WIDTH, SCREEN_HEIGHT, 0, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // Draw mini-map background.
    glColor3f(0.1f, 0.1f, 0.1f);
    glBegin(GL_QUADS);
    glVertex2i(MINI_MAP_X, MINI_MAP_Y);
    glVertex2i(MINI_MAP_X + MINI_MAP_WIDTH, MINI_MAP_Y);
    glVertex2i(MINI_MAP_X + MINI_MAP_WIDTH, MINI_MAP_Y + MINI_MAP_HEIGHT);
    glVertex2i(MINI_MAP_X, MINI_MAP_Y + MINI_MAP_HEIGHT);
    glEnd();

    // draw white border 
    glColor3f(1.0f, 1.0f, 1.0f);
    glLineWidth(2.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2i(MINI_MAP_X, MINI_MAP_Y);
    glVertex2i(MINI_MAP_X + MINI_MAP_WIDTH, MINI_MAP_Y);
    glVertex2i(MINI_MAP_X + MINI_MAP_WIDTH, MINI_MAP_Y + MINI_MAP_HEIGHT);
    glVertex2i(MINI_MAP_X, MINI_MAP_Y + MINI_MAP_HEIGHT);
    glEnd();

    // Find the sun.
    float sunX = 0.f, sunY = 0.f, sunZ = 0.f;
    for (auto* b : bodies) {
        if (b->getType() == CelestialBodyType::SUN) {
            sunX = b->getPositionX();
            sunY = b->getPositionY();
            sunZ = b->getPositionZ();
            break;
        }
    }
    float centerX = MINI_MAP_X + MINI_MAP_WIDTH * 0.5f;
    float centerY = MINI_MAP_Y + MINI_MAP_HEIGHT * 0.5f;

    // Compute maximum distance from the sun.
    float maxDist = computeMaxDistance(bodies, sunX, sunZ);
    float dynamicScale = (maxDist > 0)
        ? ((MINI_MAP_WIDTH * 0.45f) / maxDist)
        : (1 / 1000);

    glEnable(GL_SCISSOR_TEST);

    // glScissor() expects window-coords with (0,0) at the *bottom-left* of
    // the window, so we have to convert our top-left-based MINI_MAP_Y.
    GLint scissorX = MINI_MAP_X;
    GLint scissorY = SCREEN_HEIGHT - (MINI_MAP_Y + MINI_MAP_HEIGHT);
    glScissor(scissorX, scissorY, MINI_MAP_WIDTH, MINI_MAP_HEIGHT);

    /*  ↓ everything between glEnable/glDisable will now be clipped ↓  */
    // ——— draw predicted orbits ———
    glColor3f(0.5f, 0.5f, 0.5f);
    glLineWidth(1.0f);
    const int ORBIT_SEGMENTS = 60;
    for (auto* body : bodies) {
        if (body->getType() == CelestialBodyType::SUN) continue;
        float dx = body->getPositionX() - sunX;
        float dz = body->getPositionZ() - sunZ;
        float orbitRadius = sqrtf(dx * dx + dz * dz) * dynamicScale;

        glBegin(GL_LINE_LOOP);
        for (int i = 0; i < ORBIT_SEGMENTS; ++i) {
            float theta = 2.0f * M_PI * i / ORBIT_SEGMENTS;
            glVertex2f(centerX + orbitRadius * cosf(theta),
                centerY + orbitRadius * sinf(theta));
        }
        glEnd();
    }
    glDisable(GL_SCISSOR_TEST);
    
    // Draw each body.
    for (auto* body : bodies) {
        float wx = body->getPositionX() - sunX;
        float wz = body->getPositionZ() - sunZ;
        float mapX = centerX + (wx * dynamicScale);
        float mapY = centerY + (wz * dynamicScale);
        if (body->getType() == CelestialBodyType::SUN)
            glColor3f(1.f, 1.f, 0.f);
        else if (body->getType() == CelestialBodyType::PLANET) {
            if (dynamic_cast<GoldilocksPlanet*>(body))
                glColor3f(0.f, 1.f, 0.f);
            else if (dynamic_cast<GasGiantPlanet*>(body))
                glColor3f(1.f, 0.5f, 0.f);
            else
                glColor3f(0.3f, 0.5f, 1.f);
        }
        else if (body->getType() == CelestialBodyType::MOON)
            glColor3f(0.7f, 0.7f, 0.7f);
        else if (body->getType() == CelestialBodyType::ASTEROID)
            glColor3f(0.6f, 0.6f, 0.6f);
        else
            glColor3f(1.f, 1.f, 1.f);
        float r = 3.0f;
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(mapX, mapY);
        for (int i = 0; i <= 12; i++) {
            float angle = i * (2.0f * M_PI / 12);
            glVertex2f(mapX + r * cosf(angle), mapY + r * sinf(angle));
        }
        glEnd();
    }

    // cyan arrow for camera position
    {
        // compute camera offset in world coords
        float wx = camX - sunX;
        float wz = camZ - sunZ;
        // map to minimap
        float mapX = centerX + wx * dynamicScale;
        float mapY = centerY + wz * dynamicScale;

        // size of arrow
        const float arrowSize = 8.0f;
        // draw a filled triangle pointing "up" (positive Z)
        glColor3f(0.0f, 1.0f, 1.0f); // cyan
        glBegin(GL_TRIANGLES);
        // tip
        glVertex2f(mapX, mapY - arrowSize);
        // left
        glVertex2f(mapX - arrowSize * 0.6f, mapY + arrowSize * 0.6f);
        // right
        glVertex2f(mapX + arrowSize * 0.6f, mapY + arrowSize * 0.6f);
        glEnd();

        // optional: outline
        glColor3f(0, 0, 0);
        glLineWidth(1.0f);
        glBegin(GL_LINE_LOOP);
        glVertex2f(mapX, mapY - arrowSize);
        glVertex2f(mapX - arrowSize * 0.6f, mapY + arrowSize * 0.6f);
        glVertex2f(mapX + arrowSize * 0.6f, mapY + arrowSize * 0.6f);
        glEnd();
    }
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}