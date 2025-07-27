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
const int SCREEN_HEIGHT = 1030;

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

// Mini-Map location and scale
static const int MINI_MAP_WIDTH = 450;
static const int MINI_MAP_HEIGHT = 450;
static const int MINI_MAP_X = 10; // bottom-left X
static const int MINI_MAP_Y = SCREEN_HEIGHT - MINI_MAP_HEIGHT - 10; // from top if 0=top
static const float MINI_MAP_SCALE = 1.0f / 450.0f;
// adjust so your system fits in that 200x200 box


// UI control for resource nodes and provinces.
int buttonWidth = 150;
int buttonHeight = 30;
int buttonX = SCREEN_WIDTH - buttonWidth - 10; // 10 pixels from the right.
int buttonY = SCREEN_HEIGHT - buttonHeight - 10; // 10 pixels from the bottom.
int provinceButtonX = buttonX - buttonWidth - 10; // Province button to the left of resource button.

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
            // Check resource and province buttons…
            if (mouseX >= buttonX && mouseX <= buttonX + buttonWidth &&
                mouseY >= buttonY && mouseY <= buttonY + buttonHeight) {
                showResourceNodes = !showResourceNodes;
                break;
            }
            if (mouseX >= provinceButtonX && mouseX <= provinceButtonX + buttonWidth &&
                mouseY >= buttonY && mouseY <= buttonY + buttonHeight) {
                showProvinces = !showProvinces;
                break;
            }
            // Check tab.
            if (mouseX >= SCREEN_WIDTH - TAB_WIDTH && mouseY <= TAB_HEIGHT) {
                tabOpen = !tabOpen;
                break;
            }
            // Then check if click is in mini-map.
            bool inMiniMap = (mouseX >= MINI_MAP_X && mouseX <= MINI_MAP_X + MINI_MAP_WIDTH &&
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
                if (picked && sqrtf(closestDistSq) < mapThreshold) {
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
                else {
                    // If no body was found within threshold, just recenter camera
                    cameraLocked = false;
                    lockedObject = nullptr;

                    // Recenter camera to that world location 
                    // so the user can freely explore that area.
                    targetCamX = dx + sunX;
                    targetCamY = camY;  // keep same height
                    targetCamZ = dz + sunZ;
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
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Render the reset button in the top-left corner.
    int resetBtnX = 10; // 10 pixels from the left.
    int resetBtnY = 10; // 10 pixels from the top.
    int resetBtnWidth = 100;
    int resetBtnHeight = 30;

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

    // Button label.
    {
        glColor3f(1.0f, 1.0f, 1.0f); // White color.
        std::string resetText = "Reset Position";
        // Center the text.
        float textWidth = resetText.length() * 8.0f; // Approximate width per character.
        float textX = resetBtnX + (resetBtnWidth - textWidth) / 2.0f;
        float textY = resetBtnY + (resetBtnHeight - 18.0f) / 2.0f + 5.0f; // 18 is approx. height of the font.

        renderText(textX, textY, resetText, GLUT_BITMAP_HELVETICA_18);
    }

    // Render the tab.
    glColor4f(0.2f, 0.2f, 0.2f, 0.8f); // Semi-transparent dark gray.
    glBegin(GL_QUADS);
    glVertex2i(SCREEN_WIDTH - TAB_WIDTH, 0);
    glVertex2i(SCREEN_WIDTH, 0);
    glVertex2i(SCREEN_WIDTH, TAB_HEIGHT);
    glVertex2i(SCREEN_WIDTH - TAB_WIDTH, TAB_HEIGHT);
    glEnd();

    // Render a plus or minus sign on the tab.
    glColor3f(1.0f, 1.0f, 1.0f); // White color.
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    // Horizontal line.
    glVertex2i(SCREEN_WIDTH - TAB_WIDTH + 8, TAB_HEIGHT / 2);
    glVertex2i(SCREEN_WIDTH - 8, TAB_HEIGHT / 2);
    if (!tabOpen) {
        // Vertical line for plus sign.
        glVertex2i(SCREEN_WIDTH - TAB_WIDTH + TAB_WIDTH / 2, 8);
        glVertex2i(SCREEN_WIDTH - TAB_WIDTH + TAB_WIDTH / 2, TAB_HEIGHT - 8);
    }
    glEnd();

    // If tab is open, render the slider and label.
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

    // Render the resource nodes button.
    glColor4f(0.2f, 0.2f, 0.2f, 0.8f); // Semi-transparent dark gray.
    glBegin(GL_QUADS);
    glVertex2i(buttonX, buttonY);
    glVertex2i(buttonX + buttonWidth, buttonY);
    glVertex2i(buttonX + buttonWidth, buttonY + buttonHeight);
    glVertex2i(buttonX, buttonY + buttonHeight);
    glEnd();

    // Button border.
    glColor3f(1.0f, 1.0f, 1.0f); // White color.
    glLineWidth(2.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2i(buttonX, buttonY);
    glVertex2i(buttonX + buttonWidth, buttonY);
    glVertex2i(buttonX + buttonWidth, buttonY + buttonHeight);
    glVertex2i(buttonX, buttonY + buttonHeight);
    glEnd();

    // Button label.
    {
        glColor3f(1.0f, 1.0f, 1.0f); // White color.
        std::string buttonText = showResourceNodes ? "Hide Resources" : "Show Resources";
        // Center the text.
        float textWidth = buttonText.length() * 8.0f; // Approximate width per character.
        float textX = buttonX + (buttonWidth - textWidth) / 2.0f;
        float textY = buttonY + (buttonHeight - 18.0f) / 2.0f + 5.0f; // 18 is approx. height of the font.

        renderText(textX, textY, buttonText, GLUT_BITMAP_HELVETICA_18);
    }

    // Render the provinces button.
    glColor4f(0.2f, 0.2f, 0.2f, 0.8f); // Semi-transparent dark gray.
    glBegin(GL_QUADS);
    glVertex2i(provinceButtonX, buttonY);
    glVertex2i(provinceButtonX + buttonWidth, buttonY);
    glVertex2i(provinceButtonX + buttonWidth, buttonY + buttonHeight);
    glVertex2i(provinceButtonX, buttonY + buttonHeight);
    glEnd();

    // Button border.
    glColor3f(1.0f, 1.0f, 1.0f); // White color.
    glLineWidth(2.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2i(provinceButtonX, buttonY);
    glVertex2i(provinceButtonX + buttonWidth, buttonY);
    glVertex2i(provinceButtonX + buttonWidth, buttonY + buttonHeight);
    glVertex2i(provinceButtonX, buttonY + buttonHeight);
    glEnd();

    // Button label.
    {
        glColor3f(1.0f, 1.0f, 1.0f); // White color.
        std::string buttonText = showProvinces ? "Hide Provinces" : "Show Provinces";
        // Center the text.
        float textWidth = buttonText.length() * 8.0f; // Approximate width per character.
        float textX = provinceButtonX + (buttonWidth - textWidth) / 2.0f;
        float textY = buttonY + (buttonHeight - 18.0f) / 2.0f + 5.0f; // 18 is approx. height of the font.

        renderText(textX, textY, buttonText, GLUT_BITMAP_HELVETICA_18);
    }

    renderMiniMap(getAllCelestialBodies());

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
    // Dynamic scale: we want the farthest object to appear at half the mini-map width.
    float dynamicScale = (maxDist > 0.f) ? ((MINI_MAP_WIDTH / 1.0f) / maxDist) : (1.0f / 1000.0f);

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
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}