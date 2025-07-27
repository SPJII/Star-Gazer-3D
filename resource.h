// resource.h
// Defines the ResourceType enumeration and ResourceNode struct.

#ifndef RESOURCE_H
#define RESOURCE_H

// Enumeration of different resource types.
enum class ResourceType {
    Coal,
    Wood,
    Iron,
    Copper,
    Tin,
    Gold,
    Zinc,
    Lead,
    // Add more resources as needed.
};

// Struct representing a resource node on the planet.
struct ResourceNode {
    float positionX, positionY, positionZ; // Position on the planet surface.
    ResourceType type;                     // Type of resource.
};

#endif // RESOURCE_H
