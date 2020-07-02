#ifndef EROSION_HPP
#define EROSION_HPP


#include <vector>
#include <effolkronium/random.hpp>

using Random = effolkronium::random_static;

struct HeightAndGradient {
    float height;
    float gradientX;
    float gradientY;
};

class Erosion {
public:
    int seed;
    int erosionRadius = 3;
    float inertia = 0.05f;
    float sedimentCapacityFactor = 4;
    float minSedimentCapacity = 0.01f;
    float erodeSpeed = 0.3f;
    float depositSpeed = 0.3f;
    float evaporateSpeed = 0.01f;
    float gravity = 4;
    float maxDropletLifetime = 30;

    float initialWaterVolume = 1;
    float initialSpeed = 1;

    bool hasSeed = false;

    void erode(std::vector<float> *map, int mapSize, int numIterations = 1, bool resetSeed = false);

private:
    std::vector<std::vector<int> *> erosionBrushIndices = std::vector<std::vector<int> *>();
    std::vector<std::vector<float> *> erosionBrushWeights = std::vector<std::vector<float> *>();

    int currentSeed;
    int currentErosionRadius;
    int currentMapSize;

    void initialize(int mapSize, bool resetSeed);
    HeightAndGradient* calculateHeightAndGradient(std::vector<float> *nodes, int mapSize,
                                                 float posX, float posY);
    void initializeBrushIndices(int mapSize, int radius);
};


#endif
