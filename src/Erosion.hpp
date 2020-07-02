#pragma once

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
    Erosion(unsigned mapSize);
    void setSeed(int seed);
    void erode(std::vector<float> &map, unsigned numIterations = 1);

private:
    std::vector<std::vector<unsigned> *> erosionBrushIndices = std::vector<std::vector<unsigned> *>();
    std::vector<std::vector<float> *> erosionBrushWeights = std::vector<std::vector<float> *>();

    unsigned mapSize;
    int seed = 1231204;
    unsigned erosionRadius = 3;
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

    void initializeBrushIndices();
    HeightAndGradient calculateHeightAndGradient(std::vector<float> &nodes, float posX, float posY);
};
