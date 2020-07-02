#include "Erosion.hpp"
#include <cmath>
#include <iostream>
#include <fenv.h>

void Erosion::initialize(int mapSize, bool resetSeed) {
    //feenableexcept(FE_INVALID | FE_OVERFLOW);
    if (resetSeed || !hasSeed || currentSeed != seed) {
        Random::seed(seed);
        hasSeed = true;
        currentSeed = seed;
    }

    if (erosionBrushIndices.size() != 0 || currentErosionRadius != erosionRadius || currentMapSize != mapSize) {
        initializeBrushIndices(mapSize, erosionRadius);
        currentErosionRadius = erosionRadius;
        currentMapSize = mapSize;
    }
}

void Erosion::erode(std::vector<float> *map, int mapSize, int numIterations, bool resetSeed) {
    initialize(mapSize, resetSeed);

    for (int iteration = 0; iteration < numIterations; iteration++) {
        // Creates the droplet at a random X and Y on the map
        float posX = Random::get<float>(0, mapSize - 1);
        float posY = Random::get<float>(0, mapSize - 1);
        float dirX = 0;
        float dirY = 0;
        float speed = initialSpeed;
        float water = initialWaterVolume;
        float sediment = 0;

        // Simulates the droplet only up to it's max lifetime, prevents an infite loop
        for (int lifetime = 0; lifetime < maxDropletLifetime; lifetime++) {
            int nodeX = (int)posX;
            int nodeY = (int)posY;
            int dropletIndex = nodeY * mapSize + nodeX;
            // Calculates the droplet offset inside the cell
            float cellOffsetX = posX - nodeX;
            float cellOffsetY = posY - nodeY;

            // Calculate droplet's height and direction of flow with bilinear interpolation of surrounding heights
            HeightAndGradient* heightAndGradient = calculateHeightAndGradient(map, mapSize, posX, posY);
            // Update the droplet's direction and position (move position 1 unit regardless of speed)
            dirX = (dirX * inertia - heightAndGradient->gradientX * (1 - inertia));
            dirY = (dirY * inertia - heightAndGradient->gradientY * (1 - inertia));
            // Normalize direction
            float len = std::sqrt(dirX * dirX + dirY * dirY);
            if (len != 0) {
                dirX /= len;
                dirY /= len;
            }
            posX += dirX;
            posY += dirY;

            // Stop simulating droplet if it's not moving or has flowed over edge of map
            if ((dirX == 0 && dirY == 0) || posX < 0 || posX >= mapSize - 1 || posY < 0 || posY >= mapSize - 1) {
                break;
            }

            // Find the droplet's new height and calculate the deltaHeight
            auto newHeightYes = calculateHeightAndGradient(map, mapSize, posX, posY);
            float newHeight = newHeightYes->height;
            float deltaHeight = newHeight - heightAndGradient->height;

            // Calculate the droplet's sediment capacity (higher when moving fast down a slope and contains lots of water)
            float sedimentCapacity = std::max(-deltaHeight * speed * water * sedimentCapacityFactor, minSedimentCapacity);

            // If carrying more sediment than capacity, or if flowing uphill:
            if (sediment > sedimentCapacity || deltaHeight > 0) {
                // If moving uphill (deltaHeight > 0) try fill up to the current height, otherwise deposit a fraction of the excess sediment
                float amountToDeposit = (deltaHeight > 0) ? std::min (deltaHeight, sediment) : (sediment - sedimentCapacity) * depositSpeed;
                sediment -= amountToDeposit;

                // Add the sediment to the four nodes of the current cell using bilinear interpolation
                // Deposition is not distributed over a radius (like erosion) so that it can fill small pits
                map->at(dropletIndex) += amountToDeposit * (1 - cellOffsetX) * (1 - cellOffsetY);
                map->at(dropletIndex + 1) += amountToDeposit * cellOffsetX * (1 - cellOffsetY);
                map->at(dropletIndex + mapSize) += amountToDeposit * (1 - cellOffsetX) * cellOffsetY;
                map->at(dropletIndex + mapSize + 1) += amountToDeposit * cellOffsetX * cellOffsetY;
            } else {
                // Erode a fraction of the droplet's current carry capacity.
                // Clamp the erosion to the change in height so that it doesn't dig a hole in the terrain behind the droplet
                float amountToErode = std::min((sedimentCapacity - sediment) * erodeSpeed, -deltaHeight);

                // Use erosion brush to erode from all nodes inside the droplet's erosion radius
                for (int brushPointIndex = 0; brushPointIndex < erosionBrushIndices[dropletIndex]->size(); brushPointIndex++) {
                    int nodeIndex = erosionBrushIndices[dropletIndex]->at(brushPointIndex);
                    float weighedErodeAmount = amountToErode * erosionBrushWeights[dropletIndex]->at(brushPointIndex);
                    float deltaSediment = (map->at(nodeIndex) < weighedErodeAmount) ? map->at(nodeIndex) : weighedErodeAmount;
                    map->at(nodeIndex) -= deltaSediment;
                    sediment += deltaSediment;
                }
            }

            speed = std::sqrt(speed * speed + std::abs(deltaHeight) * gravity);
            water *= (1 - evaporateSpeed);

            delete heightAndGradient;
            delete newHeightYes;
        }
    }
}

HeightAndGradient* Erosion::calculateHeightAndGradient(std::vector<float> *nodes, int mapSize, float posX, float posY) {
    int coordX = (int)posX;
    int coordY = (int)posY;

    // Calculate droplet's offset inside the cell
    float x = posX - coordX;
    float y = posY - coordY;

    // Calculate heights of the nodes
    int nodeIndexNW = coordY * mapSize + coordX;
    float heightNW = nodes->at(nodeIndexNW);
    float heightNE = nodes->at(nodeIndexNW + 1);
    float heightSW = nodes->at(nodeIndexNW + mapSize);
    float heightSE = nodes->at(nodeIndexNW + mapSize + 1);

    // Calculate droplet's direction of flow with bilinear interpolation of height difference along the edges
    float gradientX = (heightNE - heightNW) * (1 - y) + (heightSE - heightSW) * y;
    float gradientY = (heightSW - heightNW) * (1 - x) + (heightSE - heightNE) * x;

    // Calculate height with bilinear interpolation of the heights of the nodes of the cell
    float height = heightNW * (1 - x) * (1 - y) + heightNE * x * (1 - y) + heightSW * (1 - x) * y + heightSE * x * y;

    HeightAndGradient* yes = new HeightAndGradient();
    yes->height = height;
    yes->gradientX = gradientX;
    yes->gradientY = gradientY;

    return yes;
}

void Erosion::initializeBrushIndices(int mapSize, int radius) {
    erosionBrushIndices.reserve(mapSize * mapSize);
    erosionBrushWeights.reserve(mapSize * mapSize);

    std::vector<int> xOffsets(radius * radius * 4);
    std::vector<int> yOffsets(radius * radius * 4);
    std::vector<float> weights(radius * radius * 4);
    float weightSum = 0;
    int addIndex = 0;

    for (int i = 0; i < (mapSize * mapSize) - 1; i++) {
        int centerX = i % mapSize;
        int centerY = i / mapSize;

        if (centerY <= radius || centerY >= mapSize - radius || centerX <= radius + 1 || centerX >= mapSize - radius) {
            weightSum = 0;
            addIndex = 0;
            for (int y = -radius; y <= radius; y++) {
                for (int x = -radius; x <= radius; x++) {
                    float sqrDst = x * x + y * y;
                    if (sqrDst < radius * radius) {
                        int coordX = centerX + x;
                        int coordY = centerY + y;

                        if (coordX >= 0 && coordX < mapSize && coordY >= 0 && coordY < mapSize) {
                            float weight = 1 - std::sqrt(sqrDst) / radius;
                            weightSum += weight;
                            weights[addIndex] = weight;
                            xOffsets[addIndex] = x;
                            yOffsets[addIndex] = y;
                            addIndex++;
                        }
                    }
                }
            }
        }

        int numEntries = addIndex;
        erosionBrushIndices.push_back(new std::vector<int>(numEntries));
        erosionBrushWeights.push_back(new std::vector<float>(numEntries));

        for (int j = 0; j < numEntries; j++) {
            erosionBrushIndices[i]->push_back((yOffsets[j] + centerY) * mapSize + xOffsets[j] + centerX);
            erosionBrushWeights[i]->push_back(weights[j] / weightSum);
        }
    }
}
