#include "Erosion.hpp"
#include <cmath>
#include <iostream>
#include <fenv.h>

Erosion::Erosion(unsigned mapSize) : mapSize(mapSize) {
    //feenableexcept(FE_INVALID | FE_OVERFLOW);
    initializeBrushIndices();
}

void Erosion::setSeed(int newSeed) {
    seed = newSeed;
    Random::seed(seed);
}

void Erosion::erode(std::vector<float> &map, unsigned numIterations) {
    for (unsigned iteration = 0; iteration < numIterations; iteration++) {
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
            HeightAndGradient heightAndGradient = calculateHeightAndGradient(map, posX, posY);

            // Update the droplet's direction and position (move position 1 unit regardless of speed)
            dirX = (dirX * inertia - heightAndGradient.gradientX * (1 - inertia));
            dirY = (dirY * inertia - heightAndGradient.gradientY * (1 - inertia));

            // Normalize direction
            float len = std::sqrt(dirX * dirX + dirY * dirY);
            if (len != 0) {
                dirX /= len;
                dirY /= len;
            }

            posX += dirX;
            posY += dirY;

            // Stop simulating droplet if it's not moving or has flowed over edge of map
            if ((dirX == 0 && dirY == 0) || posX < 0 || posX >= mapSize - 1 || posY < 0 || posY >= mapSize - 1)
                break;

            // Find the droplet's new height and calculate the deltaHeight
            float deltaHeight = calculateHeightAndGradient(map, posX, posY).height - heightAndGradient.height;

            // Calculate the droplet's sediment capacity (higher when moving fast down a slope and contains lots of water)
            float sedimentCapacity = std::max(-deltaHeight * speed * water * sedimentCapacityFactor, minSedimentCapacity);

            // If carrying more sediment than capacity, or if flowing uphill:
            if (sediment > sedimentCapacity || deltaHeight > 0) {
                // If moving uphill (deltaHeight > 0) try fill up to the current height, otherwise deposit a fraction of the excess sediment
                float amountToDeposit = (deltaHeight > 0) ? std::min (deltaHeight, sediment) : (sediment - sedimentCapacity) * depositSpeed;
                sediment -= amountToDeposit;

                // Add the sediment to the four nodes of the current cell using bilinear interpolation
                // Deposition is not distributed over a radius (like erosion) so that it can fill small pits
                map[dropletIndex] += amountToDeposit * (1 - cellOffsetX) * (1 - cellOffsetY);
                map[dropletIndex + 1] += amountToDeposit * cellOffsetX * (1 - cellOffsetY);
                map[dropletIndex + mapSize] += amountToDeposit * (1 - cellOffsetX) * cellOffsetY;
                map[dropletIndex + mapSize + 1] += amountToDeposit * cellOffsetX * cellOffsetY;
            } else {
                // Erode a fraction of the droplet's current carry capacity.
                // Clamp the erosion to the change in height so that it doesn't dig a hole in the terrain behind the droplet
                float amountToErode = std::min((sedimentCapacity - sediment) * erodeSpeed, -deltaHeight);

                // Use erosion brush to erode from all nodes inside the droplet's erosion radius
                for (unsigned brushPointIndex = 0; brushPointIndex < erosionBrushIndices[dropletIndex]->size(); brushPointIndex++) {
                    int nodeIndex = (*erosionBrushIndices[dropletIndex])[brushPointIndex];
                    float weighedErodeAmount = amountToErode * (*erosionBrushWeights[dropletIndex])[brushPointIndex];
                    float deltaSediment = (map[nodeIndex] < weighedErodeAmount) ? map[nodeIndex] : weighedErodeAmount;
                    map[nodeIndex] -= deltaSediment;
                    sediment += deltaSediment;
                }
            }

            speed = std::sqrt(speed * speed + std::abs(deltaHeight) * gravity);
            water *= (1 - evaporateSpeed);
        }
    }
}

HeightAndGradient Erosion::calculateHeightAndGradient(std::vector<float> &nodes, float posX, float posY) {
    int coordX = (int)posX;
    int coordY = (int)posY;

    // Calculate droplet's offset inside the cell
    float x = posX - coordX;
    float y = posY - coordY;

    // Calculate heights of the nodes
    int nodeIndexNW = coordY * mapSize + coordX;
    float heightNW = nodes[nodeIndexNW];
    float heightNE = nodes[nodeIndexNW + 1];
    float heightSW = nodes[nodeIndexNW + mapSize];
    float heightSE = nodes[nodeIndexNW + mapSize + 1];

    // Calculate droplet's direction of flow with bilinear interpolation of height difference along the edges
    float gradientX = (heightNE - heightNW) * (1 - y) + (heightSE - heightSW) * y;
    float gradientY = (heightSW - heightNW) * (1 - x) + (heightSE - heightNE) * x;

    // Calculate height with bilinear interpolation of the heights of the nodes of the cell
    float height = heightNW * (1 - x) * (1 - y) + heightNE * x * (1 - y) + heightSW * (1 - x) * y + heightSE * x * y;

    return (HeightAndGradient){.height = height, .gradientX = gradientX, .gradientY = gradientY};
}

void Erosion::initializeBrushIndices() {
    erosionBrushIndices.reserve(mapSize * mapSize);
    erosionBrushWeights.reserve(mapSize * mapSize);

    std::vector<int> xOffsets(erosionRadius * erosionRadius * 4);
    std::vector<int> yOffsets(erosionRadius * erosionRadius * 4);
    std::vector<float> weights(erosionRadius * erosionRadius * 4);

    float weightSum = 0;
    unsigned addIndex = 0;

    for (unsigned i = 0; i < (mapSize * mapSize) - 1; i++) {
        unsigned centerX = i % mapSize;
        unsigned centerY = i / mapSize;

        if (centerY <= erosionRadius || centerY >= mapSize - erosionRadius || centerX <= erosionRadius + 1 || centerX >= mapSize - erosionRadius) {
            weightSum = 0;
            addIndex = 0;
            for (int y = -erosionRadius; y <= (int)erosionRadius; y++) {
                for (int x = -erosionRadius; x <= (int)erosionRadius; x++) {
                    float sqrDst = x * x + y * y;
                    if (sqrDst < erosionRadius * erosionRadius) {
                        int coordX = centerX + x;
                        int coordY = centerY + y;

                        if (0 <= coordX && coordX < (int)mapSize && 0 <= coordY && coordY < (int)mapSize) {
                            float weight = 1 - std::sqrt(sqrDst) / erosionRadius;
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

        unsigned numEntries = addIndex;
        erosionBrushIndices.push_back(new std::vector<unsigned>(numEntries));
        erosionBrushWeights.push_back(new std::vector<float>(numEntries));

        for (unsigned j = 0; j < numEntries; j++) {
            erosionBrushIndices[i]->push_back((yOffsets[j] + centerY) * mapSize + xOffsets[j] + centerX);
            erosionBrushWeights[i]->push_back(weights[j] / weightSum);
        }
    }
}
