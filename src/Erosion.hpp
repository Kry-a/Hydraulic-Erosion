#pragma once

#include <type_traits>
#include <vector>
#include <string>
#include <CL/cl.h>
#include <effolkronium/random.hpp>

using Random = effolkronium::random_static;

struct point_float {
    float x, y;
};

struct point_int {
    int x, y;
};

struct HeightAndGradient {
    float height;
    struct point_float gradient;
};

enum compute_mode_e {cpu, parallel, gpu};

template <compute_mode_e COMPUTE_MODE = compute_mode_e::cpu>
class Erosion {
public:
    Erosion(unsigned mapSize);
    ~Erosion();
    void setSeed(int seed);
    void erode(std::vector<float> &map, unsigned numIterations = 1);

private:
    std::vector<float> erosionBrush;
    std::string erosion_cl;
    cl_context context;
    cl_command_queue queue;
    cl_mem map_device, erosionBrush_device, rand_device;
    cl_program program;
    cl_kernel kernel;

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
    HeightAndGradient calculateHeightAndGradient(std::vector<float> &nodes, struct point_float pos);
};
