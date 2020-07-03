#include "Erosion.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <math.h>
#include <fstream>
#include <CL/cl.hpp>
#include <omp.h>
#include <stdexcept>

template <compute_mode_e COMPUTE_MODE>
Erosion<COMPUTE_MODE>::Erosion(unsigned mapSize) :
erosionBrush(std::vector<float>()),
mapSize(mapSize)
{
    initializeBrushIndices();
}

#define MEM_SIZE 14
template <>
Erosion<compute_mode_e::gpu>::Erosion(unsigned mapSize) :
mapSize(mapSize)
{
    initializeBrushIndices();

    cl_int errcode_ret;
    // Get platform
    cl_platform_id platform_id;
    cl_uint ret_num_platforms;
    clGetPlatformIDs(1, &platform_id, &ret_num_platforms);

    // Get device
    cl_device_id device_id;
    cl_uint ret_num_devices;
    clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);

    char device_string[255];
    clGetDeviceInfo(device_id, CL_DEVICE_NAME, sizeof(device_string), &device_string, NULL);
    std::cout << "Using GPU device: " << device_string << std::endl;

    // Create context based on device
    context = clCreateContext(NULL, 1, &device_id, NULL, NULL, NULL);
    // Create command queue based on context
    queue = clCreateCommandQueue(context, device_id, 0, NULL);

    // Allocate buffers
    map_device = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * mapSize * mapSize, NULL, &errcode_ret);
    if (errcode_ret != CL_SUCCESS)
        throw std::runtime_error("Insufficient device memory (map_device)!");

    erosionBrush_device = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * erosionBrush.size(), NULL, &errcode_ret);
    if (errcode_ret != CL_SUCCESS)
        throw std::runtime_error("Insufficient device memory (erosionBrush_device)!");

    // Load and build kernel
    std::ifstream file("../src/Erosion.cl");

    if (!file.is_open())
        throw std::runtime_error("Error opening Erosion.cl!");

    std::string erosion_cl;
    char t;
    while (file.get(t))
        erosion_cl += t;

    const char *erosion_cl_c_str = erosion_cl.c_str();
    const size_t erosion_cl_length = erosion_cl.length();

    program = clCreateProgramWithSource(context, 1, &erosion_cl_c_str, &erosion_cl_length, NULL);
    clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
    kernel = clCreateKernel(program, "erode", &errcode_ret);
    if (errcode_ret != CL_SUCCESS)
        throw std::runtime_error("Error loading or building Erosion.cl!");
}

template <compute_mode_e COMPUTE_MODE>
Erosion<COMPUTE_MODE>::~Erosion() {
}

template <>
Erosion<compute_mode_e::gpu>::~Erosion() {
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseMemObject(map_device);
    clReleaseMemObject(erosionBrush_device);
    clReleaseCommandQueue(queue);
    clReleaseContext(context);
}

template <compute_mode_e COMPUTE_MODE>
void Erosion<COMPUTE_MODE>::setSeed(int newSeed) {
    seed = newSeed;
    Random::seed(seed);
}

template <compute_mode_e COMPUTE_MODE>
void Erosion<COMPUTE_MODE>::erode(std::vector<float> &map, unsigned numIterations) {
    #pragma omp parallel for schedule(static) if (COMPUTE_MODE == compute_mode_e::parallel)
    for (unsigned iteration = 0; iteration < numIterations; iteration++) {
        // Creates the droplet at a random X and Y on the map
        struct point_float pos = {.x = Random::get<float>(0, mapSize - 1), .y = Random::get<float>(0, mapSize - 1)};
        struct point_float dir = {.x = 0, .y = 0};
        float speed = initialSpeed;
        float water = initialWaterVolume;
        float sediment = 0;

        // Simulates the droplet only up to it's max lifetime, prevents an infite loop
        for (int lifetime = 0; lifetime < maxDropletLifetime; lifetime++) {
            struct point_int node = {.x = (int)pos.x, .y = (int)pos.y};

            // Calculates the droplet offset inside the cell
            struct point_float cellOffset = {.x = pos.x - node.x, .y = pos.y - node.y};

            // Calculate droplet's height and direction of flow with bilinear interpolation of surrounding heights
            HeightAndGradient heightAndGradient = calculateHeightAndGradient(map, pos);

            // Update the droplet's direction and position (move position 1 unit regardless of speed)
            dir.x = (dir.x * inertia - heightAndGradient.gradient.x * (1 - inertia));
            dir.y = (dir.y * inertia - heightAndGradient.gradient.y * (1 - inertia));

            // Normalize direction
            float len = std::sqrt(dir.x * dir.x + dir.y * dir.y);
            if (len != 0) {
                dir.x /= len;
                dir.y /= len;
            }

            pos.x += dir.x;
            pos.y += dir.y;

            // Stop simulating droplet if it's not moving or has flowed over edge of map
            if ((dir.x == 0 && dir.y == 0) || pos.x < 0 || pos.x >= mapSize - 1 || pos.y < 0 || pos.y >= mapSize - 1)
                break;

            // Find the droplet's new height and calculate the deltaHeight
            float deltaHeight = calculateHeightAndGradient(map, pos).height - heightAndGradient.height;

            // Calculate the droplet's sediment capacity (higher when moving fast down a slope and contains lots of water)
            float sedimentCapacity = std::max(-deltaHeight * speed * water * sedimentCapacityFactor, minSedimentCapacity);

            // If carrying more sediment than capacity, or if flowing uphill:
            if (sediment > sedimentCapacity || deltaHeight > 0) {
                int dropletIndex = node.y * mapSize + node.x;

                // If moving uphill (deltaHeight > 0) try fill up to the current height, otherwise deposit a fraction of the excess sediment
                float amountToDeposit = (deltaHeight > 0) ? std::min (deltaHeight, sediment) : (sediment - sedimentCapacity) * depositSpeed;
                sediment -= amountToDeposit;

                // Add the sediment to the four nodes of the current cell using bilinear interpolation
                // Deposition is not distributed over a radius (like erosion) so that it can fill small pits
                map[dropletIndex] += amountToDeposit * (1 - cellOffset.x) * (1 - cellOffset.y);
                map[dropletIndex + 1] += amountToDeposit * cellOffset.x * (1 - cellOffset.y);
                map[dropletIndex + mapSize] += amountToDeposit * (1 - cellOffset.x) * cellOffset.y;
                map[dropletIndex + mapSize + 1] += amountToDeposit * cellOffset.x * cellOffset.y;
            } else {
                // Erode a fraction of the droplet's current carry capacity.
                // Clamp the erosion to the change in height so that it doesn't dig a hole in the terrain behind the droplet
                float amountToErode = std::min((sedimentCapacity - sediment) * erodeSpeed, -deltaHeight);

                // Use erosion brush to erode from all nodes inside the droplet's erosion radius
                unsigned brushIndex = 0;
                for (int y = -erosionRadius; y <= (int)erosionRadius; y++) {
                    for (int x = -erosionRadius; x <= (int)erosionRadius; x++) {
                        int dropletIndex = (node.y + y) * mapSize + (node.x + x);
                        if (dropletIndex < 0 || dropletIndex > (int)mapSize) continue;
                        float weighedErodeAmount = amountToErode * erosionBrush[brushIndex];
                        float deltaSediment = (map[dropletIndex] < weighedErodeAmount) ? map[dropletIndex] : weighedErodeAmount;
                        map[dropletIndex + brushIndex] -= deltaSediment;
                        sediment += deltaSediment;
                        brushIndex++;
                    }
                }
            }

            speed = std::sqrt(speed * speed + std::abs(deltaHeight) * gravity);
            water *= (1 - evaporateSpeed);
        }
    }
}

template <>
void Erosion<compute_mode_e::gpu>::erode(std::vector<float> &map, unsigned numIterations) {
    cl_int errcode_ret;
    rand_device = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(struct point_float) * numIterations, NULL, &errcode_ret);
    if (errcode_ret != CL_SUCCESS)
        throw std::runtime_error("Insufficient device memory (rand_device)!");

    unsigned erosionBrushSize = erosionBrush.size();

    std::vector<struct point_float> rand;
    for(unsigned i = 0; i < numIterations; i++)
        rand.push_back((struct point_float){.x = Random::get<float>(0, mapSize - 1), .y = Random::get<float>(0, mapSize - 1)});

    // Copy from host to device
    clEnqueueWriteBuffer(queue, map_device, CL_TRUE, 0, sizeof(float) * map.size(), map.data(), 0, NULL, NULL);
    clEnqueueWriteBuffer(queue, erosionBrush_device, CL_TRUE, 0, sizeof(float) * erosionBrush.size(), erosionBrush.data(), 0, NULL, NULL);
    clEnqueueWriteBuffer(queue, rand_device, CL_TRUE, 0, sizeof(struct point_float) * numIterations, rand.data(), 0, NULL, NULL);

    // Set parameters
    clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&map_device);
    clSetKernelArg(kernel, 1, sizeof(unsigned), (void *)&mapSize);
    clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&erosionBrush_device);
    clSetKernelArg(kernel, 3, sizeof(unsigned), (void *)&erosionBrushSize);
    clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&rand_device);
    clSetKernelArg(kernel, 5, sizeof(unsigned), (void *)&numIterations);

    // Launch program
    clEnqueueTask(queue, kernel, 0, NULL, NULL);
    clFlush(queue);
    clFinish(queue);

    // Copy from device to host
    clEnqueueReadBuffer(queue, map_device, CL_TRUE, 0, sizeof(float) * map.size(), map.data(), 0, NULL, NULL);

    clReleaseMemObject(rand_device);
}

template <compute_mode_e COMPUTE_MODE>
struct HeightAndGradient Erosion<COMPUTE_MODE>::calculateHeightAndGradient(std::vector<float> &nodes, struct point_float pos) {
    struct point_int coord = {.x = (int)pos.x, .y = (int)pos.y};

    // Calculate droplet's offset inside the cell
    struct point_float offset = {.x = pos.x - coord.x, .y = pos.y - coord.y};

    // Calculate heights of the nodes
    int nodeIndexNW = coord.y * mapSize + coord.x;
    float heightNW = nodes[nodeIndexNW];
    float heightNE = nodes[nodeIndexNW + 1];
    float heightSW = nodes[nodeIndexNW + mapSize];
    float heightSE = nodes[nodeIndexNW + mapSize + 1];

    // Calculate droplet's direction of flow with bilinear interpolation of height difference along the edges
    struct point_float gradient = {.x = (heightNE - heightNW) * (1 - offset.y) + (heightSE - heightSW) * offset.y,
                                   .y = (heightSW - heightNW) * (1 - offset.x) + (heightSE - heightNE) * offset.x};

    // Calculate height with bilinear interpolation of the heights of the nodes of the cell
    float height = heightNW * (1 - offset.x) * (1 - offset.y) + heightNE * offset.x * (1 - offset.y) + heightSW * (1 - offset.x) * offset.y + heightSE * offset.x * offset.y;

    return (HeightAndGradient){.height = height, .gradient = gradient};
}

template <compute_mode_e COMPUTE_MODE>
void Erosion<COMPUTE_MODE>::initializeBrushIndices() {
    // Calculate erosion brush weights
    for (int y = -erosionRadius; y <= (int)erosionRadius; y++)
        for (int x = -erosionRadius; x <= (int)erosionRadius; x++)
            erosionBrush.push_back( erosionRadius - std::sqrt(x * x + y * y) );
}

// Guarantees that all class templates are compiled
template class Erosion<compute_mode_e::cpu>;
template class Erosion<compute_mode_e::parallel>;
template class Erosion<compute_mode_e::gpu>;
