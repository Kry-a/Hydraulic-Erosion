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

kernel void erode(global float* map, unsigned mapSize, global float* erosionBrush, unsigned erosionBrushSize, global struct point_float *rand, unsigned numIterations) {
    printf("mapSize %u\n", mapSize);
    printf("erosionBrushSize %u\n", erosionBrushSize);
    printf("numIterations %u\n", numIterations);
    for (long i = 0; i < erosionBrushSize; i++)
        printf("%.2f ", erosionBrush[i]);
}

struct HeightAndGradient calculateHeightAndGradient(float* nodes, unsigned mapSize, struct point_float pos) {
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

    
    return (struct HeightAndGradient){.height = height, .gradient = gradient};
}
