#include "Erosion.hpp"
#include "../simplex/SimplexNoise.hpp"
#include <tiffio.h>
#include <iostream>
#include <omp.h>
#include <string.h>

void writeImage(const char* name, int size, uint16_t* buffer, int sizeOfBuffer) {
    TIFF* tif = TIFFOpen(name, "w");
    if (tif) {
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, size);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, size);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, size);
        TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        TIFFWriteEncodedStrip(tif, 0, buffer, sizeOfBuffer);
        TIFFWriteDirectory(tif);
        TIFFClose(tif);
    } else {
        std::cerr << "Cannot save to file" << std::endl;
    }
}

template <compute_mode_e COMPUTE_MODE>
void generateMap(std::vector<float>& buffer, unsigned resolution) {
    const SimplexNoise noise(1.0f, 0.5f, 1.99f, 0.5f);

    #pragma omp parallel for schedule(static) if (COMPUTE_MODE == compute_mode_e::parallel)
    for (unsigned i = 0; i < resolution * resolution; i++)
        buffer[i] = ((noise.fractal(8, (float)i / resolution / resolution, (i % resolution) / (float)resolution) + 1) / 2);
}

int main(int argc, char* argv[]) {
    if (argc < 5 || argc > 6) {
        std::cerr << "Usage is ./Hydraulic-Erosion mode filename resolution iterations [seed]" << std::endl;
        return EXIT_FAILURE;
    }

    unsigned resolution = atoi(argv[3]);

    std::vector<float> map(resolution * resolution);

    if (!strcmp(argv[1], "serial")) {
        generateMap<compute_mode_e::serial>(map, resolution);

        Erosion<compute_mode_e::serial> eroder(resolution);
        if (argc == 6) eroder.setSeed(atoi(argv[5]));
        eroder.erode(map, atoi(argv[4]));
    } else if (!strcmp(argv[1], "parallel")) {
        generateMap<compute_mode_e::parallel>(map, resolution);

        Erosion<compute_mode_e::parallel> eroder(resolution);
        if (argc == 6) eroder.setSeed(atoi(argv[5]));
        eroder.erode(map, atoi(argv[4]));
    } else {
        std::cerr << "`mode` must be one of serial or parallel." << std::endl;
        return EXIT_FAILURE;
    }

    // libtiff needs it to be in uint16_t since we're saving in 16 bits
    std::vector<uint16_t> toSave(map.size());
    #pragma omp parallel for schedule(static) if (!strcmp(argv[1], "parallel"))
    for (unsigned i = 0; i < map.size(); i++)
        toSave[i] = map[i] * __UINT16_MAX__;

    writeImage(argv[2], resolution, toSave.data(), sizeof(uint16_t) * toSave.size());

    return EXIT_SUCCESS;
}
