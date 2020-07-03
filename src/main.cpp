#include "Erosion.hpp"
#include "../simplex/SimplexNoise.hpp"
#include <tiffio.h>
#include <iostream>
#include <omp.h>

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

void generateMap(std::vector<float>& buffer, unsigned resolution) {
    const SimplexNoise noise(1.0f, 0.5f, 1.99f, 0.5f);

    #pragma omp parallel for schedule(static)
    for (unsigned i = 0; i < resolution * resolution; i++)
        buffer[i] = ((noise.fractal(8, (float)i / resolution / resolution, (i % resolution) / (float)resolution) + 1) / 2);
}

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 5) {
        std::cerr << "Usage is ./Hydraulic-Erosion filename resolution iterations [seed]" << std::endl;
        return EXIT_FAILURE;
    }

    unsigned resolution = atoi(argv[2]);

    std::vector<float> map(resolution * resolution);
    generateMap(map, resolution);

    Erosion eroder = Erosion(resolution);
    if (argc == 5) eroder.setSeed(atoi(argv[4]));
    eroder.erode(map, atoi(argv[3]));

    // libtiff needs it to be in uint16_t since we're saving in 16 bits
    std::vector<uint16_t> toSave(map.size());
    #pragma omp parallel for schedule(static)
    for (unsigned i = 0; i < map.size(); i++)
        toSave[i] = map[i] * __UINT16_MAX__;

    writeImage(argv[1], resolution, toSave.data(), sizeof(uint16_t) * toSave.size());

    return EXIT_SUCCESS;
}
