#include "Erosion.hpp"
#include "simplex/SimplexNoise.hpp"
#include <tiffio.h>
#include <iostream>

void writeImage(const char* name, int size, uint16_t* buffer, int sizeOfBuffer) {
    TIFF* tif = TIFFOpen(name, "w");
    if (tif) {
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, size);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, size);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 1);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, size);
        TIFFSetField(tif, TIFFTAG_ORIENTATION, (int)ORIENTATION_TOPLEFT);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        TIFFWriteEncodedStrip(tif, 0, buffer, sizeOfBuffer);
        TIFFWriteDirectory(tif);
        TIFFClose(tif);
    }
}

std::vector<float> generateMap(int resolution) {
    std::vector<float> buf(resolution * resolution);
    const SimplexNoise noise(1.0f, 0.5f, 1.99f, 0.5f);
    for (int i = 0; i < resolution * resolution; i++) {
        //std::cout << ((noise.fractal(8, i / resolution, i % resolution) + 1) / 2 * __UINT16_MAX__) << std::endl;
        buf[i] = ((noise.fractal(8, (float)i / resolution / resolution, (i % resolution) / (float)resolution) + 1) / 2);
    }
    return buf;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Usage is ./Terrain-Generator [filename] [resolution] [iterations]" << std::endl;
        return EXIT_FAILURE;
    }

    int resolution = atoi(argv[2]);

    auto map = generateMap(resolution);
    std::cout << "Finished generating map" << std::endl;

    //std::cout << " " << map[0] << " " << map[1] << " " << map[2] << std::endl;

    Erosion eroder = Erosion();
    eroder.seed = 1231204;
    eroder.erode(&map, resolution, atoi(argv[3]), true);

    // libtiff needs it to be in uint16_t since we're saving in 16 bits
    std::vector<uint16_t> toSave(map.size());
    for (int s = 0; s < map.size(); s++) {
        toSave[s] = map.at(s) * __UINT16_MAX__;
    }
    
    writeImage(argv[1], resolution, &toSave[0], 2 * resolution * resolution);
    
    return EXIT_SUCCESS;
}
