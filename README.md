# Hydraulic-Erosion
Hydraulic Erosion simulation example made in C++ parallelized.

# Building
```sh shell-script
git clone --recursive https://github.com/arthursimas1/Hydraulic-Erosion.git
cd Hydraulic-Erosion
mkdir build
cd build
cmake ..
make
```

# Running
In the folder `build`, you can run the following command line to run the program
```
./Hydraulic-Erosion mode filename resolution iterations [seed]
```
where:
- mode: you can choose between running in `serial` or `parallel`
- filename: output filename, e.g. `eroded.tif`
- resolution: image resolution, e.g. `1024`
- iterations: erosion iterations, as known as quantity of droplets, e.g. `1000000`
- seed (optional): optionally you can provide a seed for the pseudo random number generator, e.g. `42`

# Dependencies
You may need to install the following dependencies manually in order to compile and run the program
- C++ compiler (`sudo apt install g++`)
- libtiff 4.1.0 (`sudo apt install libtiff-dev`)
- OpenMP 4.5 (already bundled in `g++`)
