# Particle Wave Runner

This project is a playground for exploring scientific computing through the
solution of numerical PDEs using the material point method (MPM). The code is
developed from the ground up with parallel execution and scalability in mind.

The ultimate demonstrator problem will focus on elastic wave propagation
(hence "Particle Wave Runner" or "pwr") using a vareity of approaches: CutMesh
MPM, High-order MPM, ect.

Use at your own risk.


## Dependencies

Compiling on WSL (Ubuntu)

```
sudo apt update
sudo apt upgrade
sudo apt install -y
```


Install `clang-format`, `cmake`, and `mpi`

```
sudo apt install clang-format
sudo apt install cmake
sudo apt install libopenmpi-dev
```


Install `vtk`

```
mkdir ~/vtk

git clone https://gitlab.kitware.com/vtk/vtk.git ~/vtk/source

cmake -S ~/vtk/source -B ~/vtk/build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DVTK_BUILD_TESTING=OFF \
  -DVTK_BUILD_EXAMPLES=OFF \
  -DVTK_BUILD_DOCUMENTATION=OFF \
  -DVTK_GROUP_ENABLE_Rendering=NO \
  -DVTK_GROUP_ENABLE_Imaging=NO \
  -DVTK_GROUP_ENABLE_Views=NO \
  -DVTK_GROUP_ENABLE_Qt=NO \
  -DVTK_GROUP_ENABLE_Web=NO \
  -DVTK_GROUP_ENABLE_MPI=NO

cmake --build ~/vtk/build -j 1

cmake --install ~/vtk/build --prefix ~/vtk/install
```

> `-j 1` is very slow, but stable for WSL. Ignore compiler warnings.


Install `Catch2`

Save preferred release into `~/Catch-X.XX.X/` directory.

From `~/Catch2-X.XX.X/` directory, run

```
cmake -B build

cmake --build build

cmake --install build --prefix install
```

> CMakeLists.txt currently expects Version 3.13.0; adjust as needed.


## Compile

From root directory, run

```
cmake -S . -B build
cmake --build build
```


## Run

From `build/` directory, run

```
mpirun -n <N> ./particle_wave_runner
```

where `<N>` is the number of MPI processes


## Format

From root directory, run

```
sh format.sh
```


## Visualization

Install `paraview`

```
sudo apt install paraview
```


## Testing

Testing utilizes Catch2 library.

From `build/tests/` directory, run

```
mpirun -n <N> ./particle_wave_runner_tests
```

where `<N>` is the number of MPI processes

> Only run using up to four MPI processes.
