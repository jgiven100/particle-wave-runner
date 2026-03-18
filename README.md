# Particle Wave Runner

## Dependencies

Compiling on WSL (Ubuntu)

```
sudo apt update
sudo apt upgrade
sudo apt install -y
```

Install `clang-format`

```
sudo apt install clang-format
```

Install `cmake`

```
sudo apt install cmake
```

Install `mpi`

```
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

## Compile

From root, run

```
cmake -S . -B build
cmake --build build
```

## Run

From build, run

```
mpirun -n <N> ./particle_wave_runner
```
where `<N>` is the number of MPI processes

## Format

From root, run

```
sh format.sh
```

## Visualization

Install `paraview` (if needed)

```
sudo apt install paraview
```

