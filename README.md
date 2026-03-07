# Particle Wave Runner

## Dependencies

Compiling on WSL (Ubuntu)

```
sudo apt update
sudo apt upgrade
sudo apt install -y
```

Install `clang-format` (if needed)

```
sudo apt install clang-format
```

Install `cmake` (if needed)

```
sudo apt install cmake
```

Install `mpi` (if needed)

```
sudo apt install libopenmpi-dev
```

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

