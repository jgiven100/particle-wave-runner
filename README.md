# particle-wave-runner

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
mkdir build && cd build
cmake ..
make
```

## Run

From build, run

```
mpirun -n <N> ./
```

## Format

From root, run

```
sh format.sh
```

