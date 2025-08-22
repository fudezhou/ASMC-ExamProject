# Linux/macOS — parallel (default)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=ON
cmake --build build -j
./build/pricer

# macOS — with Homebrew GCC (OpenMP)
CXX=g++-15 cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=ON
cmake --build build -j
./build/pricer

# Linux/macOS — force serial run
OMP_NUM_THREADS=1 ./build/pricer

# Linux/macOS — N threads
OMP_NUM_THREADS=N ./build/pricer

# Windows — parallel (default)
cmake -S . -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=ON
cmake --build build -j
.\build\pricer.exe

# Windows — force serial
set OMP_NUM_THREADS=1
.\build\pricer.exe   # cmd.exe

# or:
$env:OMP_NUM_THREADS=1; .\build\pricer.exe   # PowerShell

# Windows — N threads
set OMP_NUM_THREADS=N
.\build\pricer.exe   # cmd.exe

# or:
$env:OMP_NUM_THREADS=N; .\build\pricer.exe   # PowerShell

# Clean (ALL PLATFORMS)
cmake --build build --target clean_all

