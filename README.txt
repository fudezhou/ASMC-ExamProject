To run in serial: 

g++ -O3 -march=native src/mc_explicit.cpp src/main.cpp -Iinclude -o pricer && ./pricer


TO run in parallel:

g++ -O3 -march=native -DUSE_OMP -fopenmp \
  -I./include src/mc_explicit.cpp src/main.cpp -o pricer

export OMP_NUM_THREADS=8   # or however many cores you want
export OMP_PROC_BIND=close # optional, keeps threads near each other
./pricer

Or just:

OMP_NUM_THREADS=8 ./pricer
