## Сборка

cmake -S . -B build -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
cmake --build build -j$(nproc)

