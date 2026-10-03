## Установка
1. Склонировать репозиторий
```bash
git clone https://github.com/MikhSuv/kepler-obake.git
```
2. Установить [conda](https://www.anaconda.com/download/success?reg=skipped-miniconda)
3. Создать новое окружение
```bash
conda create --name env-name
```
4. Активировать окружение
```bash
conda activate env-name
```
5. Установить библиотеку [obake](https://github.com/bluescarni/obake)
```
conda install conda-forge::obake conda-forge::obake-devel
```
Для работы с obake созданное окружение должно быть активно.
## Сборка

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
cmake --build build -j$(nproc)
```

