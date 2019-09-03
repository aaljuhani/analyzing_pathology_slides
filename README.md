# analyzing_pathology_slides
analyzing_pathology_slides using HistomicsTK api

## Conda Virtual Env
1) Run conda create -n tall_cell python=3.6
2) source activate tall_cell


## Requirements
- pip install tqdm
- pip install numpy
- pip install HistomicsTK
- conda install -c bioconda openslide
- conda install libiconv

## Run Code on OSC
- cd /users/PAS0272/osu10259/ondemand/tall_cell/analyzing_pathology_slides
- module load python/3.6
- source activate tall_cell

### Analyse all WSI in a folder
```python slide.py --analyse_all "Path to WSIs folder"```

### Analyse single WSI
```python slide.py --analyse_WSI "Path to WSI file"```

### Extract single tile
```python slide.py --get_single_tile "Path to WSI file" xcoor_ycoor_magnification```

### Analyse single tile
```python slide.py --analyse_single_tile "Path to WSI file" xcoor_ycoor_magnification```










