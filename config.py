import numpy as np


""" ************** FOLDERS and PATH Configration ***************"""

''' folder directory of SVS images'''
FILE_DIR = "/Users/juhania/Documents/projects/tall cell/"

#FILE_DIR = "/fs/project/PAS0272/Tall_Cells/"

''' input img format'''
IMG_FORMAT = '.svs'

''' Where to save output tiles'''
#OUTPUT_DIR = "/fs/project/PAS0272/Tall_Cells/tiles/"

OUTPUT_DIR = "/Users/juhania/Documents/projects/tall cell/tiles/"

#TILES_DIR = "/fs/project/PAS0272/Tall_Cells/tiles/"

TILES_DIR = "/Users/juhania/Documents/projects/tall cell/tiles/"

FEATURES_DIR = "/Users/juhania/Documents/projects/tall cell/feat/"



""" ************** TILES Configration ***************"""

TILE_H_W = (1024,1024)

LEVEL = 8

MAGNIFICATION = 20

OVERLAP_X_Y = (0,0)

SAVE_TILES = True

# increase this to reject more
REJECT_THRESHOLD = 200

MAX_WHITE_SIZE = (TILE_H_W[0] * TILE_H_W[1]) / 2

#reference image for normalization
REF_IMG_NORM = ""


""" ************** Nuclei Segmentation Configration ***************"""
# segment foreground
FORGROUND_THRESHOLD = 60
# run adaptive multi-scale LoG filter
MIN_RADIUS = 10
MAX_RADIUS = 15

# detect and segment nuclei using local maximum clustering
LOCAL_MAX_SEARCH_RADIUS = 10

# filter out small objects
MIN_NUCLEAUS_AREA = 80

""" ************** Stains Configration ***************"""


# create stain to color map
stainColorMap = {
    'hematoxylin': [0.65, 0.70, 0.29],
    'eosin': [0.07, 0.99, 0.11],
    'dab': [0.27, 0.57, 0.78],
    'null': [0.0, 0.0, 0.0]
}

# specify stains of input image
stain_1 = 'hematoxylin'  # nuclei stain
stain_2 = 'eosin'  # cytoplasm stain
stain_3 = 'null'  # set to null of input contains only two stains

# create stain matrix
W = np.array([stainColorMap[stain_1],
              stainColorMap[stain_2],
              stainColorMap[stain_3]]).T



