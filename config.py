
''' folder directory of SVS images'''
FILE_DIR = "/fs/project/PAS0272/Tall_Cells/"

''' input img format'''
IMG_FORMAT = '.svs'

''' Where to save output tiles'''
OUTPUT_DIR = "/fs/project/PAS0272/Tall_Cells/tiles/"

DB_NAME = "tiles"

TILE_H_W = (1024,1024)

LEVEL = 8

MAGNIFICATION = 20

OVERLAP_X_Y = (0,0)

SAVE_TILES = True

# increase this to reject more
REJECT_THRESHOLD = 200

MAX_WHITE_SIZE = (TILE_H_W[0] * TILE_H_W[1]) / 2

