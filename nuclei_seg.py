import histomicstk as htk
import config as cfg

import numpy as np
import scipy as sp

import skimage.io
import skimage.measure
import skimage.color

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches



class NucleiSegmentation():
    """
    A class that utilize HistomicsTK to segment nuclei given a single tile or WSI tiles
    """
    def __init__(self):
        self.filedir = cfg.TILES_DIR
        self.




