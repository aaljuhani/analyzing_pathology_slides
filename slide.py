import py_wsi
import py_wsi.imagepy_toolkit as tk
import config as cfg
from tqdm import tqdm
import histomicstk as htk

import histomicstk.preprocessing.color_normalization as htk_cnorm
import histomicstk.segmentation as htk_seg
from PIL import Image

import large_image
import os
import numpy as np
from os import listdir
from os.path import isfile, join
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt




class WholeSlideImage():
    """
    A Class that utilize HistomicsTK to tile SVS Whole Slide Images

    """
    def __init__(self):
        self.filedir = cfg.FILE_DIR
        #self.turtle = py_wsi.Turtle(self.filedir, cfg.DB_LOCATION, cfg.DB_NAME, storage_type='disk')
        #self.files=self.turtle.files
        #self.print_slides_info()

    def print_slides_info(self):
        print("Total WSI images:    " + str(self.turtle.num_files))
        print("LMDB name:           " + str(self.turtle.db_name))
        print("File names:          " + str(self.turtle.files))

    def get_wsi_files(self):
        print("get wsi files", self.filedir)
        self.files = np.array([file for file in listdir(self.filedir)
                         if isfile(join(self.filedir, file)) and cfg.IMG_FORMAT in file])

        return self.files


    def get_slide_metadata(self, file):
        ts = large_image.getTileSource(os.path.join(cfg.FILE_DIR, file))
        print(ts.getMetadata())

        # Get the magnification associated with all levels of the image pyramid
        for i in range(ts.levels):
            print('Level-{} : {}'.format(
                i, ts.getMagnificationForLevel(level=i)))

        # compute tissue/foreground mask at low-res for whole slide images
        #im_fgnd_mask_lres, fgnd_seg_scale = cli_utils.segment_wsi_foreground_at_low_res(ts)

    def retrive_wsi_low_res(self, file):
        print("show wsi", file)
        ts = large_image.getTileSource(os.path.join(cfg.FILE_DIR, file))
        im_low_res, _ = ts.getRegion(
            scale=dict(magnification=1.25),
            format=large_image.tilesource.TILE_FORMAT_NUMPY
        )
        plt.imshow(im_low_res)
        plt.savefig(file+'.png')

    def tile_wsi(self, file , magnification, tile_width, tile_hight, tile_overlap_x, tile_overlap_y):
        ts = large_image.getTileSource(os.path.join(cfg.FILE_DIR, file))

        num_tiles = 0

        tile_means = []
        tile_areas = []

        for tile_info in tqdm(ts.tileIterator(
                #region=dict(left=5000, top=5000, width=20000, height=20000, units='base_pixels'),
                scale=dict(magnification=magnification),
                tile_size=dict(width=tile_width, height=tile_hight),
                tile_overlap=dict(x=tile_overlap_x, y=tile_overlap_y),
                format=large_image.tilesource.TILE_FORMAT_PIL
        )):
            #print(tile_info)


            if num_tiles == 100:
                print('Tile-{} = '.format(num_tiles))
                #display(tile_info)

            #img
            im_tile = np.array(tile_info['tile'])
            """ Check if the img is all white (no content in the tile)"""
            extrema = im_tile.convert("L").getextrema()
            if extrema[0] == extrema[1]:
                """ This image is one solid color so no need to save it"""
                pass
            elif (cfg.SAVE_TILES):
                self.save_tile_to_disk(cfg.OUTPUT_DIR, im_tile, (tile_info['x'], tile_info['y']),file, magnification )
            #mean rgb
            tile_mean_rgb = np.mean(im_tile[:, :, :3], axis=(0, 1))

            tile_means.append(tile_mean_rgb)
            tile_areas.append(tile_info['width'] * tile_info['height'])

            num_tiles += 1


        slide_mean_rgb = np.average(tile_means, axis=0, weights=tile_areas)

        print('Number of tiles = {}'.format(num_tiles))
        print('Slide mean color = {}'.format(slide_mean_rgb))



    def save_tile_to_disk(self, output_loc, tile, coords, file_name , magnification):
        """ Saves numpy tiles to .png files (full resolution).
            Meta data is saved in the file name.
            - output_loc       folder to save images in
            - tile             numpy image
            - coords            x, y tile coordinates
            - file_name         original source WSI name
        """

        # Construct the new PNG filename
        tile_fname = file_name + "_" + str(coords[0]) + "_" + str(coords[1]) + "_" + str(magnification) + "_"

        # Save the image.
        Image.fromarray(tile).save(output_loc + tile_fname + ".png")



if __name__ == '__main__':
    WSI = WholeSlideImage()
    wsi_files = WSI.get_wsi_files()

    for f in tqdm(WSI.files):
        WSI.get_slide_metadata(f)
        WSI.tile_wsi(f, cfg.MAGNIFICATION , 512, 512, 0 , 0)

    '''
    for f in s.files:
        s.load_wsi(f)

    for f in s.files:
        s.retrive_wsi_low_res(f)
    
    for f in s.files:
        s.tile_wsi(f)
    '''

    #s.tile_wsi("file")
