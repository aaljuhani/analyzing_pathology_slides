
import config as cfg
from tqdm import tqdm
import histomicstk as htk
import sys
import histomicstk.preprocessing.color_normalization as htk_cnorm
import histomicstk.segmentation as htk_seg
from PIL import Image
from nuclei_seg import NucleiSegmentation
import cv2
import large_image
import os
import numpy as np
from os import listdir
from os.path import isfile, join
import argparse

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt




class WholeSlideImage():
    """
    An application that utilize HistomicsTK to tile SVS images and analyze them

    """
    def __init__(self):
        self.filedir = cfg.FILE_DIR
        self.nuclei_seg = NucleiSegmentation()

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
        num_empty_tiles = 0
        num_nuclei = 0

        tile_means = []
        tile_areas = []

        for tile_info in tqdm(ts.tileIterator(
                #region=dict(left=5000, top=5000, width=20000, height=20000, units='base_pixels'),
                scale=dict(magnification=magnification),
                tile_size=dict(width=tile_width, height=tile_hight),
                tile_overlap=dict(x=tile_overlap_x, y=tile_overlap_y),
                format=large_image.tilesource.TILE_FORMAT_PIL)):
            #print(tile_info)

            if (num_tiles > 0 and num_tiles % 100 == 0):
                print('Tile - {} '.format(num_tiles))
                print('Empty Tile - {} '.format(num_empty_tiles))
                #display(tile_info)

            #img
            im_tile = np.array(tile_info['tile'])

            # To check the content of the image
            if (self.is_good_tile(im_tile)):
                if(cfg.SAVE_TILES):
                    self.save_tile_to_disk(im_tile, (tile_info['x'], tile_info['y']),file, magnification )
                #mean rgb
                tile_mean_rgb = np.mean(im_tile[:, :, :3], axis=(0, 1))
                tile_means.append(tile_mean_rgb)
                tile_areas.append(tile_info['width'] * tile_info['height'])

                num_tiles += 1

                """ Here is the call for computing morphometry features for each good tile"""
                num_nuclei += self.nuclei_seg.compute_morphometry_feat(im_tile, (tile_info['x'], tile_info['y']),file, magnification)
            else:
                """ This image is one solid color so no need to save it"""
                num_empty_tiles += 1



        slide_mean_rgb = np.average(tile_means, axis=0, weights=tile_areas)
        print('Slide: ', file)
        print('Total num of Nuclei = {}'.format(num_nuclei))
        print('Number of tiles = {}'.format(num_tiles))
        print('Number of empty tiles = {}'.format(num_empty_tiles))
        print('Slide mean color = {}'.format(slide_mean_rgb))

    def is_good_tile(self, tile):
        if(tile.shape[0] < cfg.TILE_H_W[0] or tile.shape[1] < cfg.TILE_H_W[1]):
            return False
        # tile is nparray
        img = cv2.cvtColor(tile, cv2.COLOR_BGR2GRAY)
        blur = cv2.GaussianBlur(img, (5, 5), 0)
        ret3, th3 = cv2.threshold(blur, cfg.REJECT_THRESHOLD, 255, cv2.THRESH_BINARY)
        contours, hierarchy = cv2.findContours(th3, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
        return self.get_cnt_sum(contours, 2) < cfg.MAX_WHITE_SIZE

    def get_cnt_sum(self, contours, topn):
        res = 0
        cnts = sorted(contours, key=lambda x: cv2.contourArea(x))[-topn:]
        return sum([cv2.contourArea(cnt) for cnt in cnts])


    def save_tile_to_disk(self, tile, coords, file_name , magnification):
        """ Saves numpy tiles to .png files (full resolution).
            Meta data is saved in the file name.
            - tile             numpy image
            - coords            x, y tile coordinates
            - file_name         original source WSI name
        """

        # Construct the new PNG filename
        tile_fname = file_name + "_" + str(coords[0]) + "_" + str(coords[1]) + "_" + str(magnification) + "_"

        # Save the image.
        Image.fromarray(tile).save(cfg.OUTPUT_DIR + tile_fname + ".png")

    def get_tile(self, wsi_file, tile_name):
        """
        :param tile_name:
        Should follow the same format as saved tiles:
        tile_info['x']_tile_info['y']_magnification_
        ex: 22528_4096_40_

        :return: save tile in OUTPUT_TILES_FOLDEROUTPUT_DIR

        """
        
        #parse tile name
        tile_name_list = tile_name.split("_")
        print("Tile_name_list", tile_name_list)

        x = int(tile_name_list[0])
        y = int(tile_name_list[1])
        magnification = int(tile_name_list[2])

        print("Get tile: ", wsi_file)

        # Read WSI
        ts = large_image.getTileSource(os.path.join(cfg.FILE_DIR, wsi_file))


        # get region
        im_roi, _ = ts.getRegion(
            region=dict(left=x, top=y, width=cfg.TILE_H_W[0], height=cfg.TILE_H_W[1], units='base_pixels'),
            scale=dict(magnification=magnification),
            format=large_image.tilesource.TILE_FORMAT_NUMPY
        )

        Image.fromarray(im_roi).save("roi_image.png")



if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--analyse_all', help="Provide a folder that has all WSI in svs format")
    parser.add_argument('--analyse_WSI', help="Provide a path for a single WSI")
    parser.add_argument('--get_single_tile', nargs=2,
                        help="Provide WSI path and a tile coordinates and magnification level in this format(ex: 12000_3488_40)")
    parser.add_argument('--analyse_single_tile', nargs=2,
                        help="Provide WSI path and a tile coordinates and magnification level in this format(ex: 12000_3488_40)")

    args = parser.parse_args()

    if args.analyse_all:
        print("analyse all svs files in a folder")
        cfg.FILE_DIR = args.get_single_tile

        # Code for tiling WSI and compute morpho feat
        wsi = WholeSlideImage()
        wsi_files = wsi.get_wsi_files()
        print(wsi_files)

        for f in tqdm(wsi_files):
            wsi.get_slide_metadata(f)
            wsi.tile_wsi(f, cfg.MAGNIFICATION, cfg.TILE_H_W[0], cfg.TILE_H_W[1], cfg.OVERLAP_X_Y[0], cfg.OVERLAP_X_Y[1])


    elif (args.analyse_WSI):
        print("analyse WSI")
        cfg.WSI_PATH = args.analyse_WSI
        wsi = WholeSlideImage()
        wsi.get_slide_metadata(cfg.WSI_PATH)
        wsi.tile_wsi(cfg.WSI_PATH, cfg.MAGNIFICATION, cfg.TILE_H_W[0], cfg.TILE_H_W[1], cfg.OVERLAP_X_Y[0], cfg.OVERLAP_X_Y[1])


    elif args.get_single_tile:
        print("Get Single Tile")
        WSI_PATH = args.get_single_tile[0]
        # TODO : check format
        # should be in the format of xcoor_ycoor_magnification
        tile_xcoor_ycoor_magnification = args.get_single_tile[1]

        wsi = WholeSlideImage()
        wsi.get_tile(WSI_PATH, tile_xcoor_ycoor_magnification)

    elif args.analyse_single_tile:
        # TODO: implement
        print("analyse single tile")
        WSI_PATH = args.get_single_tile = [0]
        tile_xcoor_ycoor_magnification = args.get_single_tile = [1]



