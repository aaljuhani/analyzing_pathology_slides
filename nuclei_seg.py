import histomicstk as htk
import config as cfg
import os
import numpy as np
import scipy as sp
import cv2
import skimage.io
import skimage.measure
import skimage.color
import scipy.misc
from skimage.util import img_as_float
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches



class NucleiSegmentation():
    """
    A class that utilize HistomicsTK to segment nuclei given a single tile or WSI tiles
    """
    def __init__(self):
        self.filedir = cfg.TILES_DIR
        #load ref image
        self.load_ref_img()


    def load_img(self, img_input):
        """
        :param img_input: np array  from tile wsi function
        :return:
        """
        self.img_tile = cv2.cvtColor(img_input, cv2.COLOR_RGBA2RGB)
        

    def load_ref_img(self):
        """
        :return: mean_ref , std_ref
        """
        # Load reference image for normalization
        self.ref_img_file = cfg.REF_IMG_NORM
        # load reference image for normalization
        self.im_reference = skimage.io.imread(self.ref_img_file)[:, :, :3]
        # get mean and stddev of reference image in lab space
        self.mean_ref, self.std_ref = htk.preprocessing.color_conversion.lab_mean_std(self.im_reference)


    def normalize_tile(self):
        """
        apply tile normalization based on reference image in config file
        :return: im_nmzd
        """
        # perform reinhard color normalization
        self.im_nmzd = htk.preprocessing.color_normalization.reinhard(self.img_tile, self.mean_ref, self.std_ref)
        #self.display_results(self.im_reference, "Reference Image", self.im_nmzd, "Normalized Input Image")


    def color_deconvolution(self):
        # perform standard color deconvolution
        self.im_stains = htk.preprocessing.color_deconvolution.color_deconvolution(self.im_nmzd, cfg.W).Stains
        #self.display_results(self.im_stains[:, :, 0], cfg.stain_1, self.im_stains[:, :, 1], cfg.stain_2)


    def display_results(self, img1 , title1, img2, title2):
        """
        Display the result for 2 images,
        Called by  normalize_tile and color_deconvolution
        :param img1: first img
        :param title1: title for first img
        :param img2: second img
        :param title2: title for 2nd img
        :return:
        """
        titlesize = 24

        # Display results
        plt.figure(figsize=(20, 10))

        plt.subplot(1, 2, 1)
        plt.imshow(img1)
        _ = plt.title(title1, fontsize=titlesize)

        plt.subplot(1, 2, 2)
        plt.imshow(img2)
        _ = plt.title(title2, fontsize=titlesize)
        plt.savefig("Color_Deconv.png")

    def segment_nuclei(self):
        """
        source: https://digitalslidearchive.github.io/HistomicsTK/examples/nuclei-segmentation.html
        :return:
        """
        # get nuclei/hematoxylin channel
        im_nuclei_stain = self.im_stains[:, :, 0]

        # segment foreground
        im_fgnd_mask = sp.ndimage.morphology.binary_fill_holes(
            im_nuclei_stain < cfg.FORGROUND_THRESHOLD)

        # run adaptive multi-scale LoG filter
        im_log_max, im_sigma_max = htk.filters.shape.cdog(
            im_nuclei_stain, im_fgnd_mask,
            sigma_min=cfg.MIN_RADIUS * np.sqrt(2),
            sigma_max=cfg.MAX_RADIUS * np.sqrt(2)
        )

        # detect and segment nuclei using local maximum clustering
        im_nuclei_seg_mask, seeds, maxima = htk.segmentation.nuclear.max_clustering(
            im_log_max, im_fgnd_mask, cfg.LOCAL_MAX_SEARCH_RADIUS)

        # filter out small objects
        self.im_nuclei_seg_mask = htk.segmentation.label.area_open(
            im_nuclei_seg_mask, cfg.MIN_NUCLEAUS_AREA).astype(np.int)

        # compute nuclei properties
        objProps = skimage.measure.regionprops(im_nuclei_seg_mask)

        #print('Number of nuclei = ', len(objProps))
        # number of nuclei
        return len(objProps)
        #self.display_nuclei_results(objProps )
        #plt.savefig("nuclei_results")
        


    def display_nuclei_results(self, objProps):
        # Display results
        titlesize = 24
        plt.figure(figsize=(20, 10))

        plt.subplot(1, 2, 1)
        plt.imshow(skimage.color.label2rgb(self.im_nuclei_seg_mask, self.img_tile, bg_label=0), origin='lower')
        plt.title('Nuclei segmentation mask overlay', fontsize=titlesize)

        plt.subplot(1, 2, 2)
        plt.imshow(self.img_tile)
        plt.xlim([0, self.img_tile.shape[1]])
        plt.ylim([0, self.img_tile.shape[0]])
        plt.title('Nuclei bounding boxes', fontsize=titlesize)

        for i in range(len(objProps)):
            c = [objProps[i].centroid[1], objProps[i].centroid[0], 0]
            width = objProps[i].bbox[3] - objProps[i].bbox[1] + 1
            height = objProps[i].bbox[2] - objProps[i].bbox[0] + 1

            cur_bbox = {
                "type": "rectangle",
                "center": c,
                "width": width,
                "height": height,
            }

            plt.plot(c[0], c[1], 'g+')
            mrect = mpatches.Rectangle([c[0] - 0.5 * width, c[1] - 0.5 * height],
                                       width, height, fill=False, ec='g', linewidth=2)
            plt.gca().add_patch(mrect)
            



    def compute_morphometry_feat(self, tile, coords, file_name , magnification):
        """
        This function will call the other function to compute morphometry features
        then it will append all tiles feature in one csv file named by wsi case.
        :param tile: numpy image
        :param coords: x, y tile coordinates
        :param file_name: original source WSI name
        :param magnification:
        :return:
        """
        # Construct tile fname
        tile_fname = str(coords[0]) + "_" + str(coords[1]) + "_" + str(magnification) + "_"

        self.load_img(tile)
        self.normalize_tile()
        self.color_deconvolution()
        num_nuclei = self.segment_nuclei()


        morphometry_dataframe = htk.features.compute_morphometry_features(self.im_nuclei_seg_mask)
        morphometry_dataframe['tile'] = tile_fname


        with open(os.path.join(cfg.FEATURES_DIR,file_name+"_"+str(cfg.MAGNIFICATION)+"_morpho.csv"), 'a') as f:
            morphometry_dataframe.to_csv(f, header=f.tell()==0, sep=',')
            
        return num_nuclei












