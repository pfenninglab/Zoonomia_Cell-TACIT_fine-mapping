#!/usr/bin/env python
import matplotlib
import numpy as np
import pandas as pd
import nd2, os, argparse
import pyarrow.feather as feather
from matplotlib import pyplot as plt
matplotlib.use("Agg")

import tensorflow as tf
from deepcell.utils.plot_utils import make_outline_overlay
from deepcell.applications.mesmer import Mesmer
tf.config.list_physical_devices('GPU') 

import scipy
import skimage
from skimage.measure import regionprops_table 

## load the mesmer model and prediction application
MESMER_DIR = '/home/bnphan/resources/deepcell/MultiplexSegmentation'
model = tf.keras.models.load_model(MESMER_DIR)
app = Mesmer(model)


def create_rgb_image(input_data, channel_colors):
    """Takes a stack of 1-, 2-, 3- channel data and converts it to an RGB image
    Args:
        input_data: 4D stack of images to be converted to RGB
        channel_colors: list specifying the color for each channel
    Returns:
        numpy.array: transformed version of input data into RGB version
    Raises:
        ValueError: if ``len(channel_colors)`` is not equal
            to number of channels
        ValueError: if invalid ``channel_colors`` provided
        ValueError: if input_data is not 4D, with 1 or 2 channels
    """
    if len(input_data.shape) != 4:
        raise ValueError('Input data must be 4D, '
                         'but provided data has shape {}'.format(input_data.shape))
    if input_data.shape[3] > 3:
        raise ValueError('Input data must have 1-3 channels, '
                         'but {} channels were provided'.format(input_data.shape[-1]))
    valid_channels = ['red', 'green', 'blue']
    channel_colors = [x.lower() for x in channel_colors]
    if not np.all(np.isin(channel_colors, valid_channels)):
        raise ValueError('Only red, green, or blue are valid channel colors')
    if len(channel_colors) != input_data.shape[-1]:
        raise ValueError('Must provide same number of channel_colors as channels in input_data')
    rgb_data = np.zeros(input_data.shape[:3] + (3,), dtype='float32')
    # rescale channels to aid plotting
    for img in range(input_data.shape[0]):
        for channel in range(input_data.shape[-1]):
            current_img = input_data[img, :, :, channel]
            non_zero_vals = current_img[np.nonzero(current_img)]
            # if there are non-zero pixels in current channel, we rescale
            if len(non_zero_vals) > 0:
                percentiles = np.percentile(non_zero_vals, [10, 95])
                rescaled_intensity = skimage.exposure.rescale_intensity(current_img,
                                                       in_range=(percentiles[0], percentiles[1]),
                                                       out_range='float32')
                # get rgb index of current channel
                color_idx = np.where(np.isin(valid_channels, channel_colors[channel]))
                rgb_data[img, :, :, color_idx] = rescaled_intensity
    # create a blank array for red channel
    return rgb_data


def filter_small_objects(seg, num_um, img_mpp, offset):
    """load the image to be segmented
        seg: the label mask of containing the objects
        num_um: the minimum diameter distance in microns to remove detected objects
        img_mpp: the micron per pixel scale value of the images
        offset: the integer number of objects previously used to start numbering new objects
    Returns:
        ret: the filtered segmented label mask, with numbered objects offseted
    """
    ## drop small objects with radius less than num_um
    ret = skimage.morphology.opening(seg[0,...,0], skimage.morphology.square(int(num_um/img_mpp)))
    ## relabel 
    ret = skimage.measure.label(ret)
    ret, _, _ = skimage.segmentation.relabel_sequential(
    	ret, offset = offset + 1)
    ## make same output shape
    ret = np.reshape(ret, seg.shape)
    return ret


def read_image(file):
    """load the image to be segmented
    get channels, DAPI = neun, FITC = nls-GFP, TRITC = nls-mcherry
        file: the ND2 file to be read in
    Returns:
        X: the 4-D numpy array of the image
        img_mpp: the micron per pixel value for the x dimension
    """
    f = nd2.ND2File(file) # read in image with metadata
    X = np.moveaxis(f.asarray(), 0, -1) # read to numpy array
    X = np.expand_dims(X, 0) # 
    img_mpp = f.voxel_size().x
    f.close() 
    return X, img_mpp


def segment_nuclei(X, img_mpp, app = app):
    seg1 = app.predict(X[...,[0,0]], image_mpp=img_mpp, compartment='nuclear')
    seg2 = app.predict(X[...,[1,1]], image_mpp=img_mpp, compartment='nuclear')
    seg3 = app.predict(X[...,[2,2]], image_mpp=img_mpp, compartment='nuclear')
    return seg1, seg2, seg3


def plot_segmentation(X, seg1, seg2, seg3, args):
    ## create the overlay from the 3-channel image and the segmentations
    img = create_rgb_image(X, channel_colors = ['blue', 'green', 'red']) * .7
    overlay_seg1 = make_outline_overlay(rgb_data=img, predictions=seg1)
    overlay_seg2 = make_outline_overlay(rgb_data=img, predictions=seg2)
    overlay_seg3 = make_outline_overlay(rgb_data=img, predictions=seg3)
    #$ make the plot
    fig, ax = plt.subplots(1, 3, figsize=(90, 45))
    ax[1].imshow(overlay_seg1[0, ... ])
    ax[0].imshow(overlay_seg2[0, ... ])
    ax[2].imshow(overlay_seg3[0, ... ])
    ax[1].set_title('NeuN Segmentation')
    ax[0].set_title('nls-GFP Segmentation')
    ax[2].set_title('mCherry Segmentation')
    ## save the plot
    save_fn = os.path.splitext(os.path.basename(args.file))[0] + '_segmented_objects.png'
    save_fn = os.path.join(args.out_dir, 'segmentation_plot' ,save_fn)
    os.makedirs(os.path.join(args.out_dir, 'segmentation_plot'), exist_ok = True)
    plt.savefig(save_fn, bbox_inches='tight',pad_inches = 0)
    return


def main(args):
    ## read in the image
    X, img_mpp = read_image(args.file)

    ## perform Mesmer nuclei segmentation 
    seg_neun, seg_nlsgfp, seg_mcherry = segment_nuclei(X, img_mpp = img_mpp)

    ## NeuN objects tend to be smaller than GFP or mCherry nuclei
    ## expand them to ease the colocalization
    seg_neun = skimage.segmentation.expand_labels(seg_neun, distance= int(1./img_mpp))

    ## preprocess the segmented nuclei across 3 channels
    seg_neun = filter_small_objects(seg_neun, 3, img_mpp, np.max(seg_neun))
    seg_nlsgfp = filter_small_objects(seg_nlsgfp, 3, img_mpp, np.max(seg_neun))
    seg_mcherry = filter_small_objects(seg_mcherry, 3, img_mpp, np.max(seg_nlsgfp))

    ## extract coordinates, features, and fluorescence intensities
    prop1 = ['label', 'area', 'solidity', 'centroid', 'intensity_mean', 'intensity_max']
    rename_col_intensity = { 'centroid-0': 'x',  'centroid-1': 'y',
    'intensity_mean-0': "mean_NeuN", 'intensity_mean-1':"mean_GFP",  'intensity_mean-2':'mean_mCherry',
    'intensity_max-0': "max_NeuN", 'intensity_max-1':"max_GFP", 'intensity_max-2':'max_mCherry'
    }

    df_neun= pd.DataFrame(regionprops_table(seg_neun[0,...,0], X[0,...], properties = prop1))
    df_gfp = pd.DataFrame( regionprops_table( seg_nlsgfp[0,...,0], X[0,...], properties = prop1))
    df_mcherry = pd.DataFrame(regionprops_table(seg_mcherry[0,...,0], X[0,...], properties = prop1))

    df_neun["Segmentation"] = 'NeuN'
    df_gfp["Segmentation"] = 'GFP'
    df_mcherry["Segmentation"] = 'mCherry'

    df = pd.concat([df_neun, df_gfp, df_mcherry], ignore_index = True)
    df = df.rename(columns = rename_col_intensity)

    ## extract proportion of objects overlapping NeuN masks
    neun_mask = seg_neun >0
    neun_mask = neun_mask[0,...]
    prop2 = ['label', 'intensity_mean']
    rename_col_overlap = {"intensity_mean-0": "percent_NeuN"}

    neun_and_neun = pd.DataFrame(regionprops_table(seg_neun[0,...,0], neun_mask, properties=prop2))
    nlsGFP_and_neun = pd.DataFrame(regionprops_table(seg_nlsgfp[0,...,0], neun_mask, properties=prop2))
    mCherry_and_neun = pd.DataFrame(regionprops_table(seg_mcherry[0,...,0], neun_mask, properties=prop2))

    df_overlap = pd.concat([neun_and_neun, nlsGFP_and_neun, mCherry_and_neun], ignore_index = True)
    df_overlap = df_overlap.rename(columns=rename_col_overlap)

    df2 = pd.merge(df, df_overlap, on= 'label')
    df2['file_name'] = os.path.basename(args.file)
    df2['prefix'] = args.file

    ## save the objects
    save_fn = os.path.splitext(os.path.basename(args.file))[0] + '_segmented_objects.feather'
    save_fn = os.path.join(args.out_dir, 'tables' ,save_fn)
    os.makedirs(os.path.join(args.out_dir, 'tables'), exist_ok = True)
    df2.reset_index().to_feather(save_fn)

    ## make plots
    plot_segmentation(X, seg_neun, seg_nlsgfp, seg_mcherry, args)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Segment nuclei across 3-channel fluorescence microscopy images from Nikon ND2 files.')
    parser.add_argument('--file',type=str, help='an ND2 3-channel microscopy file with fluorescence in the DAPI, FITC, and TRITC channels',required=True)
    parser.add_argument('--prefix',type=str, default='', help="Perhaps a description of the experiment.")
    parser.add_argument('--out-dir',type=str, default='.', help="Output directory to store the resultant dataframe")

    # args = parser.parse_args(['--file', 'data/raw_data/reporter_assay/images/Clarified/Slide01.tile3x3.w20x - Clarified.nd2.nd2',
    #     '--prefix', 'testing', '--out-dir', 'data/raw_data/reporter_assay'])
    args = parser.parse_args()
    main(args)