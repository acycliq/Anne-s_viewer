"""
Segments a 3D image which has shape Z x channels x Height x Width.
The channels should be two. The Red channel (channel at position 0, ie Z x 0 x Height x Width) is the dapi
and channel at position 1 is the lamin stain
"""

import skimage.io
import numpy as np
import pandas as pd
from PIL import Image, ImageOps, ImageDraw
from skimage.transform import rescale, resize, downscale_local_mean
from cellpose import models
import diplib as dip
import os
import torch
import utils
import config
import logging


logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

# use_GPU = models.use_gpu()
# print('>>> GPU activated? %d'%use_GPU)
print(os.system('nvcc --version'))
print(os.system('nvidia-smi'))
torch.cuda.empty_cache()


def extract_borders_dip(label_image, offset_x=0, offset_y=0, ignore_labels=[0]):
    """
    Takes in a label_image and extracts the boundaries. It is assumed that the background
    has label = 0
    Parameters
    ----------
    label_image: the segmentation mask, a 2-D array
    offset_x: Determines how much the boundaries will be shifted by the on the x-axis
    offset_y: Determines how much the boundaries will be shifted by the on the y-axis
    ignore_labels: list of integers. If you want to ignore some masks, put the corresponding label in this list

    Returns
    -------
    Returns a dataframa with two columns, The first one is the mask label and the second column keeps the coordinates
    of the mask boundaries
    """
    labels = sorted(set(label_image.flatten()) - {0} - set(ignore_labels))
    cc = dip.GetImageChainCodes(label_image)  # input must be an unsigned integer type
    d = {}
    for c in cc:
        if c.objectID in labels:
            # p = np.array(c.Polygon())
            p = c.Polygon().Simplify()
            p = p + np.array([offset_x, offset_y])
            p = np.uint64(p).tolist()
            p.append(p[0])  # append the first pair at the end to close the polygon
            d[np.uint64(c.objectID)] = p
        else:
            pass
    df = pd.DataFrame([d]).T
    df = df.reset_index()
    df.columns = ['label', 'coords']
    return df


def normalize_img(img, normalize=[0,100]):
    perc1, perc2 = np.percentile(img, list(normalize))
    img = img - perc1
    img /= (perc2-perc1)
    img = img * 255.0
    img = np.clip(img, 0, 255)
    img = img.astype(np.uint16)
    return img


def draw_poly(img, polys, colours, ppm):
    img2 = img.copy()
    draw = ImageDraw.Draw(img2)
    for i, poly in enumerate(polys):
        poly = ppm * np.array(poly)
        poly = poly.tolist()
        poly = [tuple(d) for d in poly]
        draw.line(poly, fill=colours[i], width=3)
        # img3 = Image.blend(img, img2, 0.4)
    return img2


def zero_page(img, flag):
    height, width = img.shape
    c = 3
    out = np.zeros([c, height, width])
    if flag == 'lamin':
        # populate only the green channel
        out[1, :, :] = img[:, :]
    elif flag == 'dapi':
        # populate all three channels to make a grayscale image
        out[0, :, :] = img[:, :]
        out[1, :, :] = img[:, :]
        out[2, :, :] = img[:, :]
    else:
        raise RuntimeError('Unknown name %s' % flag)
    return out


def prepare_page(img, flag):
    assert len(img.shape) == 2
    if flag == 'lamin':
        arr = np.zeros((*img.shape, 3))
        arr[:, :, 1] = img
        out = Image.fromarray(arr.astype(np.uint8))
    elif flag == 'dapi':
        out = Image.fromarray(img.astype(np.uint8)).convert('RGB')
    else:
        raise RuntimeError('Unknown name %s' % flag)
    return out


def draw_page(page, boundaries, img, ppm, flag):
    polys = boundaries.coords.values
    colour = boundaries['colour'].values
    if len(polys) > 0:
        res = draw_poly(page, polys, colour, ppm)
        res.save('cyto_page_2.jpg')
        rgb_arr = np.array(res)
        out = rgb_arr.transpose(2, 0, 1)  # rgb_arr is Ly, Lx, Lz
    else:
        out = zero_page(img, flag)
    return out


def populate_pages(boundaries, img, ppm, flag):
    page = prepare_page(img, flag)
    out = draw_page(page, boundaries, img, ppm, flag)
    return out.astype(np.int16)


def which_channel(stain_type):
    if stain_type == 'dapi':
        out = 0
    elif stain_type == 'lamin':
        out = 1
    else:
        raise RuntimeError('Unknown stain_type %s' % stain_type)
    return out


def draw_boundaries(masks, tif, ppm, stain_type):
    c = 3  # three channels, rgb
    Lz, _, Ly, Lx = tif.shape
    out = np.zeros([Lz, c, Ly, Lx], dtype=np.int16)
    for i, mask in enumerate(masks):
        offset_x = 0
        offset_y = 0
        boundaries = extract_borders_dip(mask.astype(np.uint32), offset_x, offset_y, [0])
        boundaries['colour'] = utils.get_colour(boundaries.label.values)
        channel_id = which_channel(stain_type)
        img = normalize_img(tif[i, channel_id, :, :])  # channel:0 is dapi, channel: 1 the cytoplasm
        out[i, :, :, :] = populate_pages(boundaries, img, ppm, stain_type)

    return out


def stack_to_images(filename, series_name, stain_type):
    dapi = skimage.io.imread(filename)
    n = dapi.shape[0]
    for i in range(n):
        page = dapi[i, :, :, :]
        img = Image.fromarray(page.astype(np.uint8), 'RGB')
        target_file = os.path.join(os.path.join(cfg['segmented_cellpose'], series_name.split('.')[0], '%s_pages') % (stain_type), 'page_%s.jpg' % str(i).zfill(3))
        if not os.path.exists(os.path.dirname(target_file)):
            os.makedirs(os.path.dirname(target_file))
        img.save(target_file)
        logger.info('page_%d.jpg saved at %s' % (i, target_file))

        scale = 0.5
        img_resized = img.resize([int(scale * s) for s in img.size], Image.ANTIALIAS)
        target_file = os.path.join(os.path.join(cfg['segmented_cellpose'], series_name.split('.')[0], '%s_pages') % (stain_type), 'resized_page_%s.jpg' % str(i).zfill(3))
        img_resized.save(target_file)
        logger.info('resized_page_%d.jpg saved at %s' % (i, target_file))


def segment(img_3D, use_stiching=False):
    # DEFINE CELLPOSE MODEL
    model = models.Cellpose(gpu=True, model_type=config.cellpose_ini['model_type'])

    if use_stiching:
        masks, flows, styles, _ = model.eval(img_3D,
                                             channels=config.cellpose_ini['channels'],
                                             batch_size=config.cellpose_ini['batch_size'],
                                             diameter=config.cellpose_ini['diameter'],
                                             do_3D=False,
                                             mask_threshold=config.cellpose_ini['mask_threshold'],
                                             flow_threshold=config.cellpose_ini['flow_threshold'],
                                             stitch_threshold=0.5)
        np.savez('masks_2D_stiched.npz', masks)
    else:
        masks, flows, styles, diams = model.eval(img_3D,
                                                 channels=config.cellpose_ini['channels'],
                                                 batch_size=config.cellpose_ini['batch_size'],
                                                 diameter=config.cellpose_ini['diameter'],
                                                 anisotropy=config.cellpose_ini['anisotropy'],
                                                 # mask_threshold=cellpose_ini['mask_threshold'],
                                                 do_3D=True)

        np.savez('masks_rescaled_anisotropy_1.0.npz', masks)
    logger.info('Masks saved to disk!')
    return masks


def reshape_img(img):
    z, c, h, w = img.shape
    if c == 1:
        out = np.zeros([z, 2, h, w])
        out[:, 0, :, :] = img[:, 0, :, :]
    else:
        out = img
    return out


def resize_img(img, ppm):
    z, Ly, Lx = img.shape
    col = np.round(Ly/ppm)
    row = np.round(Lx/ppm)
    out = skimage.transform.resize(img, (z, col, row))
    out = out * 65535
    return out.astype(np.uint16)

def purge_channels(img):
    # Keep only channel 0 (this is the dapi) and channel 1 (this is lamin)
    # The others channels, if any, have the spots and are not needed for segmentation
    return img[:, :2, :, :]


def main(ppm, series_name, settings, use_stiching=False):
    big_tif_path = os.path.join(settings['img_dir'], series_name)
    # big_tif_path = settings['big_tif']

    # For multi - channel, multi-Z tiff's, cellpose expects the image as Z x channels x Ly x Lx.
    _img = skimage.io.imread(big_tif_path) # this image has at the fist channel the cyto and at the second the nucleus
    _img = purge_channels(_img)
    # _img = _img[:, ::-1, :, :]  # reverse the order so that Red channel is the dapi and Green the cyto
    # img = resize_img(_img, ppm)
    masks = segment(_img, use_stiching)
    dapi_boundaries = draw_boundaries(masks, _img, ppm, 'dapi')
    lamin_boundaries = draw_boundaries(masks, _img, ppm, 'lamin')

    # save the 3d tif
    dapi_target_file = os.path.join(settings['segmented_cellpose'], series_name.split('.')[0], 'cell_boundaries_dapi.tif')
    if not os.path.exists(os.path.dirname(dapi_target_file)):
        os.makedirs(os.path.dirname(dapi_target_file))
    skimage.io.imsave(dapi_target_file, dapi_boundaries)
    logger.info('cell boundaries dapi saved at %s' % dapi_target_file)

    lamin_target_file = os.path.join(settings['segmented_cellpose'], series_name.split('.')[0], 'cell_boundaries_lamin.tif')
    if not os.path.exists(os.path.dirname(lamin_target_file)):
        os.makedirs(os.path.dirname(lamin_target_file))
    skimage.io.imsave(lamin_target_file, lamin_boundaries)
    logger.info('cell boundaries dapi saved at %s' % lamin_target_file)
    return dapi_target_file, lamin_target_file


if __name__ == "__main__":
    cfg = config.atto425_DY520XL_MS002
    ppm = config.ppm  # pixels per micron
    for series in cfg['series']:
        dapi_target_file, lamin_target_file = main(ppm, series, cfg, use_stiching=True)
        stack_to_images(dapi_target_file, series, 'dapi')
        logger.info('ok, dapi done')

        stack_to_images(lamin_target_file, series, 'lamin')
        logger.info('ok, lamin done')

        logger.info('ok, all done')