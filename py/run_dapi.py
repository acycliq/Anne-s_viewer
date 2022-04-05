"""
Segments a 3D dapi image which has shape Z x Height x Width.
"""

import skimage.io
import numpy as np
import pandas as pd
import cv2
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


def draw_poly(img, polys, colours):
    img2 = img.copy()
    draw = ImageDraw.Draw(img2)
    for i, poly in enumerate(polys):
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


def draw_page(page, boundaries, img, flag):
    polys = boundaries.coords.values
    colour = boundaries['colour'].values
    if len(polys) > 0:
        res = draw_poly(page, polys, colour)
        res.save('cyto_page_2.jpg')
        rgb_arr = np.array(res)
        out = rgb_arr.transpose(2, 0, 1)  # rgb_arr is Ly, Lx, Lz
    else:
        out = zero_page(img, flag)
    return out


def populate_pages(boundaries, img, flag):
    page = prepare_page(img, flag)
    out = draw_page(page, boundaries, img, flag)
    return out.astype(np.int16)


def which_channel(stain_type):
    if stain_type == 'dapi':
        out = 0
    elif stain_type == 'lamin':
        out = 1
    else:
        raise RuntimeError('Unknown stain_type %s' % stain_type)
    return out


def draw_boundaries(masks, tif, stain_type):
    c = 3  # three channels, rgb
    Lz, Ly, Lx = tif.shape
    out = np.zeros([Lz, c, Ly, Lx], dtype=np.int16)
    for i, mask in enumerate(masks):
        offset_x = 0
        offset_y = 0
        boundaries = extract_borders_dip(mask.astype(np.uint32), offset_x, offset_y, [0])
        boundaries['colour'] = utils.get_colour(boundaries.label.values)
        img = normalize_img(tif[i, :, :])  # tif is a dapi image
        out[i, :, :, :] = populate_pages(boundaries, img, stain_type)

    return out


def stack_to_images(filename, stain_type):
    dapi = skimage.io.imread(filename)
    n = dapi.shape[0]
    for i in range(n):
        page = dapi[i, :, :, :]
        img = Image.fromarray(page.astype(np.uint8), 'RGB')
        target_file = os.path.join(os.path.join(cfg['segmented_cellpose'], '%s_pages') % (stain_type), 'page_%s.jpg' % str(i).zfill(3))
        if not os.path.exists(os.path.dirname(target_file)):
            os.makedirs(os.path.dirname(target_file))
        img.save(target_file)
        logger.info('page_%d.jpg saved at %s' % (i, target_file))

        scale = 0.5
        img_resized = img.resize([int(scale * s) for s in img.size], Image.ANTIALIAS)
        target_file = os.path.join(os.path.join(cfg['segmented_cellpose'], '%s_pages') % (stain_type), 'resized_page_%s.jpg' % str(i).zfill(3))
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
    z, h, w = img.shape
    out = np.zeros([z, 2, h, w])
    out[:, 0, :, :] = img[:, :, :]
    logger.info('Image has now shape: z-by-c-by-height-by-width = [%d-by-%d-by-%d-by-%d]' %
                (out.shape[0], out.shape[1], out.shape[2], out.shape[3])
                )

    return out.astype(np.uint16)


def downsize_img(img, shrink_factor):
    z, Ly, Lx = img.shape
    col = np.round(Ly/shrink_factor)
    row = np.round(Lx/shrink_factor)
    out = skimage.transform.resize(img, (z, col, row))
    out = out * 65535
    return out.astype(np.uint16)


def restore_shape(label_image, shape):
    masks = [np.array(Image.fromarray(d).resize(shape, Image.NEAREST)) for d in label_image]
    return np.stack(masks)


def preprocess_img(img, shrink_factor):
    if shrink_factor != 1.0:
        img = downsize_img(img, shrink_factor)
    img = reshape_img(img)
    return img


def postprocess(masks, shape):
    """
    shape must be (width, height)
    """
    out = restore_shape(masks, shape)
    np.savez('masks_2D_stiched_fullsize.npz', out)
    logger.info('Full sized masks saved to disk!')
    return out


def main(shrink_factor, settings, use_stiching=False):
    big_tif_path = settings['big_tif']

    # For multi - channel, multi-Z tiff's, cellpose expects the image as Z x channels x Ly x Lx.
    img = skimage.io.imread(big_tif_path)
    z, h, w = img.shape
    _img = preprocess_img(img, shrink_factor)
    _masks = segment(_img, use_stiching)
    masks = postprocess(_masks, (w, h))
    dapi_boundaries = draw_boundaries(masks, img, 'dapi')

    # save the 3d tif
    dapi_target_file = os.path.join(settings['segmented_cellpose'], 'dapi', 'cell_boundaries_dapi.tif')
    if not os.path.exists(os.path.dirname(dapi_target_file)):
        os.makedirs(os.path.dirname(dapi_target_file))
    skimage.io.imsave(dapi_target_file, dapi_boundaries)
    logger.info('cell boundaries dapi saved at %s' % dapi_target_file)

    return dapi_target_file


if __name__ == "__main__":
    cfg = config.atto425_DY520XL_MS002_dapi_only
    shrink_factor = config.shrink_factor  # pixels per micron
    dapi_target_file = main(shrink_factor, cfg, use_stiching=True)
    stack_to_images(dapi_target_file, 'dapi')
    logger.info('ok, dapi done')

    logger.info('ok, all done')