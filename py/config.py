import os

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

shrink_factor = 1  # How much to rescale the image by

cellpose_ini = {
    # 'model_type': 'cyto', # cyto or nuclei
    # 'channels': [2, 1],  # IF YOU HAVE GREEN=cytoplasm and RED=nucleus
    'model_type': 'nuclei',  # cyto or nuclei
    'channels': [0, 0],  # IF YOU HAVE dapi only
    'diameter': 55.0,
    'batch_size': 2,
    'anisotropy': 1.0,
    'mask_threshold': -4.0,
    'flow_threshold': 4.0
}

atto425_DY520XL_MS002 = {
    # Set this to the path to your multichannel 3d tiff
    'big_tif': os.path.join('E:\\', 'data', 'Anne', '220308 50umCF seq atto425 DY520XL MS002', 'img_series', 'series_1.tif'),

    # destination folder where results will be saved at
    'segmented_cellpose': os.path.join('E:\\', 'data', 'Anne', '220308 50umCF seq atto425 DY520XL MS002', 'out'),
}

atto425_DY520XL_MS002_dapi_only = {
    # Set this to the path to your 3d dapi tiff
    'big_tif': os.path.join('E:\\', 'data', 'Anne', '220308 50umCF seq atto425 DY520XL MS002', 'matlab_output', 'background_image.tif'),

    # destination folder where results will be saved at
    'segmented_cellpose': os.path.join('E:\\', 'data', 'Anne', '220308 50umCF seq atto425 DY520XL MS002', 'out'),
}
