import os

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

ppm = 1  # Image Resolution: pixels per micron

cellpose_ini = {
    # 'model_type': 'cyto', # cyto or nuclei
    # 'channels': [2, 1],  # IF YOU HAVE GREEN=cytoplasm and RED=nucleus
    'model_type': 'nuclei',  # cyto or nuclei
    'channels': [0, 0],  # IF YOU HAVE G=cytoplasm and R=nucleus
    'diameter': 55.0,
    'batch_size': 2,
    'anisotropy': 1.0,
    'mask_threshold': -4.0,
    'flow_threshold': 4.0
}

atto425_DY520XL_MS002 = {
    'series': ['series_1.tif', 'series_2.tif', 'series_3.tif', 'series_4.tif', 'series_5.tif', 'series_4.tif', 'series_5.tif', 'series_6.tif',
               'series_7.tif', 'series_8.tif', 'series_9.tif'],

    # Set this to the path to your 3d folder with the image series
    'img_dir': os.path.join('E:\\', 'data', 'Anne', '220308 50umCF seq atto425 DY520XL MS002', 'img_series'),

    # destination folder where results will be saved at
    'segmented_cellpose': os.path.join('E:\\', 'data', 'Anne', '220308 50umCF seq atto425 DY520XL MS002', 'out'),
}

atto425_DY520XL_MS002_dapi_only = {
    # Set this to the path to your 3d dapi tiff
    'big_tif': os.path.join('E:\\', 'data', 'Anne', '220308 50umCF seq atto425 DY520XL MS002', 'matlab_output', 'background_image.tif'),

    # destination folder where results will be saved at
    'segmented_cellpose': os.path.join('E:\\', 'data', 'Anne', '220308 50umCF seq atto425 DY520XL MS002', 'out'),
}
