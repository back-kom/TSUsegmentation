from typing import List, Union
import mrcfile
import scipy
import numpy as np
from mrcfile.bzip2mrcfile import Bzip2MrcFile
from mrcfile.gzipmrcfile import GzipMrcFile
from mrcfile.mrcfile import MrcFile
from scipy.ndimage import gaussian_filter
import math

# mrc = mrcfile.open('mrcfiles/emd4297.mrc')
# test_mrc = mrcfile.new('mrcfiles/new_mrc.mrc', overwrite=True)
# test_mrc.set_data(mrc.data)
#
# a = gaussian_filter(test_mrc.data, sigma=5, mode="constant", cval=0.0, truncate=4.0)
# print(a)

