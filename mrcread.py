import numpy as np
import mrcfile
import os

with mrcfile.new('tmp.mrc') as mrc:
    mrc.set_data(np.zeros((5,5),dtype=np.int8))
    mrc.data[1:4,1:4]=10
    mrc.data
mrc.close()
