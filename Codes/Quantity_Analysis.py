import os
import mrcfile
import numpy as np


# data structure holds voxel information
class Voxel(object):
    def __init__(self, x_coordinate, y_coordinate, z_coordinate, density, region_id=-1):
        self.x_coordinate = x_coordinate
        self.y_coordinate = y_coordinate
        self.z_coordinate = z_coordinate
        self.density = density


# initialize program and ask user to input filename
def intialize():        # initialize program
    global testregion, template_region, shape, threshold    # global variables
    template_region = []
    testregion = input("choose region file:")
    tempdirectory = input("choose the template directory:")
    threshold = input("please choose the threshold:")
    directory = os.fsencode(tempdirectory)
    for file in os.listdir(directory):
        filename = os.fsencode(file)
        if filename.endswith(".mrc"):
            template_region.append(filename)
            continue

    testregion_mrc = mrcfile.open(testregion, mode='r+')
    img_matrix = np.copy(testregion_mrc.data)
    threshold = img_matrix.mean()
    nx = mrc.header.nx
    ny = mrc.header.ny
    nz = mrc.header.nz
    shape = (nx, ny, nz)
    df = pd.DataFrame({"Name":[fname]})     # dateframe to record time
