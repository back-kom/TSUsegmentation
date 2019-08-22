import mrcfile
import numpy as np
import math
from datetime import datetime
from scipy.ndimage import gaussian_filter


class Voxel(object):
    def __init__(self, density = None, regionID = None, neighborlist = None, parent = None, children = None):
        self.density = density
        self.regionID = regionID
        self.neighborList = neighborlist
        self.parent = parent
        self.children = children

    def
