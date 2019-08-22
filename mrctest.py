from typing import List, Union
import mrcfile
import scipy
import numpy as np
from mrcfile.bzip2mrcfile import Bzip2MrcFile
from mrcfile.gzipmrcfile import GzipMrcFile
from mrcfile.mrcfile import MrcFile
from scipy.ndimage import gaussian_filter
import math


# data structure holds voxel information
class Voxel(object):
    def __init__(self, x_coordinate, y_coordinate, z_coordinate, density, region_id=-1):
        self.x_coordinate = x_coordinate
        self.y_coordinate = y_coordinate
        self.z_coordinate = z_coordinate
        self.density = density
        self.regionID = region_id


# class Regionlm(object):
#     def __init__(self, rid, lmax):
#         self.rid = rid
#         self.lmax = lmax


# delete key in dictionary
def delete_key(key, dictionary):
    d = dictionary
    del d[key]
    return d


# return list of neighbor's coordinates
def neighbors(matrix: object, x: object, y: object, z: object) -> object:
    # initialize list of neighbors
    neighbor: List[tuple] = []
    # get x boundary
    row = len(matrix)
    # get y boundary value
    col = len(matrix[0])
    # get z boundary value
    dep = len(matrix[0][0])
    # loop to find neighbors coordinates, index must greater or equal to 0, and less or equal to the boundary value
    for k in range(max(0, z-1), min(row, z+2)):
        for j in range(max(0, y-1), min(col, y+2)):
            for i in range(max(0, x-1), min(dep, x+2)):
                # exclude itself
                if (i, j, k) != (x, y, z):
                    neighbor.append((i, j, k))
    return neighbor


# calculate the gradient
def gradient(mi, mj):
    distance = pow((mi.x_coordinate - mj.x_coordinate), 2) + pow((mi.y_coordinate - mj.y_coordinate), 2) + pow((mi.z_coordinate - mj.z_coordinate), 2)
    return (mj.density - mi.density)/distance


# assume all mrc files with coordinates starting from 0;
mrc = mrcfile.open('mrcfiles/emd4297.mrc', mode='r+')
# create a new mrc file writable for data change
img_matrix = np.copy(mrc.data)
print
# set threshold value as the mean density
threshold = mrc.data.mean()
nx = mrc.header.nx
ny = mrc.header.ny
nz = mrc.header.nz

tArray = []

# Read the data from mrc.d ata to voxel object and save in tArray
for z in range(0, nz):
    print('z ', z)
    for y in range(0, ny):
        for x in range(0, nx):
            density = img_matrix[x, y, z]
            if density < threshold:
                density = 0
            v = Voxel(x, y, z, density)
            img_matrix[x, y, z] = density
            tArray.append(v)

size = img_matrix.size
print(size)

# sort array decreasing order based on density; stable
tSortedArray = sorted(tArray, key=lambda voxel: voxel.density, reverse=True)

# record region number
regionNum = -1
# record region id:voxel with local maximum in corresponding region
region_to_lm: List[Voxel] = []

# create a dictionary with key as reigion id and value as list of voxels
regions = dict()
# find regions for each voxel
for i in range(0, size):
    v = tSortedArray[i]
    # only take voxels with density value larger than threshold
    if v.density > 0:
        # dictionary holds regionID and number of it, regionID:number of voxels
        regionRecord = dict()
        # Get the list of neighbors of voxel
        n = neighbors(img_matrix, v.x_coordinate, v.y_coordinate, v.z_coordinate)
        # in the neighbors, check the region id of each voxel
        for pos in n:
            # calculate the index in before sorted array to find the region id
            index = pos[0] + ny*pos[1] + ny*nz*pos[2]
            # get the region id of neighbor
            rId = tArray[index].regionID
            # if neighbors' region id exists
            if rId != -1:
                # check if this region id in the dictionary, if not assign region id with value 1, if yes then
                # increase the value of region id
                if rId in regionRecord:
                    regionRecord[rId] += 1
                else:
                    regionRecord[rId] = 1
        # check region dictionary, if it is empty then assign new region id
        if len(regionRecord) == 0:
            # increase region id record, and assign the voxel with this new region id
            regionNum += 1
            v.regionID = regionNum
            # rlm = Regionlm(regionNum, v)
            region_to_lm.insert(regionNum, v)
            regions[regionNum] = []
            regions[regionNum].append(v)
        # if region dictionary has item, assign the voxel with existing region id in region record
        # the id will be the key with maximum value in the region dictionary
        else:
            r = max(regionRecord, key=regionRecord.get)
            v.regionID = r
            regions[r].append(v)
    # if voxel density = 0, then exit for loop, as the sorted array in decrease order, the density voxel afterwards
    # are all zero, no need to check
    else:
        break

# img_matrix = gaussian_filter(img_matrix, sigma=1, mode="constant", cval=0.0, truncate=1.0)
# emdrg1 = mrcfile.new('mrcfiles/emd4297g1.mrc', overwrite=True)
# emdrg1.set_data(img_matrix)
# emdrg1.close()

# mrc.close()
# print(len(region_to_lm))

#set the number of smoothing steps
tStep = 1
count = 0
# reorder M in increasing order based on density,
region_to_lm.reverse()  # mSorted=sorted(region_to_lm, key=lambda voxel: voxel.density)
# calculate the gaussian filtered image, modify sigma, truncate value to optimize
while count < tStep:
    #img_matrix = gaussian_filter(img_matrix, sigma=5, mode="constant", cval=0.0, truncate=4.0)
    p = len(region_to_lm)
    q = int((1 + p)/2)
    # update density after smoothing
    # for i in range(0, p):
    #     xc = region_to_lm[i].x_coordinate
    #     yc = region_to_lm[i].y_coordinate
    #     zc = region_to_lm[i].z_coordinate
    #     region_to_lm[i].density = img_matrix[xc, yc, zc]


    for i in range(0, q):
        for j in range(i+1, p):
            # dictionary to keep all gradients for mi
            gra = dict()
            gra[j] = gradient(region_to_lm[i], region_to_lm[j])
        # find steepest ascent, the max value of g
        k = max(gra, key=gra.get)
        k = region_to_lm[k].regionID
        # merge region i to region k
        temp = region_to_lm[i].regionID
        for v in regions[temp]:
            v.regionID = k
        regions[k].extend(regions[temp])
        # remove checked regions
        regions = delete_key(temp, regions)
    del region_to_lm[0:q]
    count += 1

print(len(regions))

shape = (nx, ny, nz)
emdr0 = mrcfile.new('mrcfiles/emd4297r0.mrc', overwrite=True)
emdr0.set_data(np.zeros(shape, dtype=np.float32))
for v in regions[0]:
    emdr0.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density

emdr1 = mrcfile.new('mrcfiles/emd4297r1.mrc', overwrite=True)
emdr1.set_data(np.zeros(shape, dtype=np.float32))
for v in regions[1]:
    emdr1.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density

emdr2 = mrcfile.new('mrcfiles/emd4297r2.mrc', overwrite=True)
emdr2.set_data(np.zeros(shape, dtype=np.float32))
for v in regions[2]:
    emdr2.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density

emdr3 = mrcfile.new('mrcfiles/emd4297r3.mrc', overwrite=True)
emdr3.set_data(np.zeros(shape, dtype=np.float32))
for v in regions[3]:
    emdr3.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density

emdr4 = mrcfile.new('mrcfiles/emd4297r4.mrc', overwrite=True)
emdr4.set_data(np.zeros(shape, dtype=np.float32))
for v in regions[4]:
    emdr4.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density

emdr5 = mrcfile.new('mrcfiles/emd4297r5.mrc', overwrite=True)
emdr5.set_data(np.zeros(shape, dtype=np.float32))
for v in regions[5]:
    emdr5.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density

emdr6 = mrcfile.new('mrcfiles/emd4297r6.mrc', overwrite=True)
emdr6.set_data(np.zeros(shape, dtype=np.float32))
for v in regions[6]:
    emdr6.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density

emdr7 = mrcfile.new('mrcfiles/emd4297r7.mrc', overwrite=True)
emdr7.set_data(np.zeros(shape, dtype=np.float32))
for v in regions[7]:
    emdr7.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density

emdr8 = mrcfile.new('mrcfiles/emd4297r8.mrc', overwrite=True)
emdr8.set_data(np.zeros(shape, dtype=np.float32))
for v in regions[8]:
    emdr8.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density


emdr0.close()
emdr1.close()
emdr2.close()
emdr3.close()
emdr4.close()
emdr5.close()
emdr6.close()
emdr7.close()
emdr8.close()
mrc.close()

# tArray[0]=tArray[0]+(j,)
# print(tArray)
