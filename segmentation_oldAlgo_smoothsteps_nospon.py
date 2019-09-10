from typing import List
import mrcfile
import numpy as np
from scipy.ndimage import gaussian_filter
import math
from datetime import datetime
import pandas as pd


# data structure holds voxel information
class Voxel(object):
    def __init__(self, x_coordinate, y_coordinate, z_coordinate, density, region_id=-1):
        self.x_coordinate = x_coordinate
        self.y_coordinate = y_coordinate
        self.z_coordinate = z_coordinate
        self.density = density
        self.regionID = region_id


# initialize program and ask user to input filename
def intialize():        # initialize program
    global mrc, img_matrix, shape, threshold, nx, ny, nz, df    # global variables
    fname = 'mrcfiles/emd4297.mrc'
    mrc = mrcfile.open(fname, mode='r+')
    img_matrix = np.copy(mrc.data)
    # img_matrix = gaussian_filter(img_matrix, sigma=1, mode="constant", cval=0.0, truncate=4.0)
    # img_matrix = gaussian_filter(img_matrix, sigma=1, mode="constant", cval=0.0, truncate=4.0)
    threshold = img_matrix.mean()
    nx = mrc.header.nx
    ny = mrc.header.ny
    nz = mrc.header.nz
    shape = (nx, ny, nz)
    # df = pd.DataFrame({"Name":[fname]})     # dateframe to record time


def delete_key(key, dictionary):        # delete key in dictionary
    d = dictionary
    del d[key]
    return d


def neighbors(matrix: object, x: object, y: object, z: object) -> object:   # return list of neighbor's coordinates
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
                if matrix[i, j, k] >= threshold:
                    # exclude itself
                    if (i, j, k) != (x, y, z):
                        neighbor.append((i, j, k))
    return neighbor


def gradient(mi, mj):       # calculate the gradient
    distance = pow((mi.x_coordinate - mj.x_coordinate), 2) + pow((mi.y_coordinate - mj.y_coordinate), 2) + pow((mi.z_coordinate - mj.z_coordinate), 2)
    distance = math.sqrt(distance)
    return (mj.density - mi.density)/distance


def readData(matrix):       # Read the data from mrc.data to voxel object and save in tArray
    tArray = []             # list holds all the voxels with density>=threshold
    # regionx = []          # list holds all the voxels with density<threshold
    row = matrix.shape[0]
    col = matrix.shape[1]
    dep = matrix.shape[2]
    temp_img = np.copy(matrix)
    for z in range(0, dep):
        for y in range(0, col):
            for x in range(0, row):
                density = temp_img[x, y, z]
                if density < threshold:
                    # vx = Voxel(x, y, z, density)
                    # regionx.append(vx)
                    density = 0
                v = Voxel(x, y, z, density)
                temp_img[x, y, z] = density
                tArray.append(v)
    return tArray


def gen_regions(matrix,tArray):      # generate regions for each block
    # t1 = datetime.now()
    # sort array decreasing order based on density; stable
    tSortedArray = sorted(tArray, key=lambda voxel: voxel.density, reverse=True)
    # t2 = datetime.now()
    # delta = t2 - t1         # record sort time
    # df['sorted']=(delta)
    # print("sorted done: " + str(delta))
    # record region number
    regionNum = -1
    # record region id:voxel with local maximum in corresponding region
    region_to_lm: List[Voxel] = []
    # create a dictionary with key as reigion id and value as list of voxels
    regions = dict()
    # find regions for each voxel
    size = matrix.size
    # t1 = datetime.now()
    for i in range(0, size):
        # print("number: " + str(i))
        v = tSortedArray[i]
        # only take voxels with density value larger than threshold
        if v.density > 0:
            # dictionary holds regionID and number of it, regionID:number of voxels
            regionRecord = dict()
            # Get the list of neighbors of voxel
            n = neighbors(matrix, v.x_coordinate, v.y_coordinate, v.z_coordinate)
            # in the neighbors, check the region id of each voxel
            for pos in n:
                # calculate the index in before sorted array to find the region id
                index = pos[0] + ny * pos[1] + ny * nz * pos[2]
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
    # t2 = datetime.now()
    # delta = t2 - t1
    # df['get_region'] = delta
    return [region_to_lm, regions]


def merge_region(M_regions, regions):
    # t1 = datetime.now()
    # reorder M in increasing order based on density,
    temp = sorted(M_regions, key=lambda voxel:voxel.density)     # mSorted=sorted(region_to_lm, key=lambda voxel: voxel.density)
    M_regions = temp
    # t2 = datetime.now()
    # delta = t2 - t1
    # df['region sort']=delta
    # print("region sorte: "+ str(delta))
    # calculate the gaussian filtered image, modify sigma, truncate value to optimize
    count = 0
    imgmatrix = img_matrix
    while len(M_regions) > 15:
        # t1 = datetime.now()
        imgmatrix = gaussian_filter(imgmatrix, sigma=1, mode="constant", cval=0.0, truncate=4.0)
        # t2 = datetime.now()
        # delta = t2 - t1
        # smooth = 'smooth' + str(count)
        # df[smooth]=delta
        # # print(smooth + ": " + str(delta))
        p = len(M_regions)
        q = int((1 + p) / 2)
        # t1 = datetime.now()
        # # update density after smoothing
        for i in range(0, p):
            xc = M_regions[i].x_coordinate
            yc = M_regions[i].y_coordinate
            zc = M_regions[i].z_coordinate
            M_regions[i].density = imgmatrix[xc, yc, zc]
        #
        # t2 = datetime.now()
        # delta = t2 - t1
        # densityupdate = 'density updated' + str(count)
        # df[densityupdate]=delta
        # print(densityupdate + ": " + str(delta))
        # t1 = datetime.now()
        for i in range(0, q):
            # t3 = datetime.now()
            gra = dict()
            for j in range(i + 1, p):
                # dictionary to keep all gradients for mi
                gra[j] = gradient(M_regions[i], M_regions[j])
            # find steepest ascent, the max value of g
            # print(gra)
            k = max(gra, key=gra.get)
            # t4 = datetime.now()
            # print(k)
            k = M_regions[k].regionID
            # merge region i to region k
            temp = M_regions[i].regionID
            for v in regions[temp]:
                v.regionID = k
            regions[k].extend(regions[temp])
            # remove checked regions
            regions = delete_key(temp, regions)
            # t5 = datetime.now()
        del M_regions[0:q]

    num_region = len(M_regions)
    if num_region > 8:
        num_pair = num_region - 8
        M_regions, regions = merge_pair(M_regions, regions, num_pair)
    return [M_regions, regions]

def merge_pair(M_regions, regions, num_pair):
    num_region = len(M_regions)
    # M=[[0 for j in range(num_region)] for i in range(num_region)]
    M = []
    for i in range(num_region):
        steepascentlst = []
        for j in range(num_region):
            if (j <= i):
                steepascentlst.append(0)
            else:
                gra = gradient(M_regions[i], M_regions[j])
                steepascentlst.append(gra)
        M.append(steepascentlst)

    pairs = []
    for i in range(num_pair):
        m = (0,0)
        for j in range(len(M)):
            for k in range(len(M[j])):
                if M[j][k] > M[m[0]][m[1]]:
                    m = (j,k)
        pairs.append(m)
        toreset = min(m[0],m[1])
        for j in range(len(M)):
            for k in range(len(M[j])):
                if j == toreset or k == toreset:
                    M[j][k]=0
    todelete = []
    for item in pairs:
        a, b = item[0], item[1]
        if a < b:
            for v in regions[a]:
                v.regionID = b
            regions[b].extend(regions[a])
            regions = delete_key(a, regions)
            todelete.append(M_regions[a])
        else:
            for v in regions[b]:
                v.regionID = a
            regions[a].extend(regions[b])
            regions = delete_key(b, regions)
            todelete.append(M_regions[b])
    temp_M = []
    for item in M_regions:
        if item not in todelete:
            temp_M.append(item)
    M_regions = temp_M
    return [M_regions, regions]


def outputregion(regions, shape):       # output segment regions
    rs = len(regions)
    for key in regions:
        fname='emdr'+str(key)+'.mrc'
        # generate mrcfile for each region
        mrc_new = mrcfile.new('mrcfilestest/emd4297smoothstepssingle_nospon_half/{}'.format(fname), overwrite=True)
        mrc_new.set_data(np.zeros(shape, dtype=np.float32))
        mrc_new.voxel_size = mrc.voxel_size
        for v in regions[key]:
            mrc_new.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density
        mrc_new.close()


intialize()
# t1 = datetime.now()
tArray = readData(img_matrix)
# t2 = datetime.now()
# delta = t2-t1
# df['read data']=delta
region_to_lm, regions = gen_regions(img_matrix, tArray)
M_regions, regions = merge_region(region_to_lm, regions)
# df.to_csv('emd4297smoothtwicesingleoutput.csv')
outputregion(regions, shape)
mrc.close()
print('done')
