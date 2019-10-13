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
    fname = input("choose file name: ")
    fname = "mrcfiles/" + fname
    mrc = mrcfile.open(fname, mode='r+')
    img_matrix = np.copy(mrc.data)
    # for i in range(3):
    #     img_matrix = gaussian_filter(img_matrix, sigma=1, mode="constant", cval=0.0, truncate=4.0)
    threshold = img_matrix.mean()
    nx = mrc.header.nx
    ny = mrc.header.ny
    nz = mrc.header.nz
    shape = (nx, ny, nz)
    df = pd.DataFrame({"Name":[fname]})     # dateframe to record time
    df['region number'] = "4"

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


def gradient(mi, mj, param):       # calculate the gradient
    distance = math.sqrt(pow((mi.x_coordinate - mj.x_coordinate), 2) + pow((mi.y_coordinate - mj.y_coordinate), 2) + pow((mi.z_coordinate - mj.z_coordinate), 2))
    # distance **= 1/2
    param **= 2/3
    return abs(mj.density - mi.density)*param/distance


def max_distance_dif(M_regions):
    max_diff = 0
    for i in range(len(M_regions)-1):
        for j in range(i+1, len(M_regions)):
            distance = math.sqrt(pow((M_regions[i].x_coordinate - M_regions[j].x_coordinate), 2) + pow((M_regions[i].y_coordinate - M_regions[j].y_coordinate), 2) + pow((M_regions[i].z_coordinate - M_regions[j].z_coordinate), 2))
            if distance > max_diff:
                max_diff = distance
    max_diff **= 1/2
    return max_diff

def max_density_dif(M_regions):
    max_den_diff = 0
    for i in range(len(M_regions)-1):
        for j in range(i+1, len(M_regions)):
            den_diff = M_regions[i].density - M_regions[j].density
            if den_diff > max_den_diff:
                max_den_diff = den_diff
    return max_den_diff


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
    t1 = datetime.now()
    # sort array decreasing order based on density; stable
    tSortedArray = sorted(tArray, key=lambda voxel: voxel.density, reverse=True)
    t2 = datetime.now()
    delta = (t2 - t1).total_seconds()         # record sort time
    df['sorted']=(delta)
    # print("sorted done: " + str(delta))
    # record region number
    regionNum = -1
    # record region id:voxel with local maximum in corresponding region
    region_to_lm: List[Voxel] = []
    # create a dictionary with key as reigion id and value as list of voxels
    regions = dict()
    # find regions for each voxel
    size = matrix.size
    t1 = datetime.now()
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
    t2 = datetime.now()
    delta = (t2 - t1).total_seconds()
    df['get_region'] = delta
    return [region_to_lm, regions]


def merge_region(M_regions, regions, imgmatrix, n_regions):
    # t1 = datetime.now()
    # reorder M in increasing order based on density,
    # temp = sorted(M_regions, key=lambda voxel:voxel.density)     # mSorted=sorted(region_to_lm, key=lambda voxel: voxel.density)
    # M_regions = temp
    # t2 = datetime.now()
    # delta = t2 - t1
    # df['region sort']=delta
    # print("region sorte: "+ str(delta))
    # calculate the gaussian filtered image, modify sigma, truncate value to optimize
    count = 0
    # imgmatrix = img_matrix
    for i in range(2):
        t1 = datetime.now()
        imgmatrix = gaussian_filter(imgmatrix, sigma=1, mode="constant", cval=0.0, truncate=4.0)
        t2 = datetime.now()
        smooth = 'smooth' + str(i)
        delta = (t2 - t1).total_seconds()
        df[smooth] = delta
    t1 = datetime.now()
    for i in range(0, len(M_regions)):
        xc = M_regions[i].x_coordinate
        yc = M_regions[i].y_coordinate
        zc = M_regions[i].z_coordinate
        M_regions[i].density = imgmatrix[xc, yc, zc]
    t2 = datetime.now()
    delta = (t2 - t1).total_seconds()
    df['density updated'] = delta
    t1 = datetime.now()
    M_regions = sorted(M_regions, key=lambda voxel: voxel.density)
    t2 = datetime.now()
    delta = (t2 - t1).total_seconds()
    df['region sort']=delta
    max_density_diff = M_regions[-1].density - M_regions[0].density
    # print(max_density_diff)
    max_distance_diff = max_distance_dif(M_regions)
    # print(max_distance_diff)
    normal_parameter = max_distance_diff / max_density_diff
    t1 = datetime.now()
    while len(M_regions) > 2 * n_regions -1:
        # t1 = datetime.now()
        # #
        # t2 = datetime.now()
        # delta = t2 - t1
        # smooth = 'smooth' + str(count)
        # df[smooth]=delta
        # # print(smooth + ": " + str(delta))
        p = len(M_regions)
        q = int((1 + p) / 2)
        # t1 = datetime.now()
        # # update density after smoothing
        # for i in range(0, p):
        #     xc = M_regions[i].x_coordinate
        #     yc = M_regions[i].y_coordinate
        #     zc = M_regions[i].z_coordinate
        #     M_regions[i].density = imgmatrix[xc, yc, zc]
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
                gra[j] = gradient(M_regions[i], M_regions[j],normal_parameter)
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
            t5 = datetime.now()
        del M_regions[0:q]

    num_region = len(M_regions)
    if num_region > n_regions:
        num_pair = num_region - n_regions
        M_regions, regions = merge_pair(M_regions, regions, num_pair, normal_parameter)
    t2 = datetime.now()
    delta = (t2 - t1).total_seconds()
    df['merge']=delta
    return [M_regions, regions]

def merge_pair(M_regions, regions, num_pair, normal_param):
    num_region = len(M_regions)
    # M=[[0 for j in range(num_region)] for i in range(num_region)]
    M = []
    for i in range(num_region):
        steepascentlst = []
        for j in range(num_region):
            if (j <= i):
                steepascentlst.append(0)
            else:
                gra = gradient(M_regions[i], M_regions[j], normal_param)
                steepascentlst.append(gra)
        M.append(steepascentlst)
    print(num_pair)
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
    print(pairs)
    todelete = []
    for item in pairs:
        a, b = item[0], item[1]
        # print(a, b)
        if a < b:
            a, b = M_regions[a].regionID, M_regions[b].regionID
            # print("region id")
            # print(a, b)
            for v in regions[a]:
                v.regionID = b
            regions[b].extend(regions[a])
            regions = delete_key(a, regions)
            todelete.append(M_regions[item[0]])
        else:
            a, b = M_regions[a].regionID, M_regions[b].regionID
            # print("region id")
            # print(a, b)
            for v in regions[b]:
                v.regionID = a
            regions[a].extend(regions[b])
            regions = delete_key(b, regions)
            todelete.append(M_regions[item[1]])
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
        mrc_new = mrcfile.new('mrcfilestest/emd4297/oldalgo_smooth_2_param23_r4_rep1/{}'.format(fname), overwrite=True)
        mrc_new.set_data(np.zeros(shape, dtype=np.float32))
        mrc_new.voxel_size = mrc.voxel_size
        for v in regions[key]:
            mrc_new.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density
        mrc_new.close()


intialize()
n = int(input("n-- number of regions final: "))
t1 = datetime.now()
tArray = readData(img_matrix)
t2 = datetime.now()
delta = (t2-t1).total_seconds()
df['read data']=delta
region_to_lm, regions = gen_regions(img_matrix, tArray)
M_regions, regions = merge_region(region_to_lm, regions, img_matrix, n)
df.to_csv('mrcfilestest/emd4297/oldalgo_smooth_2_param23_r4_rep1/oldalgo_r4_rep1.csv')
outputregion(regions, shape)
mrc.close()
print('done')
