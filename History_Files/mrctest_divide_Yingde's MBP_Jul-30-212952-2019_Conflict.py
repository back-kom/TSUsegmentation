from typing import List, Union
import mrcfile
import scipy
import numpy as np
from mrcfile.bzip2mrcfile import Bzip2MrcFile
from mrcfile.gzipmrcfile import GzipMrcFile
from mrcfile.mrcfile import MrcFile
from scipy.ndimage import gaussian_filter
import math
from datetime import datetime
import pandas as pd


# data structure holds voxel information
class Voxel(object):
    def __init__(self, x_coordinate, y_coordinate, z_coordinate, density, region_id=-1, cube_id=None):
        self.x_coordinate = x_coordinate
        self.y_coordinate = y_coordinate
        self.z_coordinate = z_coordinate
        self.density = density
        self.regionID = region_id
        self.cubeID = cube_id


# class Regionlm(object):
#     def __init__(self, rid, lmax):
#         self.rid = rid
#         self.lmax = lmax

# initialize program and ask user to input filename
def intialize():        # initialize program
    global mrc, img_matrix, shape, nx, ny, nz, df, unit, global_regionid, cube_id
    fname = input("choose mrc file:")
    mrc = mrcfile.open(fname, mode='r+')
    img_matrix = np.copy(mrc.data)
    # threshold = img_matrix.mean()
    nx = mrc.header.nx
    ny = mrc.header.ny
    nz = mrc.header.nz
    shape = (nx, ny, nz)
    unit = int(math.sqrt(nx))
    df = pd.DataFrame({"Name":[fname]})
    global_regionid = -1


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
                # exclude itself
                if (i, j, k) != (x, y, z):
                    neighbor.append((i, j, k))
    return neighbor


def gradient(mi, mj):       # calculate the gradient
    distance = pow((mi.x_coordinate - mj.x_coordinate), 2) + pow((mi.y_coordinate - mj.y_coordinate), 2) + pow((mi.z_coordinate - mj.z_coordinate), 2)
    return (mj.density - mi.density)/distance


def readData(matrix, threshold, cubeid):       # Read the data from mrc.data to voxel object and save in tArray
    tArray = []
#    regionx = []
    row = matrix.shape[0]
    col = matrix.shape[1]
    dep = matrix.shape[2]
    temp_img = np.copy(matrix)
    for z in range(0, dep):
        #print("z value:" + str(z))
        for y in range(0, col):
            for x in range(0, row):
                density = temp_img[x, y, z]
                print(x, y, z, density)
                if density < threshold:
                    # vx = Voxel(x, y, z, density, id)
                    # regionx.append(vx)
                    density = 0
                v = Voxel(x, y, z, density, -1, cubeid)
                temp_img[x, y, z] = density
                tArray.append(v)
    return tArray


def gen_regions(matrix,tArray, g_regionid):      # generate regions for each block
    t1 = datetime.now()
    # sort array decreasing order based on density; stable
    tSortedArray = sorted(tArray, key=lambda voxel: voxel.density, reverse=True)
    t2 = datetime.now()
    delta = t2 - t1
    df['sort']=(delta)
    #print("sorted done: " + str(delta))
    global global_regionid
    global_regionid = g_regionid
    # record region number
    regionNum = global_regionid
    # record region id:voxel with local maximum in corresponding region
    region_to_lm: List[Voxel] = []
    # create a dictionary with key as reigion id and value as list of voxels
    regions = dict()
    # find regions for each voxel
    size = matrix.size
    temp_ny = matrix.shape[1]
    temp_nz = matrix.shape[2]
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
                index = pos[0] + temp_ny * pos[1] + temp_ny * temp_nz * pos[2]
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
            global_regionid = regionNum
        # if voxel density = 0, then exit for loop, as the sorted array in decrease order, the density voxel afterwards
        # are all zero, no need to check
        else:
            break
    t2 = datetime.now()
    delta = t2 - t1
    df['get_region'] =delta
    # print("generate regions: " + str(delta))
    return [region_to_lm, regions]


def merge_region(region_to_lm, regions, temp_image):
    t1 = datetime.now()
    # reorder M in increasing order based on density,
    region_to_lm.reverse()     # mSorted=sorted(region_to_lm, key=lambda voxel: voxel.density)
    t2 = datetime.now()
    delta = t2 - t1
    df['region sort']=delta
    #print("region sorted: "+ str(delta))
    # calculate the gaussian filtered image, modify sigma, truncate value to optimize
    count = 0
    imgmatrix = temp_image
    while len(region_to_lm) > 2:
        #print("ok" + str(count))
        t1 = datetime.now()
        imgmatrix = gaussian_filter(imgmatrix, sigma=1, mode="constant", cval=0.0, truncate=4.0)
        t2 = datetime.now()
        delta = t2 - t1
        smooth = 'smooth ' + str(count)
        df[smooth]=delta
        #print(smooth + ": " + str(delta))
        p = len(region_to_lm)
        q = int((1 + p) / 2)
        t1 = datetime.now()
        # update density after smoothing
        for i in range(0, p):
            xc = region_to_lm[i].x_coordinate
            yc = region_to_lm[i].y_coordinate
            zc = region_to_lm[i].z_coordinate
            if xc>220 and yc > 220 and zc > 220:
                print(xc, yc, zc)
            region_to_lm[i].density = imgmatrix[xc, yc, zc]

        t2 = datetime.now()
        delta = t2 - t1
        densityupdate = 'density updated' + str(count)
        df[densityupdate]=delta
        #print(densityupdate + ": " + str(delta))
        t1 = datetime.now()
        for i in range(0, q):
            gra = dict()
            for j in range(i + 1, p):
                # dictionary to keep all gradients for mi
                gra[j] = gradient(region_to_lm[i], region_to_lm[j])
            # find steepest ascent, the max value of g
            # print(gra)
            k = max(gra, key=gra.get)
            # print(k)
            k = region_to_lm[k].regionID
            # merge region i to region k
            temp = region_to_lm[i].regionID
            for v in regions[temp]:
                v.regionID = k
            regions[k].extend(regions[temp])
            # remove checked regions
            regions = delete_key(temp, regions)
        del region_to_lm[0:q]
        t2 = datetime.now()
        delta = t2 - t1
        merge = 'merged'+str(count)
        df[merge]=delta
        #print(merge + ": " + str(delta))
        count += 1
        #print("count:" + str(count))
    #  print("number of regions:" + str(len(region_to_lm)))
    # print("merge finished")
    # for v in region_to_lm:
    #     print(v.regionID)
    # for key in regions:
    #     for v in regions[key]:
    #         print(v.regionID, v.x_coordinate, v.y_coordinate, v.z_coordinate)
    return [region_to_lm, regions]


def reverse_coordinate(temp_M_region, temp_regions, i_step, j_step, k_step):
    for v in temp_M_region:
        v.x_coordinate += i_step
        v.y_coordinate += j_step
        v.z_coordinate += k_step
        if v.x_coordinate > 220 and v.y_coordinate > 220 and v.z_coordinate > 220:
            print("m0 update")
            print(v.x_coordinate, v.y_coordinate, v.z_coordinate)
    rs = len(temp_regions)
    # print("number of regions:" + str(rs))

    if rs > 0:
        for key in temp_regions:
            length = len(temp_regions[key])
            list_region = temp_regions[key]
            for i in range(1,length):
                v = list_region[i]
                #print(v.x_coordinate, v.y_coordinate, v.z_coordinate, v.regionID)
                v.x_coordinate += i_step
                v.y_coordinate += j_step
                v.z_coordinate += k_step
                if v.x_coordinate>220 and v.y_coordinate>220 and v.z_coordinate>220:
                    print("r0 update")
                    print(v.x_coordinate, v.y_coordinate, v.z_coordinate)
        return [temp_M_region, temp_regions]


def outputregion(regions, shape):       # output segment regions
    # rs = len(regions)
    for key in regions:
        fname='emdr'+str(i)+'.mrc'
        mrc_new = mrcfile.new('mrcfilestest/emd0339/{}'.format(fname), overwrite=True)
        mrc_new.set_data(np.zeros(shape, dtype=np.float32))
        mrc_new.voxel_size = mrc.voxel_size
        for v in regions[key]:
            print(v.x_coordinate, v.y_coordinate, v.z_coordinate, v.density)
            mrc_new.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density
        #print(mrc_new.data)
        mrc_new.close()


# intialize()
# t1 = datetime.now()
# tArray= readData(img_matrix)
# t2 = datetime.now()
# delta = t2-t1
# df['read data']=delta
# print("read data: "+ str(delta))
# region_to_lm, regions = gen_regions(img_matrix, tArray)
# M_regions, regions = merge_region(region_to_lm, regions)
# df.to_csv('mrcfilestest/emd0339/output.csv')
# outputregion(regions, shape)
# mrc.close()
#
intialize()
cube_id = -1
regions = dict()
region_to_lm = []
for k in range(0, nz-unit, unit):
    for j in range(0, ny-unit, unit):
        for i in range(0, nx-unit, unit):
            cube_id += 1
            i_boundary, j_boundary, k_boundary = min(i+unit,nx), min(j+unit,ny), min(k+unit,nz)
            temp_matrix = img_matrix[i:i_boundary, j:j_boundary, k:k_boundary]
            #print(temp_matrix)
            threshold = temp_matrix.mean()
            #print(threshold)
            tArray = readData(temp_matrix,threshold,cube_id)
            temp_region_to_lm, temp_regions = gen_regions(temp_matrix,tArray,global_regionid)
            # for v in temp_region_to_lm:
            #     print(v.regionID)
            # rlen = len(temp_regions)
            # for i in range(rlen):
            #     for v in temp_regions[i]:
            #         print(v.regionID)
            temp_M_region, temp_regions = merge_region(temp_region_to_lm, temp_regions, temp_matrix)
            #print(temp_regions)
            temp_M_region, temp_regions = reverse_coordinate(temp_M_region, temp_regions, i, j, k)
            regions.update(temp_regions)
            region_to_lm.extend(temp_M_region)

print("length of local max: " + str(len(region_to_lm)))
print("number of regions:" + str(len(regions)))

region_to_lm, regions = merge_region(region_to_lm, regions, img_matrix)
outputregion(regions, shape)
print("done")

# tArray[0]=tArray[0]+(j,)
# print(tArray)





# def splitblock(matrix, unit):       # split matrix into smaller 3d blocks
#     width = range(0, matrix.shape[0],unit)
#     height = range(0, matrix.shape[1],unit)
#     depth = range(0, matrix.shape[2], unit)
#     subs = []
#     shared_v = []
#     for d in depth:
#         for h in height:
#             for w in width:
#                 subs.append(matrix[d:d+unit, h:h+unit, w:w+unit])
#                 if (d != 0 and d % unit == 0) or (h !=0 and h % unit == 0) or (w !=0 and w % unit == 0):
#                     v = Voxel(w, h, d, density = img_matrix[w, h, d])
#                     shared_v.append(v)
#     return [subs, shared_v]











# # assume all mrc files with coordinates starting from 0;
# mrc = mrcfile.open('mrcfiles/emd0339.mrc', mode='r+')
# # create a new mrc file writable for data change
# img_matrix = np.copy(mrc.data)
#
# # set threshold value as the mean density
# threshold = mrc.data.mean()
# nx = mrc.header.nx
# ny = mrc.header.ny
# nz = mrc.header.nz
# unit= int(math.sqrt(nx))
# splited = splitblock(img_matrix,unit)
# print(splited)









# def final_merge(combined)

# tArray = []
# regionx = []
# Read the data from mrc.data to voxel object and save in tArray



# for z in range(0, nz):
#     print('z ', z)
#     for y in range(0, ny):
#         for x in range(0, nx):
#             density = img_matrix[x, y, z]
#             if density < threshold:
#                 vx = Voxel(x, y, z, density)
#                 regionx.append(vx)
#                 density = 0
#             v = Voxel(x, y, z, density)
#             img_matrix[x, y, z] = density
#             tArray.append(v)

# size = img_matrix.size
# print(size)
# splited = splitblock(img_matrix,unit)
# print(splited)
# sort array decreasing order based on density; stable
# tSortedArray = sorted(tArray, key=lambda voxel: voxel.density, reverse=True)

# # record region number
# regionNum = -1
# # record region id:voxel with local maximum in corresponding region
# region_to_lm: List[Voxel] = []
#
# # create a dictionary with key as reigion id and value as list of voxels
# regions = dict()
# # find regions for each voxel
# for i in range(0, size):
#     v = tSortedArray[i]
#     # only take voxels with density value larger than threshold
#     if v.density > 0:
#         # dictionary holds regionID and number of it, regionID:number of voxels
#         regionRecord = dict()
#         # Get the list of neighbors of voxel
#         n = neighbors(img_matrix, v.x_coordinate, v.y_coordinate, v.z_coordinate)
#         # in the neighbors, check the region id of each voxel
#         for pos in n:
#             # calculate the index in before sorted array to find the region id
#             index = pos[0] + ny*pos[1] + ny*nz*pos[2]
#             # get the region id of neighbor
#             rId = tArray[index].regionID
#             # if neighbors' region id exists
#             if rId != -1:
#                 # check if this region id in the dictionary, if not assign region id with value 1, if yes then
#                 # increase the value of region id
#                 if rId in regionRecord:
#                     regionRecord[rId] += 1
#                 else:
#                     regionRecord[rId] = 1
#         # check region dictionary, if it is empty then assign new region id
#         if len(regionRecord) == 0:
#             # increase region id record, and assign the voxel with this new region id
#             regionNum += 1
#             v.regionID = regionNum
#             # rlm = Regionlm(regionNum, v)
#             region_to_lm.insert(regionNum, v)
#             regions[regionNum] = []
#             regions[regionNum].append(v)
#         # if region dictionary has item, assign the voxel with existing region id in region record
#         # the id will be the key with maximum value in the region dictionary
#         else:
#             r = max(regionRecord, key=regionRecord.get)
#             v.regionID = r
#             regions[r].append(v)
#     # if voxel density = 0, then exit for loop, as the sorted array in decrease order, the density voxel afterwards
#     # are all zero, no need to check
#     else:
#         break

# img_matrix = gaussian_filter(img_matrix, sigma=1, mode="constant", cval=0.0, truncate=1.0)
# emdrg1 = mrcfile.new('mrcfiles/emd4297g1.mrc', overwrite=True)
# emdrg1.set_data(img_matrix)
# emdrg1.close()

# mrc.close()
# print(len(region_to_lm))

#set the number of smoothing steps
# tStep = 10
# count = 0
# # reorder M in increasing order based on density,
# region_to_lm.reverse()  # mSorted=sorted(region_to_lm, key=lambda voxel: voxel.density)
# # calculate the gaussian filtered image, modify sigma, truncate value to optimize
# while count < tStep:
#     img_matrix = gaussian_filter(img_matrix, sigma=1, mode="constant", cval=0.0, truncate=4.0)
#     print("smoothed")
#     p = len(region_to_lm)
#     q = int((1 + p)/2)
#     # update density after smoothing
#     for i in range(0, p):
#         xc = region_to_lm[i].x_coordinate
#         yc = region_to_lm[i].y_coordinate
#         zc = region_to_lm[i].z_coordinate
#         region_to_lm[i].density = img_matrix[xc, yc, zc]
#     print("update density")
#     for i in range(0, q):
#         gra = dict()
#         for j in range(i+1, p):
#             # dictionary to keep all gradients for mi
#             gra[j] = gradient(region_to_lm[i], region_to_lm[j])
#         # find steepest ascent, the max value of g
#         # print(gra)
#         k = max(gra, key=gra.get)
#         # print(k)
#         k = region_to_lm[k].regionID
#         # merge region i to region k
#         temp = region_to_lm[i].regionID
#         for v in regions[temp]:
#             v.regionID = k
#         regions[k].extend(regions[temp])
#         # remove checked regions
#         regions = delete_key(temp, regions)
#     del region_to_lm[0:q]
#     count += 1
#     print(count)
#
# print(len(regions))
# rs=len(regions)


# tArray[0]=tArray[0]+(j,)
# print(tArray)