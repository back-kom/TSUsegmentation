import mrcfile
import numpy as np
import math
from datetime import datetime
from scipy.ndimage import gaussian_filter


# data structure holds voxel information
class Voxel(object):
    def __init__(self, x, y, z, density, region_id=-1, nlist=None):
        self.x_coordinate = x
        self.y_coordinate = y
        self.z_coordinate = z
        self.density = density
        self.regionID = region_id
        self.nlist=nlist

# class Voxel(object):
#     def __init__(self, x, y, z, density, region_id= -1):
#         self.x_coordinate = x
#         self.y_coordinate = y
#         self.z_coordinate = z
#         self.density = density
#         self.regionID = region_id
#     def updaterId(self, x, y, z, rId):
#         if self.x_coordinate == x and self.y_coordinate == y and self.z_coordinate == z:
#             self.regionID = rId


def three_d_array(x, y, z):
    return [[[None for k in range(z)] for j in range(y)]for i in range(x)]


class DNode(object):
    def __init__(self, data = None, prev = None, next = None, ):
        self.data = data
        self.prev = prev
        self.next = next


class DLinkedList(object):
    def __init__(self):
        self.head = None
        self.tail = None
        self.size = 0

    def AppendToHead(self, data):
        new_node = DNode(data=data)
        if self.head:
            new_node.next = self.head
            self.head.prev = new_node
            self.head = new_node
        else:
            self.head = new_node
            self.tail = new_node
        self.size += 1

    def AppendToTail(self, data):
        new_node = DNode(data=data)
        if self.tail:
            new_node.prev = self.tail
            self.tail.next = new_node
            self.tail = new_node
        else:
            self.tail = new_node
            self.head = new_node
        self.size +=1

    def Search(self, key):
        current = self.head
        while current and current.data != key:
            current = current.next
        return current

    def RemoveFromHead(self):
        x = self.head
        if self.head:
            if self.head == self.tail:
                self.head = None
                self.tail = None
            else:
                self.head = self.head.next
                self.head.prev = None
        self.size -= 1
        return x

    def RemoveFromTail(self):
        x = self.tail
        if self.tail:
            if self.head == self.tail:
                self.head = None
                self.tail = None
            else:
                self.tail = self.tail.prev
                self.tail.next = None
        self.size -= 1
        return x


# data structure for regions
class Tree(object):
    def __init__(self, root):
        self.root = root
        self.children = []
        self.size = 1

    def add_child(self, node):
        # assert isinstance(node, Tree)
        self.children.append(node)
        self.size += 1

    def add_children(self, children):
        if children is not None:
            for child in children:
                self.add_child(child)

    def get_size(self):
        return self.size

    def get_root(self):
        return self.root

# class ListNeighbor(object):
#     def __init__(self):


def readData(matrix):       # Read the data from mrc.data to voxel object and save in vList
    vList = []
    regionx = []
    row = matrix.shape[0]
    col = matrix.shape[1]
    dep = matrix.shape[2]
    for z in range(dep):
        for y in range(col):
            for x in range(row):
                density = matrix[x, y, z]
                v = Voxel(x, y, z, density)
                if density < threshold:
                    regionx.append(v)
                else:
                    vList.append(v)
    # return [vList, regionx]
    return vList

# initialize program and ask user to input filename
def initialize():        # initialize program
    global mrc, img_matrix, nx, ny, nz, size, img_3d #, unit
    fname = input("choose mrc file:")
    mrc = mrcfile.open(fname, mode='r+')
    img_matrix = np.copy(mrc.data)
    nx = mrc.header.nx
    ny = mrc.header.ny
    nz = mrc.header.nz
    size = img_matrix.size
    img_3d = three_d_array(None,nx,ny,nz)
    # checkpoint
    print("number of total voxels: %d" % (size))
    #unit = int(math.sqrt(nx))


def smoothing(matrix, sig=1, cv=0.0, trunc=4.0):    #gaussian filter
    return gaussian_filter(matrix, sigma=sig, mode="constant", cval=cv, truncate=trunc)


def neighbors(matrix, voxel):       # return list of neighbor's coordinates
    # initialize list of neighbors
    neighbor= []
    # get x boundary
    row = len(matrix)
    # get y boundary value
    col = len(matrix[0])
    # get z boundary value
    dep = len(matrix[0][0])
    # get x, y, z coordinates
    x, y, z = voxel.x_coordinate, voxel.y_coordinate, voxel.z_coordinate
    # loop to find neighbors coordinates, index must greater or equal to 0, and less or equal to the boundary value
    for k in range(max(0, z-1), min(dep, z+2)):
        for j in range(max(0, y-1), min(col, y+2)):
            for i in range(max(0, x-1), min(row, x+2)):
                # check if density is less than threshold
                if matrix[i, j, k] >= threshold:
                    # exclude itself
                    if(i, j, k) != (x, y, z):
                        neighbor.append((i, j, k))
                        
    return neighbor


def getRegions(matrix, vList):
    mregion = []
    regionNum = -1
    t1 = datetime.now()
    vSortedList = sorted(vList, key=lambda voxel: voxel.density, reverse=True)
    t2 = datetime.now()
    delta = t2 - t1
    print("time cost of sort is : %f" % delta.total_seconds())
    print("sorted done, and the number of voxels above threshold is %d" % (len(vSortedList)))
    c = 0
    for v in vSortedList:
        regionRecord = dict()
        if v.density >= threshold:
            c += 1
            vi = v.x_coordinate + nx * v.y_coordinate + nx * ny * v.z_coordinate
            nb = neighbors(matrix, v)
            for pos in nb:
                index = pos[0] + nx*pos[1] + nx*ny*pos[2]
                rId = vList[index].regionID
                if rId != -1:
                    if rId in regionRecord:
                        regionRecord[rId] += 1
                    else:
                        regionRecord[rId] = 1

            if len(regionRecord) == 0:
                regionNum += 1
                v.regionID = regionNum
                vList[vi].regionID = regionNum
                tree = Tree(root = v)
                mregion.insert(regionNum, tree)
            else:
                r = max(regionRecord, key = regionRecord.get)
                v.regionID = r
                vList[vi].regionID = r
                mregion[r].add_child(v)
        else:
            break

    print("number of voxels above threshold: %d" % c)
    return mregion


def gradient(mi, mj):       # calculate the gradient
    distance = pow((mi.x_coordinate - mj.x_coordinate), 2) + pow((mi.y_coordinate - mj.y_coordinate), 2) + pow((mi.z_coordinate - mj.z_coordinate), 2)
    return (mj.density - mi.density)/distance

t1 = datetime.now()
initialize()
t2 = datetime.now()
delta = t2 - t1
print("time cost of initialize is : %f" % delta.total_seconds())
img_matrix = smoothing(img_matrix)
print(img_matrix.mean())
img_matrix = smoothing(img_matrix)
threshold = img_matrix.mean()
print(threshold)
vList = readData(img_matrix)
t1 = datetime.now()
mregion = getRegions(img_matrix,vList)
t2 = datetime.now()
delta = t2 - t1
print("time cost of initial M0 is : %f" % delta.total_seconds())
print("number of regions at first before merge: %d" % (len(mregion)))

# tStep = 14
# count = 0
# mregion.reverse()
#
# while count < tStep:
#     t1 = datetime.now()
#     img_matrix = smoothing(img_matrix)
#     print("smoothed %d times" % (count+1))
#     p = len(mregion)
#     q = int((1+p)/2)
#     for t in mregion:
#         v = t.root
#         xc = v.x_coordinate
#         yc = v.y_coordinate
#         zc = v.z_coordinate
#         v.density = img_matrix[xc, yc, zc]
#     print("update density %d times" % (count+1))
#     for i in range(q+1):
#         gra = dict()
#         for j in mregion[1:-1]:
#             j = mregion.index(j)
#             gra[j]=gradient(mregion[0].root, mregion[j].root)
#         g = max(gra, key=gra.get)
#         k = mregion[g].root.regionID
#         mregion[k].add_child(mregion[0].root)
#         mregion[k].children.extend(mregion[0].children)
#         mregion.pop(0)
#     count += 1
#     t2 = datetime.now()
#     delta = t2 - t1
#     print("time cost of merge is : %f" % delta.total_seconds())
#     print("merged %d times" % (count))
#
# rs = len(mregion)
# print("number of regions: %d" % rs)
# shape = (nx, ny, nz)
#
# for i in range(0,rs):
#     fname='emdr'+str(i)+'.mrc'
#     mrc_new = mrcfile.new('mrcfilestest/{}'.format(fname), overwrite=True)
#     mrc_new.set_data(np.zeros(shape, dtype=np.float32))
#     mrc_new.voxel_size = mrc.voxel_size
#     t = mregion[i]
#     r = t.root
#     childlist = t.children
#     mrc_new.data[t.x_coordinate, t.y_coordinate, t.z_coordinate] = t.density
#     for v in childlist:
#         mrc_new.data[v.x_coordinate, v.y_coordinate, v.z_coordinate] = v.density
#     mrc_new.close()
#
# mrc.close()





