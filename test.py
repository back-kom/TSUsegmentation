# class Voxel(object):
#     def __init__(self, x_coordinate, y_coordinate, z_coordinate, density, region_id=-1, cube_id=None):
#         self.x_coordinate = x_coordinate
#         self.y_coordinate = y_coordinate
#         self.z_coordinate = z_coordinate
#         self.density = density
#         self.regionID = region_id
#         self.cubeID = cube_id
#
#     def get_coordiate(self):
#         return (self.x_coordinate, self.y_coordinate, self.z_coordinate)
#
#
# v = Voxel(1,2,3,0.5)
# print(v.get_coordiate())


import mrcfile
import numpy as np
import math


fname1 = 'mrcfiles/emd4297.mrc'
# fname3 = 'SourceFiles/emd4297segger6steps/emd_4297_region_1741.mrc'
# fname2 = 'mrcfiles/emd0414.mrc'
# fname4 = 'SourceFiles/emd4297segger6steps/emd_4297_region_1743.mrc'
# fname5 = 'SourceFiles/emd4297segger6steps/emd_4297_region_1750.mrc'
# fname6 = 'SourceFiles/emd4297segger6steps/emd_4297_region_1752.mrc'
# fname7 = 'SourceFiles/emd4297segger6steps/emd_4297_region_1756.mrc'
# fname8 = 'SourceFiles/emd4297segger6steps/emd_4297_region_1757.mrc'
# fname9 = 'SourceFiles/emd4297segger6steps/emd_4297_region_1758.mrc'
# fname10 = 'SourceFiles/emd4297segger6steps/emd_4297_region_1759.mrc'
#
#
mrc1 = mrcfile.open(fname1, mode='r+')
# mrc2 = mrcfile.open(fname2, mode='r+')
# mrc3 = mrcfile.open(fname3, mode='r+')
# mrc4 = mrcfile.open(fname4, mode='r+')
# mrc5 = mrcfile.open(fname5, mode='r+')
# mrc6 = mrcfile.open(fname6, mode='r+')
# mrc7 = mrcfile.open(fname7, mode='r+')
# mrc8 = mrcfile.open(fname8, mode='r+')
# mrc9 = mrcfile.open(fname9, mode='r+')
# mrc10 = mrcfile.open(fname10, mode='r+')
#
# print(mrc1.header)
# print(mrc2.header)
# print(mrc3.header)
# print(mrc4.header)
# print(mrc5.header)
# print(mrc6.header)
# print(mrc7.header)
# print(mrc8.header)
# print(mrc9.header)
# print(mrc10.header)

img_matrix = np.copy(mrc1.data)
threshold = img_matrix.mean()
nx = mrc1.header.nx
print("nx: " + str(nx))
ny = mrc1.header.ny
nz = mrc1.header.nz
shape = (nx, ny, nz)
unit = 13
print("n: " + str(unit))

for k in range(0, nz, unit):
    for j in range(0, ny, unit):
        for i in range(0, nx, unit):
            i_boundary, j_boundary, k_boundary = min(i+unit,nx), min(j+unit,ny), min(k+unit,nz)
            temp_matrix = img_matrix[i:i_boundary, j:j_boundary, k:k_boundary]
            print("x of cube: " + str(len(temp_matrix)))
            print("y of cube: " + str(len(temp_matrix[0])))
            print("z of cube: " + str(len(temp_matrix[0][0])))
