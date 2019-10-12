import os
import mrcfile
import numpy as np
import pandas as pd


This_Folder = os.path.dirname(os.path.abspath(__file__))
# data structure holds voxel information
# class Voxel(object):
#     def __init__(self, x_coordinate, y_coordinate, z_coordinate, density, region_id=-1):
#         self.x_coordinate = x_coordinate
#         self.y_coordinate = y_coordinate
#         self.z_coordinate = z_coordinate
#         self.density = density


# get the set of voxels from the test region
def get_region_set(fname, threshold):
    # initialize the set
    regionset = set()
    fname = os.path.join(This_Folder, fname)
    tempregion_mrc = mrcfile.open(fname, mode='r+')
    temp_img = np.copy(tempregion_mrc.data)
    row = temp_img.shape[0]
    col = temp_img.shape[1]
    dep = temp_img.shape[2]
    # print(row, col, dep)
    # read data
    for z in range(0, dep):
        for y in range(0, col):
            for x in range(0, row):
                density = temp_img[x, y, z]
                # check density, if >= threshold, then update to 1, and add voxel to the set
                if density >= threshold:
                    density = 1
                    # v = Voxel(x, y, z, density)
                    # print(x,y,z,density)
                    regionset.add((x,y,z,density))
    return regionset


# get dictionary of templates with pairs of "filename:set of voxels"
def get_regiondict(dname, threshold):
    template_regionsetdict = {}
    dname = os.path.join(This_Folder, dname)
    for file in os.listdir(dname):
        if file.endswith(".mrc"):
            fname = file
            file = dname + '/' + file
            templateregionset = get_region_set(file, threshold)
            template_regionsetdict.update({fname:templateregionset})
    return template_regionsetdict


output_path = "Output"
testdir = input("choose the test region directory:")
testdir = "mrcfilestest/"+ testdir
tempdir = input("choose the templates directory:")
tempdir = "mrcfilestest/"+ tempdir
file = input("choose the original mrc file:")
file = "mrcfiles/" + file
mrc = mrcfile.open(file)
# threshold = float(input("input the threshold value:"))
threshold = mrc.data.mean()
testregiondict = get_regiondict(testdir, threshold)

templateregiondict = get_regiondict(tempdir, threshold)
df = pd.DataFrame()

for testname, testregionset in testregiondict.items():
    testname = "dc_"+ testname
    column_lst = []
    output_list = []
    for templatename, templateregionset in templateregiondict.items():
        templatename = "non_dc_" + templatename
        column_lst.append(templatename)
        intersection = testregionset.intersection(templateregionset)
        # get the number of elements in test region set and template region sets
        num_intersection = len(intersection)
        num_templateset = len(templateregionset)
        print(num_intersection, num_templateset)
        percentage = format((num_intersection / num_templateset),'.2%')
        output_list.append(percentage)
        # print(percentage)
    df[testname]=output_list
df = df.transpose()
df = pd.DataFrame(df.values, df.index, column_lst)
print(df)


df.to_csv(os.path.join(output_path, r'emd0414smoothtwicesparam23_6_overlap.csv'))


