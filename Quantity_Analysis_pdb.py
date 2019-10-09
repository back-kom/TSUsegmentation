import os
import mrcfile
import numpy as np
import pandas as pd
from Bio.PDB.PDBParser import PDBParser


def initialize():
    global This_Folder
    This_Folder = os.path.dirname(os.path.abspath(__file__))


def pdb_ca_coord_dict():
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("test", "pdbfiles/emd4297_6fq5.pdb")
    ca_dict = dict()
    for model in structure.get_list():
        for chain in model.get_list():
            c_id = chain.get_id()
            # print(c_id)
            coord_set = set()
            for residue in chain.get_list():
                if residue.has_id("CA"):
                    ca = residue["CA"]
                    # if ca.get_bfactor() > 50.0:
                    coordinate = ca.get_coord()
                    x_co, y_co, z_co = int(round(coordinate[0])), int(round(coordinate[1])), int(round(coordinate[2]))
                    coordinate = (x_co, y_co, z_co)
                    coord_set.add(coordinate)
            ca_dict[c_id]=coord_set
    return ca_dict


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
                    # density = 1
                    # v = Voxel(x, y, z, density)
                    # print(x,y,z,density)
                    regionset.add((x,y,z))
    return regionset


def get_regiondict(dname, threshold):
    template_regionsetdict = dict()
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
threshold = float(input("input the threshold value:"))
testregiondict = get_regiondict(testdir, threshold)

templateregiondict = get_regiondict(tempdir, threshold)
df = pd.DataFrame()



