import numpy as np
import mrcfile
from pathlib import Path
import os

mrc = mrcfile.open('emd4297.mrc')
print(mrc.header)
# Print the whole header using built in command
mrc.print_header
print('\n\n\n')
#show the physical dimensions of the volume
print("nx = " + str(mrc.header.nx))
print("ny = " + str(mrc.header.ny))
print("nz = " + str(mrc.header.nz))

# show voxel size
print("voxel size = " + str(mrc.voxel_size))

# show data shape
print("data shape = " + str(mrc.data.shape))


# show volume info
print("Is a volume = " + str(mrc.is_volume()))

# show one data element
print("Data 10,10,10" + str(mrc.data[120,120,120]) )
print(str(mrc.header.__mod__))
for z in range (0,mrc.header.nz):
    for y in range( 0, mrc.header.ny):
        for x in range(0, mrc.header.nx):
            print ("Data " + str(x) + ", " + str(y) + ", " + str(z) +" = " + str(mrc.data[x,y,z]))