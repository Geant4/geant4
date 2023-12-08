import sys
import struct
import math
import numpy as np

# lists for super resolution (vol_work) and regular resolution (vol_result)
vol_work = []
vol_result = []

############################################################################
# size image 128x128 -> 500 microm -> resol x,y : 3.90625 microm      #
# nb slices 128 -> 500 microm -> resol z : 3.90625 microm            #
############################################################################

sizex = 128
sizey = 128
sizez = 128
resolx = 3.90625
resoly = 3.90625
resolz = 3.90625
superres = 8  # factor of super resolution

# sphere 1 (radius) - outer sphere
r1 = 196

# sphere 2 (radius) - inner sphere
r2 = 171

# density values for sphere 1 and sphere 2
type = 2  # type for constructing a STIM or PIXE phantom
# type = 1, STIM phantom, density value for STIM in 0.01 g/cm3
# type =2, PIXE phantom, density value for PIXE in 0.000001 g/cm3, namely microgram/cm3
value1 = 0
value2 = 0

if type == 1:
    value1 = 108
    value2 = 0
elif type == 2:
    value1 = 54000
    value2 = 0

# center of two spheres
x0 = 0.0
y0 = 0.0
z0 = -1.953125

# size in super resolution by voxel
super_sizex = sizex * superres
super_sizey = sizey * superres
super_sizez = sizez * superres
center_shift = sizex / 2  # translation of half scan


def make_sphere_center(r, x, y, z, value):
    rsample = (r / resolx) * superres  # radius in super resolution by voxel
    xsample = (x / resolx) * superres
    ysample = (y / resoly) * superres
    zsample = (z / resolz) * superres
    center_shift_sample = center_shift * superres  # translation of the center for x (i) and y (j) axis

    print(rsample, xsample, ysample, zsample)
    number = 0

    for k in range(0, super_sizez):
        for j in range(0, super_sizey):
            for i in range(0, super_sizex):
                ii = i + 0.5
                jj = j + 0.5
                kk = k + 0.5
                res = pow((ii - xsample - center_shift_sample), 2) / (rsample * rsample) + \
                      pow((jj - ysample - center_shift_sample), 2) / (rsample * rsample) + \
                      pow((kk - zsample - center_shift_sample), 2) / (rsample * rsample)  # z-axis correction done
                # if the point (ii, jj, kk) is in the sphere, we attribute the voxel (i,j, k) value
                if (res <= 1.0):
                    vol_work[i + j * super_sizex + k * super_sizex * super_sizey] = value
                    number += 1

    print(number)


# initialisation for two tables vol_result (128*128*128), vol_work (128*superres)*(128*superres)*(128*superres)
def initialize():
    for k in range(0, sizez):
        for j in range(0, sizey):
            for i in range(0, sizex):
                vol_result.append(0.0)
    for k in range(0, super_sizez):
        for j in range(0, super_sizey):
            for i in range(0, super_sizex):
                vol_work.append(0.0)


# Calculate the vol_result based on vol_work
def undersample():
    x = 0
    y = 0
    z = 0
    for k in range(0, super_sizez, superres):
        for j in range(0, super_sizey, superres):
            for i in range(0, super_sizex, superres):
                # print ("***",i,j,k)
                total = 0
                for kk in range(0, superres):
                    for jj in range(0, superres):
                        for ii in range(0, superres):
                            total = total + vol_work[
                                i + ii + (j + jj) * super_sizex + (k + kk) * (super_sizex * super_sizey)]
                vol_result[x + y * sizex + z * (sizex * sizey)] = total / (superres * superres * superres)
                # print ("###",vol_result[x+y*sizex+z*(sizex*sizey)])
                x += 1
            y += 1
            x = 0
        z += 1
        y = 0


# save the total volume of super resolution
def save_whole_work(file):
    fd = open(file, "wb")
    for i in range(0, len(vol_work)):
        fd.write(struct.pack("f", vol_work[i]))
    fd.close()


# save one slice in super resolution vol_work
#  0<=slice<super_sizez (128*8=1024)
def save_workslice(file, slice):
    fd = open(file, "wb")
    for j in range(0, super_sizey):
        for i in range(0, super_sizex):
            fd.write(struct.pack("f", vol_work[i + j * super_sizex + slice * super_sizex * super_sizey]))
    fd.close()


# save the total result volume
def save_whole_result(file):
    fd = open(file, "wb")
    for i in range(0, len(vol_result)):
        fd.write(struct.pack("f", vol_result[i]))
    fd.close()


# save one slie in regular resolution vol_result
#  0<=slice<128
def save_slice(file, slice):
    fd = open(file, "wb")
    for j in range(0, sizey):
        for i in range(0, sizex):
            fd.write(struct.pack("f", vol_result[i + j * sizex + slice * sizex * sizey]))

    fd.close()


print("------intialisation------")
initialize()
print("------end intialisation------")

print("------begin sphere 1------")
make_sphere_center(r1, x0, y0, z0, value1)
print("------end sphere 1------")

print("------begin sphere 2------")
make_sphere_center(r2, x0, y0, z0, value2)
print("------end sphere 2------")

print("------begin undersample------")
undersample()
print("------end undersample------")

print("The length of vol_work: ", len(vol_work))
print("The length of vol_result: ", len(vol_result))

# print ("------begin save slice 63 ------")
# save_slice("./slice_63.dat",63)
# print ("------end save slice 63 ------")

save_whole_result("./vol_result.dat")
save_whole_work("./vol_work.dat")
print("------end save work------")
