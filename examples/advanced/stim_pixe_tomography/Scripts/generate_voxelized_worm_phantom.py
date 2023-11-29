import sys
import struct
import math

# tables for super resolution (vol_work) and regular resolution (vol_result)
vol_work = []
vol_result = []

############################################################################
# size image 128x128 -> 76.464 microm -> resol x,y : 0,597375 microm      #
# nb slices 128 -> 201,217 microm -> resol z : 1,57200781 microm            #
# slice of interest : 11 -> 18,07 microm                                      #
############################################################################

sizex = 128
sizey = 128
sizez = 128
double_sizez = 256  # for computation only
resolx = 0.597375
resoly = 0.597375
resolz = 1.57200781
superres = 8

slice = 11

# ellipsoide 1 (semi-axes, center, rotation z, value) - skin
a1 = 20.61
b1 = 21.42
c1 = 187.82
x1 = 0.0
y1 = 0.0
z1 = 0.0
rz1 = 0.0
value1 = 49.73

# ellipsoide 2 (semi-axes, center, rotation z, value) - body
a2 = 18.61
b2 = 19.01
c2 = 186.64
x2 = -0.39
y2 = 0.0
z2 = 0.0
rz2 = 0.0
value2 = 40.85

# ellipsoide 3 (semi-axes, center, rotation z, value) - core 1
a3 = 1.95
b3 = 3.23
c3 = 4.32
x3 = 1.97
y3 = -7.09
z3 = 18.07
rz3 = 0.0
value3 = 66.26

# ellipsoide 4 (semi-axes, center, rotation z, value) - core 2
a4 = 2.08
b4 = 2.46
c4 = 4.32
x4 = 8.27
y4 = -3.15
z4 = 18.07
rz4 = 0.0
value4 = 60.23

# ellipsoide 5 (semi-axes, center, rotation z, value) - intestine
a5 = 3.67
b5 = 16.48
c5 = 28.68
x5 = 1.25
y5 = 0.61
z5 = 0
rz5 = -58.99
value5 = 54.06

# ellipsoide 6 ((emi-axes, center, rotation z, value) - region titane
a6 = 1.62
b6 = 1.95
c6 = 1.62
x6 = 6.25
y6 = 3.61
z6 = 18.07
rz6 = 0.0
value6 = 75.14

# size in super resolution by voxel
super_sizex = sizex * superres
super_sizey = sizey * superres
super_sizez = double_sizez * superres  # we will save only half of the z
center_shift = sizex / 2  # for x and y axis
center_shift_z = double_sizez / 2


def make_ellipse_center(a, b, c, x, y, z, rz, value):
    asample = (a / resolx) * superres
    bsample = (b / resoly) * superres
    csample = (c / resolz) * superres
    # if (rz != 0.0):
    rzsample = math.radians(rz)  # angle in radians
    xsample = (x / resolx) * superres
    ysample = (y / resoly) * superres
    zsample = (z / resolz) * superres
    center_shift_sample = center_shift * superres  # translation of center for x and y axis
    center_shift_z_sample = center_shift_z * superres  # translation of center for z axis

    print(asample, bsample, csample, xsample, ysample, zsample)
    number = 0  ## debug

    # a loop in axes of voxels
    for k in range(0, super_sizez):
        for j in range(0, super_sizey):
            for i in range(0, super_sizex):
                ii = i + 0.5
                jj = j + 0.5
                kk = k + 0.5
                res = pow(((ii - xsample - center_shift_sample) * math.cos(rzsample) + (
                            jj - ysample - center_shift_sample) * math.sin(rzsample)), 2) / (asample * asample) + \
                      pow(((ii - xsample - center_shift_sample) * math.sin(rzsample) - (
                                  jj - ysample - center_shift_sample) * math.cos(rzsample)), 2) / (bsample * bsample) + \
                      ((kk - zsample - center_shift_z_sample) * (kk - zsample - center_shift_z_sample)) / (
                                  csample * csample)  # z-axis correction done
                # print(res)
                # if the voxel belongs to the ellipsoide, set the value
                if (res <= 1.0):
                    vol_work[i + j * super_sizex + k * super_sizex * super_sizey] = value
                    number += 1

    print(number)


# initialisation
def initialize():
    for k in range(0, double_sizez):
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
                # on se place sur v1 et on recupere les valeurs de densite des 8 voxels du voisinage qui vont correspondre a 1 voxel de l'image finale
                # v1 = vol_work[i+j*super_sizex+k*(super_sizex*super_sizey)]
                # v2 = vol_work[i+j*super_sizex+k*(super_sizex*super_sizey) + 1]
                # v3 = vol_work[i+j*super_sizex+k*(super_sizex*super_sizey) + super_sizex]
                # v4 = vol_work[i+j*super_sizex+k*(super_sizex*super_sizey) + super_sizex + 1]
                # v5 = vol_work[i+j*super_sizex+k*(super_sizex*super_sizey) + super_sizex*super_sizey]
                # v6 = vol_work[i+j*super_sizex+k*(super_sizex*super_sizey) + 1 + super_sizex*super_sizey]
                # v7 = vol_work[i+j*super_sizex+k*(super_sizex*super_sizey) + super_sizex + super_sizex*super_sizey]
                # v8 = vol_work[i+j*super_sizex+k*(super_sizex*super_sizey) + super_sizex + 1 + super_sizex*super_sizey]
                # vol_result[x+y*sizex+z*(sizex*sizey)] = (v1+v2+v3+v4+v5+v6+v7+v8)/8.0
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


# Save the total volume (including the negative parts of the ellipsoids)
def save_whole_work(file):
    fd = open(file, "wb")
    for i in range(0, len(vol_work)):
        fd.write(struct.pack("f", vol_work[i]))


# Save half the volume (including only the positive parts of the ellipsoids)
def save_half_work(file):
    fd = open(file, "wb")
    for i in range(int(len(vol_work) / 2), len(vol_work)):
        fd.write(struct.pack("f", vol_work[i]))


# save one slice in super resolution vol_work
#  0<=slice<super_sizez
def save_workslice(file, slice):
    fd = open(file, "wb")
    for j in range(0, super_sizey):
        for i in range(0, super_sizex):
            fd.write(struct.pack("f", vol_work[i + j * super_sizex + slice * super_sizex * super_sizey]))


# Save the total result volume (including the negative parts of the ellipsoids)
def save_whole_result(file):
    fd = open(file, "wb")
    for i in range(0, len(vol_result)):
        fd.write(struct.pack("f", vol_result[i]))


# Save half of the result volume (including only the positive parts of the ellipsoids)
def save_half_result(file):
    fd = open(file, "wb")
    for i in range(int(len(vol_result) / 2), len(vol_result)):
        fd.write(struct.pack("f", vol_result[i]))


# save one slie in regular resolution vol_result
#  0<=slice<128
def save_slice(file, slice):
    fd = open(file, "wb")
    for j in range(0, sizey):
        for i in range(0, sizex):
            fd.write(struct.pack("f", vol_result[i + j * sizex + slice * sizex * sizey]))


print("------intialisation------")
initialize()
print("------end intialisation------")

print("------begin ellipse 1------------")
make_ellipse_center(a1, b1, c1, x1, y1, z1, rz1, value1)
print("------end ellipse 1--------------")

print("------begin ellipse 2------------")
make_ellipse_center(a2, b2, c2, x2, y2, z2, rz2, value2)
print("------end ellipse 2--------------")

print("------begin ellipse 3------------")
make_ellipse_center(a3, b3, c3, x3, y3, z3, rz3, value3)
print("------end ellipse 3--------------")

print("------begin ellipse 4------------")
make_ellipse_center(a4, b4, c4, x4, y4, z4, rz4, value4)
print("------end ellipse 4--------------")

print("------begin ellipse 5------------")
make_ellipse_center(a5, b5, c5, x5, y5, z5, rz5, value5)
print("------end ellipse 5--------------")

print("------begin ellipse 6------------")
make_ellipse_center(a6, b6, c6, x6, y6, z6, rz6, value6)
print("------end ellipse 6--------------")

print("------begin undersample------")
undersample()
print("------end undersample------")

print("The length of vol_work: ", len(vol_work))
print("The length of vol_result: ", len(vol_result))

# print ("------begin save slice 139 ------")
save_slice("./slice_139.dat", 139)
# print ("------end save slice 139 ------")


save_half_result("./vol_result.dat")
# save_whole_work("./vol_work.dat")
print("------end save work------")
