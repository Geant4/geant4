#!/usr/bin/env python
# Python program explaining  
# load() function  
  
import numpy as geek 
import matplotlib.pyplot as plt 
  
b = geek.load('hits.npy') 
  
print("size of b is:",b.size)
print("shape of b is:",b.shape)
print("ndim of b is:",b.ndim)
print(b)
print(b[0,0,0])
print(b[0,3,0])
print("b is printed from hits.npy")
x = b[:,0,0]
y = b[:,0,1]
z = b[:,0,2]
lamb = b[:,2,3]
print(x.shape)
print(y.shape)
print(lamb.shape)

#axs = plt.subplots(2, 2, figsize=(5, 5))
plt.subplot(131)
plt.ylabel('#Photons')
plt.xlabel('x-position/mm')
plt.hist(x,50,(-200,200.))
plt.subplot(132)
plt.xlabel('x-position/mm')
plt.ylabel('y-position/mm')
plt.scatter(x, y)
plt.subplot(133)
plt.xlabel('lambda/nm')
plt.ylabel('#Photons')
plt.hist(lamb,50,(50,800.))
#plt.hist(z)
plt.show()
