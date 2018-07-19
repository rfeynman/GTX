#!/usr/bin/python

from math import *
import os
import sys 
from numpy import *
from scipy import *
from scipy.linalg import *
from numpy.linalg import *

fo = open("_b.agr", "w")

B  = 1
b1 = 0.999
b2 = 0

gap = 0.1
b1 = 2/gap
print b1

y =  0.005

for i in range(-500, 501) :
	z = i/100.
	f = b1*z + b2*(z**2-y**2)
	h = y*(b1 + 2 * b2 * z)
	b = 1 / (1 + exp(b1*z))
	by =  B * (1.+exp(f)*cos(h)) / ( 1 + 2*exp(f)*cos(h) + exp(2*f))
	bz = -B * (exp(f)*sin(h)) / ( 1 + 2*exp(f)*cos(h) + exp(2*f))
	print >> fo, z, by

print >> fo, "&"

for i in range(-500, 501) :
	z = i/100.
	f = b1*z + b2*(z**2-y**2)
	h = y*(b1 + 2 * b2 * z)
	b = 1 / (1 + exp(b1*z))
	by =  B * (1.+exp(f)*cos(h)) / ( 1 + 2*exp(f)*cos(h) + exp(2*f))
	bz = -B * (exp(f)*sin(h)) / ( 1 + 2*exp(f)*cos(h) + exp(2*f))
	print >> fo, z, bz


fo.close()
os.system("xmgrace -geometry 780x630+600 -fixed 510 380 -noask _b.agr &")
