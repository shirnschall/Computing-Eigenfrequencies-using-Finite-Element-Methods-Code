from materials import *
from geometry import *
from evsolver import *
from ngsolve.la import EigenValues_Preconditioner
import netgen.gui
import netgen.geom2d as geom2d
from netgen.csg import *
from netgen.geom2d import unit_square
from ngsolve import *
import math
import scipy.linalg
from scipy import random
from netgen.NgOCC import *
from time import sleep

import matplotlib.pyplot as plt

evsolver=EVSolver()

#Geometry(material, geometryPath, maxh, order, dirichlet, shift, low_order_space)
tuningfork=Geometry(steel, "../3d-models/fork/tuning-fork.stp", 2, 2, [1,2,3], 0, True)

tuningfork.generateMesh()
Draw(tuningfork.mesh)

#pre(type, inverse)
tuningfork.pre("multigrid", "sparsecholesky")
tuningfork.assemble()


ev,evec,res = evsolver.lobpcg(tuningfork, 10, 500, 10e-6)




input("Press any key to continue..")