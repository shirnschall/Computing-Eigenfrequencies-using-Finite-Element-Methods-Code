from materials import *

import netgen.gui
import netgen.geom2d as geom2d
from netgen.csg import *
from netgen.geom2d import unit_square
from ngsolve import *
import math
import scipy.linalg
from scipy import random
from netgen.NgOCC import *

class Geometry():
    material=steel
    geometry="" #3d model
    maxh=1
    order=0
    dirichlet=0
    low_order_space=True
    mesh=0
    fes=0
    u=0   #TrialFunction
    v=0   #TestFunction
    A=0   #StiffnessMatrix
    M=0   #MassMatrix
    shift=0
    precond=False   #preconditioner
    def __init__(self, material, geometryPath, maxh, order, dirichlet, shift, low_order_space):
        self.shift=shift
        self.maxh=maxh
        self.order=order
        self.dirichlet=dirichlet
        self.low_order_space=low_order_space
        if(geometryPath!=""):
            self.material=material
            self.geometry=OCCGeometry(geometryPath)


    def generateMesh(self):
        self.mesh = Mesh(self.geometry.GenerateMesh(maxh=self.maxh))
        self.mesh.Curve(2)
        self.fes = VectorH1(self.mesh, order=self.order, dirichlet=self.dirichlet, low_order_space=self.low_order_space)
        self.u = self.fes.TrialFunction()
        self.v = self.fes.TestFunction()
        # StiffnessMatrix depending on material type:
        if (self.material.type == "linear"):
            self.A = BilinearForm(self.fes)
            self.A += 2 * self.material.mu * InnerProduct(1 / 2 * (grad(self.u) + grad(self.u).trans),
                                                     1 / 2 * (grad(self.v) + grad(self.v).trans)) * dx
            self.A += self.material.lam * div(self.u) * div(self.v) * dx
            self.A += self.shift * self.material.rho * self.u * self.v * dx
        # else if(1==1):
        # other material types

        self.M = BilinearForm(self.fes)
        self.M += self.material.rho * self.u * self.v * dx

    def pre(self,type,inverse=""):
        self.precond = Preconditioner(self.A, type, inverse=inverse)

    def draw(self):
        Draw(mesh)

    def assemble(self):
        self.A.Assemble()
        self.M.Assemble()

