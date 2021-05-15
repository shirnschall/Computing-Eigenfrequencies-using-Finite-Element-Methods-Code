from ngsolve.la import EigenValues_Preconditioner
from geometry import *


import netgen.gui
import netgen.geom2d as geom2d
from netgen.csg import *
from netgen.geom2d import unit_square
from ngsolve import *
import math
import scipy.linalg
from scipy import random
from netgen.NgOCC import *
from numpy import linalg as la
import time

# finding eigenvalues of   A u = lambda C u
# condition number = largest/smallest
class EVSolver():

    def __init__(self):
        pass

    def getF(self,l,shift=0):
        return (sqrt(abs(l - shift)) / (2 * math.pi))
        
        
    def lobpsd(self, geo , num, maxNumIterations=50, eps= 10e-6):
        lam = EigenValues_Preconditioner(mat=geo.A.mat, pre=geo.precond)

        u = GridFunction(geo.fes, multidim=num)

        #using multivectors for better performance
        uvecs = MultiVector(u.vec, num)
        vecs = MultiVector(u.vec, 2 * num)

        for v in vecs[0:num]:
            v.SetRandom()
        uvecs[:] = geo.precond * vecs[0:num]
        lams = Vector(num * [1])
        numIterations=0
        res=[1]*(maxNumIterations)

        #weird python do while loop
        while True:
            numIterations+=1
            vecs[0:num] = geo.A.mat * uvecs - (geo.M.mat * uvecs).Scale(lams)
            vecs[num:2 * num] = geo.precond * vecs[0:num]

            # T-norm res
            r = InnerProduct(vecs[num], vecs[0])
            for i in range(1, num):
                tmp = InnerProduct(vecs[num + i], vecs[i])
                if (r < tmp):
                    r = tmp
            res[numIterations-1] = r

            vecs[0:num] = uvecs

            vecs.Orthogonalize()

            asmall = InnerProduct(vecs, geo.A.mat * vecs)
            msmall = InnerProduct(vecs, geo.M.mat * vecs)

            ev, evec = scipy.linalg.eigh(a=asmall, b=msmall)
            prev = lams
            lams = Vector(ev[0:num])
            print(numIterations, ":", [self.getF(l, geo.shift) for l in lams])
            print("res:", res[numIterations-1],"\n")
            
            uvecs[:] = vecs * Matrix(evec[:, 0:num])


            if(abs(res[numIterations-1]) < eps or numIterations>=maxNumIterations):
                break

        for j in range(num):
            u.vecs[j][:] = 0.0
            u.vecs[j].data += uvecs[j]

        Draw(u, geo.mesh, "mode")
        SetVisualization(deformation=True)
        return ev,evec,res


    def lobpcg(self, geo , num, maxNumIterations=50, eps= 10e-6):
        lam = EigenValues_Preconditioner(mat=geo.A.mat, pre=geo.precond)

        u = GridFunction(geo.fes, multidim=num)

        # using multivectors for better performance
        uvecs = MultiVector(u.vec, num)
        vecs = MultiVector(u.vec, 2 * num)

        for v in vecs[0:num]:
            v.SetRandom()
        uvecs[:] = geo.precond * vecs[0:num]
        lams = Vector(num * [1])
        numIterations = 0
        res=[1]*(maxNumIterations)

        # weird python do while loop
        while True:
            numIterations += 1
            vecs[0:num] = geo.A.mat * uvecs[0:num] - (geo.M.mat * uvecs[0:num]).Scale(lams)
            vecs[num:2 * num] = geo.precond * vecs[0:num]

            # T-norm res
            r = InnerProduct(vecs[num], vecs[0])
            for i in range(1, num):
                tmp = InnerProduct(vecs[num + i], vecs[i])
                if (r < tmp):
                    r = tmp
            res[numIterations-1] = r

            vecs[0:num] = uvecs[0:num]

            vecs.Orthogonalize()

            asmall = InnerProduct(vecs, geo.A.mat * vecs)
            msmall = InnerProduct(vecs, geo.M.mat * vecs)

            ev, evec = scipy.linalg.eigh(a=asmall, b=msmall)
            prev = lams
            lams = Vector(ev[0:num])
            print(numIterations, ":", [self.getF(l, geo.shift) for l in lams])
            print("res:", res[numIterations-1], "\n")

            if (numIterations==1):
                tmp = MultiVector(u.vec, 2 * num)
                tmp[0:2*num] = vecs
                vecs = MultiVector(u.vec, 3 * num)
                vecs[0:2*num] = tmp[0:2 * num]

            uvecs[0:num] = vecs * Matrix(evec[:, 0:num])

            #todo: use span{w^i,x^i,p^i} instead of span{w^i,x^i,x^{i-1}} for better stability
            vecs[2 * num:3 * num] = vecs[0:num]

            if (abs(res[numIterations-1]) < eps or numIterations >= maxNumIterations):
                break

        for j in range(num):
            u.vecs[j][:] = 0.0
            u.vecs[j].data += uvecs[j]

        Draw(u)
        #SetVisualization(deformation=True)
        return ev, evec, res
