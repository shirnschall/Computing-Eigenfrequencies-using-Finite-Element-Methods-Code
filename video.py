for i in range(20):
    vecs[0:num] = a.mat * uvecs - (m.mat * uvecs).Scale(lams)
    vecs[num:2*num] = pre * vecs[0:num]
    vecs[0:num] = uvecs

    vecs.Orthogonalize()

    asmall = InnerProduct (vecs, a.mat * vecs)
    msmall = InnerProduct (vecs, m.mat * vecs)

    ev,evec = scipy.linalg.eigh(a=asmall, b=msmall)
    lams = Vector(ev[0:num])
    print (i, ":", [l/math.pi**2 for l in lams])

    uvecs[:] = vecs * Matrix(evec[:,0:num])

