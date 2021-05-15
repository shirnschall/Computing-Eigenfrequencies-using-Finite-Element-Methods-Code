def lobpcg(mesh,fes,num,A,M,pre,maxNumIterations):
    u = fes.TrialFunction()
    v = fes.TestFunction()

    lam = EigenValues_Preconditioner(mat=A.mat, pre=pre)

    u = GridFunction(fes, multidim=num)

    # using multivectors for better performance
    uvecs = MultiVector(u.vec, num)
    vecs = MultiVector(u.vec, 2 * num)

    for v in vecs[0:num]:
        v.SetRandom()
    uvecs[:] = pre * vecs[0:num]
    lams = Vector(num * [1])
    res = Vector((maxNumIterations) * [1])
    numIterations = 0

    # weird python do while loop
    while True:
        numIterations += 1
        vecs[0:num] = A.mat * uvecs[0:num] - (M.mat * uvecs[0:num]).Scale(lams)
        vecs[num:2 * num] = pre * vecs[0:num]

        #T-norm res
        r = InnerProduct(vecs[num], vecs[0])
        for i in range(1, num):
            tmp = InnerProduct(vecs[num + i], vecs[i])
            if (r < tmp):
                r = tmp
        res[numIterations-1] = r

        vecs[0:num] = uvecs[0:num]

        vecs.Orthogonalize()

        asmall = InnerProduct(vecs, A.mat * vecs)
        msmall = InnerProduct(vecs, M.mat * vecs)

        ev, evec = scipy.linalg.eigh(a=asmall, b=msmall)
        prev = lams
        lams = Vector(ev[0:num])
        print(numIterations, ":", [l for l in lams])

        if (numIterations==1):
            tmp = MultiVector(u.vec, 2 * num)
            tmp[0:2*num] = vecs
            vecs = MultiVector(u.vec, 3 * num)
            vecs[0:2*num] = tmp[0:2 * num]

        uvecs[0:num] = vecs * Matrix(evec[:, 0:num])

        #todo: use span{w^i,x^i,p^i} instead of span{w^i,x^i,x^{i-1}} for better stability
        vecs[2 * num:3 * num] = vecs[0:num]

        print("res:", res[numIterations], "\n")

        if (res[numIterations] < 10e-16 or numIterations >= maxNumIterations):
            break

    for j in range(num):
        u.vecs[j][:] = 0.0
        u.vecs[j].data += uvecs[j]

    Draw(u, mesh, "mode")
    SetVisualization(deformation=True)
    return res,numIterations