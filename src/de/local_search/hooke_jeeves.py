import numpy as np


def hooke(protein, rho=0.5, eps=1e-04, itermax=100):
    nvars = protein.total_number_of_angles
    startpt = protein.angles
    f = protein

    return hooke_(nvars, startpt, rho, eps, itermax, f)


def hooke_(nvars, startpt, rho, eps, itermax, f):
    newx = startpt.copy()
    xbefore = startpt.copy()

    delta = np.zeros(nvars)

    for i in range(0, nvars):
        if (startpt[i] == 0.0):
            delta[i] = rho
        else:
            delta[i] = rho * abs(startpt[i])

    funevals = 0
    steplength = rho
    iters = 0
    fbefore = f(newx, nvars)
    funevals = funevals + 1
    newf = fbefore

    while (iters < itermax and eps < steplength):
        iters = iters + 1

        for i in range(0, nvars):
            newx[i] = xbefore[i]

        newf, newx, funevals = best_nearby(delta, newx, fbefore, nvars, f, funevals)
        keep = True

        while (newf < fbefore and keep):

            for i in range(0, nvars):
                if (newx[i] <= xbefore[i]):
                    delta[i] = - abs(delta[i])
                else:
                    delta[i] = abs(delta[i])
                tmp = xbefore[i]
                xbefore[i] = newx[i]
                newx[i] = newx[i] + newx[i] - tmp

            fbefore = newf
            newf, newx, funevals = best_nearby(delta, newx, fbefore, nvars, f, funevals)
            if (fbefore <= newf):
                break
            keep = False

            for i in range(0, nvars):
                if (0.5 * abs(delta[i]) < abs(newx[i] - xbefore[i])):
                    keep = True
                    break

        if (eps <= steplength and fbefore <= newf):
            steplength = steplength * rho
            for i in range(0, nvars):
                delta[i] = delta[i] * rho

    endpt = xbefore.copy()

    return iters, endpt


def best_nearby(delta, point, prevbest, nvars, f, funevals):
    z = point.copy()

    minf = prevbest

    for i in range(0, nvars):
        z[i] = point[i] + delta[i]

        ftmp = f(z, nvars)
        funevals = funevals + 1

        if ftmp < minf:
            minf = ftmp
        else:
            delta[i] = - delta[i]
            z[i] = point[i] + delta[i]
            ftmp = f(z, nvars)
            funevals = funevals + 1

            if (ftmp < minf):
                minf = ftmp
            else:
                z[i] = point[i]

    point = z.copy()
    newbest = minf

    return newbest, point, funevals
