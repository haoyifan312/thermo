import numpy as np


def staverman_guggenheim(rq, x, pow=1):
    """
    Combinatorial part of UNIQUAC gamma
    """
    rI = np.power(rq[0,:],pow)
    qI = np.power(rq[1,:],pow)
    # phiI/xI and thetaI/xI to avoid divide by 0
    phiI = rI/np.sum(x*rI)
    thetaI = qI/np.sum(x*qI)
    z = 10
    lI = z/2.0*(rI-qI)-(rI-1.0)
    sumxIlI = np.sum(x*lI)
    lnGamC = np.log(phiI) + z/2.0*qI*np.log(thetaI/phiI) + lI - phiI*sumxIlI
    return lnGamC


def get_tauij(p, T):
    """
    tauij = exp(a + b/T + clnT + dT + e/T^2)
    T - K
    """
    if T <= 0.0:
        raise Exception("UNIQUAC tau calculation with infeasible T")
    a = p[0] + p[1]/T + p[2]*np.log(T) + p[3]*T + p[4]/T/T
    return np.exp(-a/T)