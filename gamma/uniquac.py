import numpy as np


def staverman_guggenheim(r, q, x, pow=1):
    """
    Combinatorial part of UNIQUAC gamma
    r - volume parameter
    q - surface area parameter
    x - molefraction 
    these three numpy 1-D array should have the same size
    """
    # check size
    if r.size != x.size or q.size != x.size:
        raise Exception('SG model input x, r, q are not in the same size')
        
    rI = np.power(r, pow)
    qI = np.power(q, pow)
    # phiI/xI and thetaI/xI to avoid divide by 0
    phiI = rI/np.sum(x*rI)
    thetaI = qI/np.sum(x*qI)
    z = 10
    lI = z/2.0*(rI-qI)-(rI-1.0)
    sumxIlI = np.sum(x*lI)
    lnGamC = np.log(phiI) + z/2.0*qI*np.log(thetaI/phiI) + lI - phiI*sumxIlI
    return lnGamC


def get_tauij(T, *p):
    """
    tauij = exp-(a + b/T + clnT + dT + e/T^2)
    T - K
    """
    if T <= 0.0:
        raise Exception("UNIQUAC tau calculation with infeasible T")
    a = p[0] + p[1]/T + p[2]*np.log(T) + p[3]*T + p[4]/T/T
    return np.exp(-a/T)


class UNIQUAC():
    def __init__(self, size):
        self.size = size

    def __repr__(self):
        return f'{self.__class__.__name__}({self.size})'

    @property
    def x(self):
        "mole fraction array"
        return self._x

    @x.setter
    def x(self, x):
        if x.size != self.size:
            raise Exception('UNIQUAC x size inconsistency')
        self._x = x

    @property
    def r(self):
        "volume parameter array"
        return self._r

    @r.setter
    def r(self, r):
        if r.size != self.size:
            raise Exception('UNIQUAC r size inconsistency')
        self._r = r

    @property
    def q(self):
        "surface area parameter array"
        return self._q

    @q.setter
    def q(self, q):
        if q.size != self.size:
            raise Exception('UNIQUAC size inconsistency')
        self._q = q

    @property
    def tauij(self):
        "UNIQUAC interaction parameter"
        return self._tauij

    @tauij.setter
    def tauij(self, tauij):
        if tauij.shape(0) != self.size or \
                tauij.shape(1) != self.size:
            raise Exception('UNIQUAC size inconsistency')
        self._tauij = tauij
            
    def compute(self):
        if self.size <= 0 or \
                self.x.size != self.size or \
                self.r.size != self.size or \
                self.q.size != self.size:
            raise Exception('UNIQUAC is not correctly setup')

        # compute combinatorial gamma from Stavernman-Guggenheim
        gammaC = staverman_guggenheim(self.r, self.q, self.x)



        