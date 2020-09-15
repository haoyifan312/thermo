import numpy as np


def staverman_guggenheim(r: np.array, q: np.array, x: np.array, pow=1):
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
    return np.log(phiI) + z/2.0*qI*np.log(thetaI/phiI) + lI - phiI*sumxIlI


def get_tauij(T: float, p: np.array):
    """
    tauij = exp-(a + b/T + clnT + dT + e/T^2)
    T - K, float
    """
    if T <= 0.0:
        raise Exception("UNIQUAC tau calculation with infeasible T")
    a = p[0] + p[1]/T + p[2]*np.log(T) + p[3]*T + p[4]/T/T
    return np.exp(a)


class UNIQUAC():
    def __init__(self, size: int, r: np.array, q: np.array, tauij: np.ndarray):
        self.size = size
        self.x = np.empty(size, dtype=float)
        self.r = r
        self.q = q
        self.tauij = tauij
        assert r.size == size
        assert q.size == size
        assert tauij.shape[0] == size
        assert tauij.shape[1] == size

    def __repr__(self):
        return f'{self.__class__.__name__}({self.size})'

    def normalize_x(self):
        """
        normalize x vector
        :return: void
        """
        total = self.x.sum()
        assert total > 0.0
        self.x = self.x / total

    def compute(self):
        """
        compute UNIQUAC activity coefficients
        :return: np.array(size)
        """
        self.normalize_x()
        # compute combinatorial gamma from Stavernman-Guggenheim
        gammaC = staverman_guggenheim(self.r, self.q, self.x)

        # start calculate residual term
        # sumi_xi*qi
        sumqx = (self.x * self.q).sum()
        # sumj_xj*qj*tauji
        sumj_xj_qj_tauji = np.empty(self.size, dtype=float)
        for i in range(self.size):
            sumj_xj_qj_tauji[i] = (self.x * self.q * self.tauij[:, i]).sum()
        # last term in index i
        last_term = np.empty(self.size, dtype=float)
        for i in range(self.size):
            last_term[i] = (self.q * self.x * self.tauij[i, :] / sumj_xj_qj_tauji).sum()
        # add up
        gammaR = self.q * (1.0 - np.log(sumj_xj_qj_tauji/sumqx) - last_term)
        return gammaC + gammaR
        