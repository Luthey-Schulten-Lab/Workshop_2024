# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False
# cython: initializedcheck=False

import numpy as np
cimport numpy as cnp

cdef class RHS:
    cdef int nsp, nrxn
    cdef int[:,:] stoch
    cdef int[:,:] dep
    cdef double[:] dys
    cdef double[:] ks
    
    def __init__(self, stochIn, depIn, ksIn):
        self.nsp = stochIn.shape[1]
        self.nrxn = stochIn.shape[0]
        self.stoch = stochIn
        self.dep = depIn
        self.ks = ksIn
        self.dys = np.zeros(self.nsp, dtype=np.float64)
        
    cdef void rhs(self, double[:] ys):
        cdef int i,j,k
        cdef double acc
        for i in range(self.nsp):
            self.dys[i] = 0
            
        for i in range(self.nrxn):
            acc = self.ks[i]
            for j in range(self.nsp):
                k = self.dep[i,j]
                while k > 0:
                    acc *= ys[j]
                    k -= 1
            for j in range(self.nsp):
                self.dys[j] += self.stoch[i,j]*acc
                
    def evaluate(self, ys, *args):
        self.rhs(ys)
        return np.array(self.dys)
