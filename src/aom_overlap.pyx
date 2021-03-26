cimport cython
import numpy as np


cdef extern from 'mulliken_functions.h':
    double overlap(double X1, double Y1, double Z1, double X2,
                   double Y2, double Z2, double mu1, double mu2,
                   int type1, int type2)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double _AOM_overlap_calculation(
    int istart, int istop, int jstart, int jstop, double[:] x, double[:] y,
        double[:] z, int[:] STO_id_array, int[:] STO_type_array,
        double[:] STO_mu_array, double[:,:] STO_matrix
):
    cdef double S = 0.0
    cdef int i, j, locali, localj

    for i in range(istart,istop):
        locali=(STO_type_array[i]-2)%4
        for j in range(jstart,jstop):
            localj=(STO_type_array[j]-2)%4
            if locali>0 and localj>0:
                if STO_id_array[i]!=STO_id_array[j]:
                    S=S+STO_matrix[STO_id_array[i]-1, locali]*STO_matrix[STO_id_array[j]-1, localj]\
                        *overlap(x[STO_id_array[i]-1],y[STO_id_array[i]-1],z[STO_id_array[i]-1],
                                 x[STO_id_array[j]-1],y[STO_id_array[j]-1],z[STO_id_array[j]-1],
                                 STO_mu_array[i],STO_mu_array[j],
                                 STO_type_array[i],STO_type_array[j])
                else:
                    if STO_type_array[i]==STO_type_array[j]:
                        S=S+STO_matrix[STO_id_array[i]-1, locali]*STO_matrix[STO_id_array[j]-1, localj]
    return S


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double[:,:] _calculate_overlap_S_matrix(
    double[:] x, double[:] y, double[:] z, int STOs, int[:] STO_id_array,
        int[:] STO_type_array, double[:] STO_mu_array
):
    cdef double[:,:] Smatrix = np.identity(STOs)
    cdef int i, j
    for i in range(STOs):
        for j in range(STOs):
            if STO_id_array[i]!=STO_id_array[j] and i <= j:
                Smatrix[i, j] = overlap(x[STO_id_array[i]-1],y[STO_id_array[i]-1],z[STO_id_array[i]-1],x[STO_id_array[j]-1],y[STO_id_array[j]-1],z[STO_id_array[j]-1],STO_mu_array[i],STO_mu_array[j],STO_type_array[i],STO_type_array[j])

    return Smatrix


def AOM_overlap_calculation(
    int istart, int istop, int jstart, int jstop, double[:] x, double[:] y,
        double[:] z, int[:] STO_id_array, int[:] STO_type_array,
        double[:] STO_mu_array, double[:,:] STO_matrix
):
    return _AOM_overlap_calculation(istart, istop, jstart, jstop, x, y,
        z, STO_id_array, STO_type_array, STO_mu_array, STO_matrix)


def calculate_overlap_S_matrix(
    double[:] x, double[:] y, double[:] z, int STOs, int[:] STO_id_array,
        int[:] STO_type_array, double[:] STO_mu_array
):
    s_matrix = np.array(_calculate_overlap_S_matrix(
        x, y, z, STOs, STO_id_array, STO_type_array, STO_mu_array
    ))
    return s_matrix + s_matrix.T - np.diag(np.diag(s_matrix))
