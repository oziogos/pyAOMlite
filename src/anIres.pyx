cdef extern from 'anIres_functions.h':
    double anIres_func(double X1, double Y1, double Z1, double a1, int lx1,
                  int ly1, int lz1, double X2, double Y2, double Z2,
                  double a2, int lx2, int ly2, int lz2)


def anIres(double X1, double Y1, double Z1, double a1, int lx1, int ly1,
           int lz1, double X2, double Y2, double Z2, double a2, int lx2,
           int ly2,int lz2):

    return anIres_func(X1, Y1, Z1, a1, lx1, ly1, lz1, X2, Y2, Z2, a2, lx2, ly2, lz2)
