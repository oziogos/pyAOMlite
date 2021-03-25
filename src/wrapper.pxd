cdef extern from 'anIres_functions.h':
    double anIres(double X1, double Y1, double Z1, double a1, int lx1,
                  int ly1, int lz1, double X2, double Y2, double Z2,
                  double a2, int lx2, int ly2, int lz2)

cdef extern from 'mulliken_functions.h':
    double overlap(double X1, double Y1, double Z1, double X2,
                   double Y2, double Z2, double mu1, double mu2,
                   int type1, int type2)
