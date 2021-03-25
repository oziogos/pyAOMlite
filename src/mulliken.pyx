cimport wrapper

def overlap(double X1, double Y1, double Z1, double X2, double Y2, double Z2,
            double mu1, double mu2, int type1, int type2):
    return wrapper.overlap(X1, Y1, Z1, X2, Y2, Z2, mu1, mu2, type1, type2)
