/* -----------------------------------------------------------------------------

 Unified STO projection / AOM overlap code

 Beta version: 1.0
 1-Mar-2021

 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk

 For more information, examine the README file in the parent directory.

-------------------------------------------------------------------------------- */

int equal(double A,double B);
long int factorial(int n);
void cross_product(double ax,double ay,double az,double bx,double by,double bz,double *crossx,double *crossy,double *crossz);
void MakeLocalCoordsSystems(double X1,double Y1,double Z1,double X2,double Y2,double Z2,double *exb_X,double *exb_Y,double *exb_Z,double *eyb_X,double *eyb_Y,double *eyb_Z,double *ezb_X,double *ezb_Y,double *ezb_Z);
double Afunction(int k,double p);
double Bfunction(int k,double p,double t);
double pvalue(double X1,double Y1,double Z1,double X2,double Y2,double Z2,double mu1,double mu2);
double tvalue(double mu1,double mu2);
double Smulliken_s_s(double p,double t,int type1,int type2);
double Smulliken_s_psigma(double p,double t,int type1,int type2);
double Smulliken_psigma_psigma(double p,double t,int type1,int type2);
double Smulliken_ppi_ppi(double p,double t,int type1,int type2);
double overlap(double X1,double Y1,double Z1,double X2,double Y2,double Z2,double mu1,double mu2,int type1,int type2);
