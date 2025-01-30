#include "itrf.h"

extern void get_parameter_itrf2020_to_itrf2014(double* T, double* R, double *D, double t)
{
    /* from ITRF2020 to ITRF2000 both 2015.0 epoch */
    double pT[3] = { -1.4 / 1000.0,-0.9 / 1000.0,  1.4 / 1000.0 };
    double vT[3] = {  0.0 / 1000.0,-0.1 / 1000.0,  0.2 / 1000.0 };
    double pD = -0.42 * 1.0e-9;
    double vD =  0.00 * 1.0e-9;
    double pR[3] = { 0 * MAS2R, 0 * MAS2R, 0 * MAS2R };
    double vR[3] = { 0 * MAS2R, 0 * MAS2R, 0 * MAS2R };

    double dt = t - 2015.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;

    R[0] = pR[0] + vR[0] * dt;
    R[1] = pR[1] + vR[1] * dt;
    R[2] = pR[2] + vR[2] * dt;

    *D = pD + vD * dt;

    return;
}

extern void get_parameter_itrf2020_to_itrf2008(double* T, double* R, double *D, double t)
{
    /* from ITRF2020 to ITRF2000 both 2015.0 epoch */
    double pT[3] = {  0.2 / 1000.0, 1.0 / 1000.0,  3.3 / 1000.0 };
    double vT[3] = {  0.0 / 1000.0,-0.1 / 1000.0,  0.1 / 1000.0 };
    double pD = -0.29 * 1.0e-9;
    double vD =  0.00 * 1.0e-9;
    double pR[3] = { 0 * MAS2R, 0 * MAS2R, 0 * MAS2R };
    double vR[3] = { 0 * MAS2R, 0 * MAS2R, 0 * MAS2R };

    double dt = t - 2015.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;

    R[0] = pR[0] + vR[0] * dt;
    R[1] = pR[1] + vR[1] * dt;
    R[2] = pR[2] + vR[2] * dt;

    *D = pD + vD * dt;

    return;
}

extern void get_parameter_itrf2020_to_itrf2000(double* T, double* R, double *D, double t)
{
    /* from ITRF2020 to ITRF2000 both 2015.0 epoch */
    double pT[3] = { -0.2 / 1000.0, 0.8 / 1000.0, -34.2 / 1000.0 };
    double vT[3] = {  0.1 / 1000.0, 0.0 / 1000.0, - 1.7 / 1000.0 };
    double pD = 2.25 * 1.0e-9;
    double vD = 0.11 * 1.0e-9;
    double pR[3] = { 0 * MAS2R, 0 * MAS2R, 0 * MAS2R };
    double vR[3] = { 0 * MAS2R, 0 * MAS2R, 0 * MAS2R };

    double dt = t - 2015.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;

    R[0] = pR[0] + vR[0] * dt;
    R[1] = pR[1] + vR[1] * dt;
    R[2] = pR[2] + vR[2] * dt;

    *D = pD + vD * dt;

    return;
}

extern void get_parameter_itrf2020_to_itrf1996(double* T, double* R, double *D, double t)
{
    /* from ITRF2020 to ITRF2000 both 2015.0 epoch */
    double pT[3] = {  6.5 / 1000.0,-3.9 / 1000.0, -77.9 / 1000.0 };
    double vT[3] = {  0.1 / 1000.0,-0.6 / 1000.0, - 3.1 / 1000.0 };
    double pD = 3.98 * 1.0e-9;
    double vD = 0.12 * 1.0e-9;
    double pR[3] = { 0 * MAS2R, 0 * MAS2R, 0.36 * MAS2R };
    double vR[3] = { 0 * MAS2R, 0 * MAS2R, 0.02 * MAS2R };

    double dt = t - 2015.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;

    R[0] = pR[0] + vR[0] * dt;
    R[1] = pR[1] + vR[1] * dt;
    R[2] = pR[2] + vR[2] * dt;

    *D = pD + vD * dt;

    return;
}

extern void coordinate_transformation_to(double *xyz_src, double *xyz_to, double *T, double *R, double D)
{
	xyz_to[0] = xyz_src[0] + T[0] + D * xyz_src[0] - R[2] * xyz_src[1] + R[1] * xyz_src[2];
	xyz_to[1] = xyz_src[1] + T[1] + D * xyz_src[1] + R[2] * xyz_src[0] - R[0] * xyz_src[2];
	xyz_to[2] = xyz_src[2] + T[2] + D * xyz_src[2] - R[1] * xyz_src[0] + R[0] * xyz_src[1];
	return;
}

/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
static double dot_(const double *a, const double *b, int n)
{
    double c=0.0;
    
    while (--n>=0) c+=a[n]*b[n];
    return c;
}
/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
*-----------------------------------------------------------------------------*/
extern void ecef2pos_(const double *r, double *pos, double a, double f)
{
    double e2=f*(2.0-f),r2=dot_(r,r,2),z,zk,v=a,sinp;
    
    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=a/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
    pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
    pos[2]=sqrt(r2+z*z)-v;
}