#ifndef _ITRF_H_
#define _ITRF_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

#ifndef PI
#define PI          3.1415926535897932  /* pi */
#endif

#ifndef D2R
#define D2R         (PI/180.0)          /* deg to rad */
#endif

#ifndef R2D
#define R2D         (180.0/PI)          /* rad to deg */
#endif

#ifndef RE_WGS84
#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#endif

#ifndef FE_WGS84
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */
#endif 

#ifndef RE_GRS80
#define RE_GRS80    6378137.0           /* earth semimajor axis (WGS84) (m) */
#endif

#ifndef FE_GRS80
#define FE_GRS80    (1.0/298.257222101) /* earth flattening (WGS84) */
#endif

#ifndef MAS2R
#define MAS2R       (0.001/3600.0*D2R)
#endif

/* transformation parameters are from https://itrf.ign.fr/en/solutions/transformations */
void get_parameter_itrf2020_to_itrf2014(double* T, double* R, double *D, double t);
void get_parameter_itrf2020_to_itrf2008(double* T, double* R, double *D, double t);
void get_parameter_itrf2020_to_itrf2000(double* T, double* R, double *D, double t);
void get_parameter_itrf2020_to_itrf1996(double* T, double* R, double *D, double t);

/* transformation */
void coordinate_transformation_to(double *xyz_src, double *xyz_to, double *T, double *R, double D);

/* ecef xyz to latitude,longitude and height => ITRF need to use RE_GRS80 and FE_GRS80 */
void ecef2pos_(const double *r, double *pos, double a, double f);


#ifdef __cplusplus
}
#endif


#endif
