
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>

#include "itrf.h"


struct nrcan_ppp_t
{
  std::string fname;
  std::string stime;
  std::string etime;
  std::string coord_name;
  double amb_fix_rate;
  double epoch;
  double lat;    /* deg */
  double lon;    /* deg */
  double ht;     /* m */
  double xyz[3]; /* ecef xyz m */
  double sigma95_xyz[3];
  double blh[3];
  double sigma95_NEU[3]; /* accuracy (95%) in NEU [m] */
  nrcan_ppp_t()
  {
    amb_fix_rate = 0;
    epoch = 0;
    lat = 0;
    lon = 0;
    ht = 0;
    xyz[0] = xyz[1] = xyz[2] = 0;
    sigma95_xyz[0] = sigma95_xyz[1] = sigma95_xyz[2] = 0;
    blh[0] = blh[1] = blh[2] = 0;
    sigma95_NEU[0] = sigma95_NEU[1] = sigma95_NEU[2] = 0;
  }
};

static void extract_coordinate(const char *fname)
{
  FILE *fLOG = fopen(fname, "rt");

  if (fLOG == NULL)
    return;

  char buffer[255] = {0};

  double xyz[3] = {0};
  double rms_xyz[3] = {0};
  int idx[3] = {0};
  double mp12 = 0;
  double mp21 = 0;
  double mp15 = 0;
  double mp51 = 0;
  double mp17 = 0;
  double mp71 = 0;
  int iod_slip = 0;
  int mp_slip = 0;
  int oslps = 0;
  int numobs = 0;
  int numepoch = 0;
  int expepoch = 0;
  char oslip_str[50] = {0};
  char name[255] = {0};
  char datestr[255] = {0};
  strcpy(buffer, fname);
  char *temp = strrchr(buffer, '-');
  if (temp)
    strcpy(name, temp + 1);
  temp = strrchr(name, '.');
  if (temp)
    temp[0] = '\0';
  strcpy(datestr, fname);
  temp = strrchr(datestr, '-');
  if (temp)
    temp[0] = '\0';

  nrcan_ppp_t ppp;
  while (fLOG && !feof(fLOG))
  {
    fgets(buffer, sizeof(buffer), fLOG);
    temp = strrchr(buffer, '\n');
    if (temp)
      temp[0] = '\0';
    temp = strrchr(buffer, '\r');
    if (temp)
      temp[0] = '\0';
    if (strlen(buffer) >= 4 && buffer[0] == 'R' && buffer[1] == 'N' && buffer[2] == 'X')
    {
      temp = strrchr(buffer, '.');
      if (temp)
        temp[0] = '\0';
      temp = strrchr(buffer, '-');
      if (temp)
        ppp.fname = std::string(temp + 1);
      else
        ppp.fname = std::string(buffer + 4);
    }
    if (strlen(buffer) >= 4 && buffer[0] == 'B' && buffer[1] == 'E' && buffer[2] == 'G')
    {
      ppp.stime = std::string(buffer + 4);
    }
    if (strlen(buffer) >= 4 && buffer[0] == 'E' && buffer[1] == 'N' && buffer[2] == 'D')
    {
      ppp.etime = std::string(buffer + 4);
    }
    if (strlen(buffer) >= 4 && buffer[0] == 'I' && buffer[1] == 'A' && buffer[2] == 'R')
    {
      ppp.amb_fix_rate = atof(buffer + 4);
    }
    if (strlen(buffer) >= 83 && buffer[0] == 'P' && buffer[1] == 'O' && buffer[2] == 'S' && buffer[6] == 'X')
    {
      ppp.coord_name = std::string(buffer + 8).substr(0, 5);
      int yy = atoi(buffer + 14);
      int doy = atoi(buffer + 17);
      ppp.epoch = 2000.0 + yy + doy / 365.0;
      ppp.xyz[0] = atof(buffer + 44);
      ppp.sigma95_xyz[0] = atof(buffer + 73);
    }
    if (strlen(buffer) >= 83 && buffer[0] == 'P' && buffer[1] == 'O' && buffer[2] == 'S' && buffer[6] == 'Y')
    {
      ppp.coord_name = std::string(buffer + 8).substr(0, 5);
      int yy = atoi(buffer + 14);
      int doy = atoi(buffer + 17);
      ppp.epoch = 2000.0 + yy + doy / 365.0;
      ppp.xyz[1] = atof(buffer + 44);
      ppp.sigma95_xyz[1] = atof(buffer + 73);
    }
    if (strlen(buffer) >= 83 && buffer[0] == 'P' && buffer[1] == 'O' && buffer[2] == 'S' && buffer[6] == 'Z')
    {
      ppp.coord_name = std::string(buffer + 8).substr(0, 5);
      int yy = atoi(buffer + 14);
      int doy = atoi(buffer + 17);
      ppp.epoch = 2000.0 + yy + doy / 365.0;
      ppp.xyz[2] = atof(buffer + 44);
      ppp.sigma95_xyz[2] = atof(buffer + 73);
    }
    if (strlen(buffer) >= 83 && buffer[0] == 'P' && buffer[1] == 'O' && buffer[2] == 'S' && buffer[4] == 'L' && buffer[5] == 'A' && buffer[6] == 'T')
    {
      ppp.coord_name = std::string(buffer + 8).substr(0, 5);
      int yy = atoi(buffer + 14);
      int doy = atoi(buffer + 17);
      ppp.epoch = 2000.0 + yy + doy / 365.0;
      char dd_str[7] = { 0 };
      strncpy(dd_str, buffer + 44, 6);
      dd_str[6] = '\0';
      temp = strchr(dd_str, '-');
      int is_neg = 0;
      if (temp)
      {
          is_neg = 1;
          temp[0] = ' ';
      }
      int dd = atoi(dd_str);// atoi(buffer + 44);
      int mm = atoi(buffer + 50);
      double ss = atof(buffer + 53);
      if (is_neg)
      {
        ppp.blh[0] = -(dd + mm / 60.0 + ss / 3600.0);
      }
      else
      {
        ppp.blh[0] = (dd + mm / 60.0 + ss / 3600.0);
      }
      ppp.sigma95_NEU[0] = atof(buffer + 73);
    }
    if (strlen(buffer) >= 83 && buffer[0] == 'P' && buffer[1] == 'O' && buffer[2] == 'S' && buffer[4] == 'L' && buffer[5] == 'O' && buffer[6] == 'N')
    {
      ppp.coord_name = std::string(buffer + 8).substr(0, 5);
      int yy = atoi(buffer + 14);
      int doy = atoi(buffer + 17);
      ppp.epoch = 2000.0 + yy + doy / 365.0;
      char dd_str[7] = { 0 };
      strncpy(dd_str, buffer + 44, 6);
      dd_str[6] = '\0';
      temp = strchr(dd_str, '-');
      int is_neg = 0;
      if (temp)
      {
          is_neg = 1;
          temp[0] = ' ';
      }
      int dd = atoi(dd_str);// atoi(buffer + 44);
      int mm = atoi(buffer + 50);
      double ss = atof(buffer + 53);
      if (is_neg)
      {
        ppp.blh[1] = -(dd + mm / 60.0 + ss / 3600.0);
      }
      else
      {
        ppp.blh[1] = (dd + mm / 60.0 + ss / 3600.0);
      }
      ppp.sigma95_NEU[1] = atof(buffer + 73);
    }
    if (strlen(buffer) >= 83 && buffer[0] == 'P' && buffer[1] == 'O' && buffer[2] == 'S' && buffer[4] == 'H' && buffer[5] == 'G' && buffer[6] == 'T')
    {
      ppp.coord_name = std::string(buffer + 8).substr(0, 5);
      int yy = atoi(buffer + 14);
      int doy = atoi(buffer + 17);
      ppp.epoch = 2000.0 + yy + doy / 365.0;
      ppp.blh[2] = atof(buffer + 44);
      ppp.sigma95_NEU[2] = atof(buffer + 73);
    }
  }
  if (fabs(ppp.xyz[0]) < 0.001 || fabs(ppp.xyz[1]) < 0.001 || fabs(ppp.xyz[2]) < 0.001)
  {

  }
  else
  {
    /* output solution with valid xyz */
    double blh[3] = { 0 };
    ecef2pos_(ppp.xyz, blh, RE_GRS80, FE_GRS80);
    FILE* fITRF_BLH = fopen("itrf20_blh.txt", "w");
    if (fITRF_BLH)
    {
        fprintf(fITRF_BLH, "%.9f,%.9f,%.4f,%s,%.2f,%.4f,%.4f,%.4f,%s(%6.2f),%s,%s\n", ppp.blh[0], -ppp.blh[1], ppp.blh[2], ppp.fname.c_str(), ppp.amb_fix_rate, ppp.sigma95_NEU[0], ppp.sigma95_NEU[1], ppp.sigma95_NEU[2], ppp.coord_name.c_str(), ppp.epoch, ppp.stime.c_str(), ppp.etime.c_str());
        fclose(fITRF_BLH);
    }
    FILE* fITRF_XYZ = fopen("itrf20_xyz.txt", "w");
    if (fITRF_XYZ)
    {
        fprintf(fITRF_XYZ, "%.4f,%.4f,%.4f,%s,%.2f,%.4f,%.4f,%.4f,%s(%6.2f),%s,%s\n", ppp.xyz[0], ppp.xyz[1], ppp.xyz[2], ppp.fname.c_str(), ppp.amb_fix_rate, ppp.sigma95_xyz[0], ppp.sigma95_xyz[1], ppp.sigma95_xyz[2], ppp.coord_name.c_str(), ppp.epoch, ppp.stime.c_str(), ppp.etime.c_str());
        fclose(fITRF_XYZ);
    }
    FILE* fITRF_2015 = fopen("cfg_itrf2020_2015.txt", "w");
    if (fITRF_2015)
    {
        fprintf(fITRF_2015, "3\n2\n%.2f\n24\n5\n%s\n2\n%.2f\n%s\n0\n0\n", 2015.0, "itrf20_blh.txt", ppp.epoch, "itrf20_2015.txt");
        fprintf(fITRF_2015, "convert ITRF20 coordinate from solution epoch to reference epoch (2015.0), also get the ITRF20 velocity");
        fclose(fITRF_2015);
    }
    FILE* fITRF_2024 = fopen("cfg_itrf2020_2024.txt", "w");
    if (fITRF_2024)
    {
        fprintf(fITRF_2024, "3\n2\n%.2f\n24\n5\n%s\n2\n%.2f\n%s\n0\n0\n", 2024.0, "itrf20_blh.txt", ppp.epoch, "itrf20_2024.txt");
        fprintf(fITRF_2024, "convert ITRF20 coordinate from solution epoch to reference epoch (2024.0), also get the ITRF20 velocity, used to convert ETRS89(ETRF2020)(2024.0) ");
        fclose(fITRF_2024);
    }
    FILE* fITRF_WGS84 = fopen("cfg_itrf2020_wgs84.txt", "w");
    if (fITRF_WGS84)
    {
        fprintf(fITRF_WGS84, "4\n%s\n24\n10\n2\n%.2f\n2\n%0.2f\n4\n%s\n0\n0\n", "wgs84_2139.txt", ppp.epoch, floor(ppp.epoch) + 0.5, "itrf20_xyz.txt");
        fprintf(fITRF_WGS84, "convert ITRF2020 at solution epoch to WGS84(2139) at the mid-year epoch %.2f", floor(ppp.epoch) + 0.5);
        fclose(fITRF_WGS84);
    }
    FILE* fITRF_NAD83 = fopen("cfg_itrf2020_nad83.txt", "w");
    if (fITRF_NAD83)
    {
        fprintf(fITRF_NAD83, "4\n%s\n24\n1\n2\n%.2f\n2\n%0.2f\n4\n%s\n0\n0\n", "nad83_2010.txt", ppp.epoch, 2010.0, "itrf20_xyz.txt");
        fprintf(fITRF_NAD83, "convert ITRF2020 at solution epoch to NAD83(2011) 2010.0");
        fclose(fITRF_NAD83);
    }
    /* GDA2020=ITRF2014(2020.0) */
    FILE* fITRF_GDA2020 = fopen("cfg_itrf2020_gda2020.txt", "w");
    if (fITRF_GDA2020)
    {
        fprintf(fITRF_GDA2020, "4\n%s\n24\n23\n2\n%.2f\n2\n%0.2f\n4\n%s\n0\n0\n", "gda2020.txt", ppp.epoch, 2020.0, "itrf20_xyz.txt");
        fprintf(fITRF_GDA2020, "convert ITRF2020 at solution epoch to GDA2020 reference epoch 2020.0 (same as ITRF2014)");
        fclose(fITRF_GDA2020);
    }
    /* NZGD2000=ITRF96(2000.0) */
    FILE* fITRF_NZGD2000 = fopen("cfg_itrf2020_nzgd2000.txt", "w");
    if (fITRF_NZGD2000)
    {
        fprintf(fITRF_NZGD2000, "4\n%s\n24\n18\n2\n%.2f\n2\n%0.2f\n4\n%s\n0\n0\n", "nzgd2000.txt", ppp.epoch, 2000.0, "itrf20_xyz.txt");
        fprintf(fITRF_NZGD2000, "convert ITRF2020 at solution epoch to NZGD2000 reference epoch 2000.0 (same as ITRF96)");
        fclose(fITRF_NZGD2000);
    }
    /* TUREF ITRF96(2005.0) */
    FILE *fITRF_TUREF = fopen("cfg_itrf2020_turef.txt", "w");
    if (fITRF_TUREF)
    {
        fprintf(fITRF_TUREF, "4\n%s\n24\n18\n2\n%.2f\n2\n%0.2f\n4\n%s\n0\n0\n", "turef.txt", ppp.epoch, 2005.0, "itrf20_xyz.txt");
        fprintf(fITRF_TUREF, "convert ITRF2020 at solution epoch to TUREF reference epoch 2005.0 (same as ITRF96)");
        fclose(fITRF_TUREF);
    }    
    /* INDIAN ITRF2008(2005.0) */
    FILE *fITRF_INDIAN = fopen("cfg_itrf2020_indian.txt", "w");
    if (fITRF_INDIAN)
    {
        fprintf(fITRF_INDIAN, "4\n%s\n24\n22\n2\n%.2f\n2\n%0.2f\n4\n%s\n0\n0\n", "indian.txt", ppp.epoch, 2005.0, "itrf20_xyz.txt");
        fprintf(fITRF_INDIAN, "convert ITRF2020 at solution epoch to INDIAN reference epoch 2005.0 (same as ITRF2008)");
        fclose(fITRF_INDIAN);
    }    
    /* WGS84(G873)=ITRF94(1996.0) => need another step to SLD99 */
    FILE *fITRF_WGS84_G873 = fopen("cfg_itrf2020_wgs84_g873.txt", "w");
    if (fITRF_WGS84_G873)
    {
        fprintf(fITRF_WGS84_G873, "4\n%s-%s\n24\n17\n2\n%.2f\n2\n%0.2f\n4\n%s\n0\n0\n", ppp.fname.c_str(), "wgs84_g873.txt", ppp.epoch, 1996.0, "itrf20_xyz.txt");
        fprintf(fITRF_WGS84_G873, "convert ITRF2020 at solution epoch to WGS84(G873) reference epoch 1996.0 (same as ITRF94)");
        fclose(fITRF_WGS84_G873);
    }    
  }
  if (fLOG) fclose(fLOG);
  FILE* fOUT = fopen("coord-sum.csv", "a");
  if (fOUT)
  {
      fprintf(fOUT, "%14.4f,%14.4f,%14.4f,%14.9f,%14.9f,%14.4f,%15s,%7.2f,%10.4f,%10.4f,%10.4f,%8s(%6.2f),%s,%s,%s\n", ppp.xyz[0], ppp.xyz[1], ppp.xyz[2], ppp.blh[0], ppp.blh[1], ppp.blh[2], ppp.fname.c_str(), ppp.amb_fix_rate, ppp.sigma95_xyz[0], ppp.sigma95_xyz[1], ppp.sigma95_xyz[2], ppp.coord_name.c_str(), ppp.epoch, ppp.stime.c_str(), ppp.etime.c_str(), fname);

      if (1)
      {
          /* from ITRF2020 to ITRF2014 both 2015.0 epoch */
          double T[3] = { 0 };
          double R[3] = { 0 };
          double D = 0;
          double xyz[3] = { 0 };
          double blh[3] = { 0 };
          get_parameter_itrf2020_to_itrf2014(T, R, &D, ppp.epoch);
          coordinate_transformation_to(ppp.xyz, xyz, T, R, D);
          ecef2pos_(xyz, blh, RE_GRS80, FE_GRS80);
          fprintf(fOUT, "%14.4f,%14.4f,%14.4f,%14.9f,%14.9f,%14.4f,%15s,%7.2f,%10.4f,%10.4f,%10.4f,%8s(%6.2f),%s,%s,%s\n", xyz[0], xyz[1], xyz[2], blh[0] * R2D, blh[1] * R2D, blh[2], ppp.fname.c_str(), ppp.amb_fix_rate, ppp.sigma95_xyz[0], ppp.sigma95_xyz[1], ppp.sigma95_xyz[2], "ITRF2014", ppp.epoch, ppp.stime.c_str(), ppp.etime.c_str(), fname);
      }

      if (1)
      {
          /* from ITRF2020 to ITRF2000 both 2015.0 epoch */
          double T[3] = { 0 };
          double R[3] = { 0 };
          double D = 0;
          double xyz[3] = { 0 };
          double blh[3] = { 0 };
          get_parameter_itrf2020_to_itrf2000(T, R, &D, ppp.epoch);
          coordinate_transformation_to(ppp.xyz, xyz, T, R, D);
          ecef2pos_(xyz, blh, RE_GRS80, FE_GRS80);
          fprintf(fOUT, "%14.4f,%14.4f,%14.4f,%14.9f,%14.9f,%14.4f,%15s,%7.2f,%10.4f,%10.4f,%10.4f,%8s(%6.2f),%s,%s,%s\n", xyz[0], xyz[1], xyz[2], blh[0] * R2D, blh[1] * R2D, blh[2], ppp.fname.c_str(), ppp.amb_fix_rate, ppp.sigma95_xyz[0], ppp.sigma95_xyz[1], ppp.sigma95_xyz[2], "ITRF2000", ppp.epoch, ppp.stime.c_str(), ppp.etime.c_str(), fname);
      }

      fclose(fOUT);
  }
  return;
}

int main(int argc, const char *argv[])
{
  if (argc > 1)
    extract_coordinate(argv[1]);
  return 0;
}
