//***********************************************************************************************************
//     Check_PixeEventFile.C
// Root command file
// Use it by typing in the command line of Root terminal: root Check_PixeEventFile.C
//
//
// More information is available in UserGuide
// Created by Z.LI LP2i Bordeaux 2022
//***********************************************************************************************************

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <vector>
// using namespace std;

#define PI 3.14159265f

// Define a structure to read and write each event in the required binary format
struct PixeEvent
{
  uint16_t energy_10eV;
  uint16_t pixelIndex;
  uint16_t sliceIndex;
  uint8_t projectionIndex;
};
struct ParticleInfo
{
  float energy_keV;
  float mx;
  float my;
  float mz;
};
struct RunInfo
{
  // uint_16t
  uint8_t projectionIndex;  // 1 byte
  uint16_t sliceIndex;  //
  uint16_t pixelIndex;
  uint32_t nbParticle;  // 4 bytes int
};

double DegreeToRadian(double degree)
{
  return (PI * degree / 180.);
}

struct Point
{
  double m_x;
  double m_y;
  double m_z;
};
bool IsDetected(Point poi1, Point poi2, double theta)
{
  double a = (poi1.m_x * poi2.m_x + poi1.m_y * poi2.m_y + poi1.m_z * poi2.m_z)
             / sqrt(poi1.m_x * poi1.m_x + poi1.m_y * poi1.m_y + poi1.m_z * poi1.m_z)
             / sqrt(poi2.m_x * poi2.m_x + poi2.m_y * poi2.m_y + poi2.m_z * poi2.m_z);
  if (a > 1.0) a = 1;
  if (a < -1.0) a = -1;
  double r = acos(a);
  if (r > theta)
    return false;
  else
    return true;
}

void Check_PixeEventFile()
{
  FILE* input2 =
    fopen("../build/PixeEvent_std_AtExit_Detector135_Aperture70_50Projections.DAT", "rb");
  PixeEvent ppp;
  int proj = -1;
  while (fread(&ppp, 7, 1, input2)) {
    if (ppp.projectionIndex != proj) {
      printf("__ProjectionIndex=%d\n", ppp.projectionIndex);
      proj = ppp.projectionIndex;
    }
    // printf("__ProjectionIndex=%d, SliceIndex=%d, PixelIndex=%d, Energy_10eV=%d\n",
    // ppp.projectionIndex, ppp.sliceIndex, ppp.pixelIndex, ppp.energy_10eV);
  }
  fclose(input2);
}
