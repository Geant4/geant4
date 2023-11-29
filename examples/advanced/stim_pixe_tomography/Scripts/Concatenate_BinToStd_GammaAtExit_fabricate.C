//***********************************************************************************************************
//     Concatenate_BinToStd_GammaAtExit_fabricate.C
// Root command file
// Use it by typing in the command line of Root terminal: root
// Concatenate_BinToStd_GammaAtExit_fabricate.C
//
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

void Concatenate_BinToStd_GammaAtExit_fabricate()
{
  //***********************************************************************
  //**************************Detection parameters (begin)*****************
  //***********************************************************************

  const int nbProjection = 100;
  const int nbSlice = 1;
  const int nbPixel = 128;
  double totalAngleSpan = 180.;  // in degree

  double angleOfDetector =
    135.;  // angle of detector relative to the incident direction of the primary protons //
  double distanceObjectDetector = 22.;  // 22 mm
  double radiusOfDetector = 5.;  // 5 mm
  // double theta = atan(radiusOfDetector/distanceObjectDetector); //half apex angle of the right
  // circular cone in radian double theta = 14.726*TMath::DegToRad();    // in radian
  double theta = 70 * TMath::DegToRad();  // in radian
  // double theta = 70*TMath::DegToRad();    // in radian
  // double theta = DegreeToRadian(70);

  int P_interrupt = 1;  // Projection of interruption

  //***********************************************************************
  //**************************Detection parameters (end)*******************
  //***********************************************************************

  // assuming there is one interruption
  FILE* input1 = fopen("../RT7_GDP_1Projs_1Slice_128Pixels_2000000_4MeV/GammaAtExit.dat", "rb");
  FILE* out =
    fopen("../RT7_GDP_1Projs_1Slice_128Pixels_2000000_4MeV/PixeEvent_std_AtExit.DAT", "wb");
  // FILE* temp;
  // temp =fopen("temp.DAT","wb");

  if (input1 == NULL) {
    printf("error for opening the input GammaAtExit.dat file\n");
    return;
  }

  RunInfo runInfo;
  PixeEvent pixeEvent;
  Point centerOfDetector;
  Point gammaMomentum;
  long long count1 = 0;
  long long count2 = 0;
  int runID = -1;  // index of simulations, namely runID, starting from 0
  std::vector<PixeEvent> eventVec;

  // ************************************************************(begin)
  // **********************READ FIRST FILE***********************
  // ************************************************************
  while (fread(&runInfo, sizeof(RunInfo), 1, input1)) {
    runID++;
    runInfo.projectionIndex = runID / (nbSlice * nbPixel);
    int remain = runID % (nbSlice * nbPixel);
    runInfo.sliceIndex = remain / nbPixel;
    runInfo.pixelIndex = remain % nbPixel;
    if (runInfo.projectionIndex == P_interrupt) {
      runID--;
      break;
    }

    int nbParticle = runInfo.nbParticle;

    std::vector<ParticleInfo> gammaAtExit(nbParticle);
    fread(&gammaAtExit[0], sizeof(ParticleInfo), nbParticle, input1);

    //***********************************************************************
    //**************************Print information (begin)********************
    //***********************************************************************

    printf("-1--runId %d, ProjectionIndex=%d, SliceIndex=%d, PixelIndex=%d, nbParticle = %d\n",
           runID, runInfo.projectionIndex, runInfo.sliceIndex, runInfo.pixelIndex, nbParticle);

    //***********************************************************************
    //**************************Print information (end)**********************
    //***********************************************************************

    // angleOfDetector+totalAngleSpan/nbProjection*runInfo.projectionIndex means the angle between
    // source direction and detector, which should be constant when source is rotating
    double ra =
      DegreeToRadian(angleOfDetector + totalAngleSpan / nbProjection * runInfo.projectionIndex);
    centerOfDetector.m_x = distanceObjectDetector * cos(ra);
    centerOfDetector.m_y = distanceObjectDetector * sin(ra);
    centerOfDetector.m_z = 0;

    for (int i = 0; i < nbParticle; ++i) {
      // gamma selection: energy should be lower than 4095*10eV = 49.45 keV
      if (gammaAtExit[i].energy_keV >= 40.95 || gammaAtExit[i].energy_keV <= 0.9)
        continue;  // gamma selection

      gammaMomentum.m_x = gammaAtExit[i].mx;
      gammaMomentum.m_y = gammaAtExit[i].my;
      gammaMomentum.m_z = gammaAtExit[i].mz;

      if (!IsDetected(centerOfDetector, gammaMomentum, theta))
        continue;
      else {
        pixeEvent.energy_10eV = floor(100 * gammaAtExit[i].energy_keV + 0.5);
        pixeEvent.projectionIndex = runInfo.projectionIndex;
        pixeEvent.sliceIndex = runInfo.sliceIndex;
        pixeEvent.pixelIndex = runInfo.pixelIndex;

        eventVec.push_back(pixeEvent);
        count1++;
      }
    }
  }
  printf("---------------Number of PixeEvent in the first file: %lld------------------------\n",
         count1);
  fclose(input1);
  // fclose(temp);

  // ************************************************************(end)
  // **********************READ FIRST FILE***********************
  // ************************************************************

  // ************************************************************(begin)
  // **********************READ SECOND FILE**********************
  // ************************************************************
  // temp =fopen("temp.DAT","rb");

  PixeEvent pp;
  PixeEvent p;

  for (int i = 0; i < nbProjection; ++i) {
    int size = eventVec.size();
    for (int j = 0; j < size; ++j) {
      p = eventVec[j];
      pp.energy_10eV = p.energy_10eV;
      pp.projectionIndex = p.projectionIndex + i;
      pp.sliceIndex = p.sliceIndex;  // index of slices should be reset, starting from 0
      pp.pixelIndex = p.pixelIndex;
      // printf("__ProjectionIndex=%d, SliceIndex=%d, PixelIndex=%d, Energy_10eV=%d\n",
      // pp.projectionIndex, pp.sliceIndex, pp.pixelIndex, pp.energy_10eV);
      fwrite(&pp, 7, 1, out);
    }
  }

  // ************************************************************(end)
  // **********************READ SECOND FILE**********************
  // ************************************************************

  // fclose(temp);
  fclose(out);

  // Recheck the output file in case
  FILE* input2 =
    fopen("../RT7_GDP_1Projs_1Slice_128Pixels_2000000_4MeV/PixeEvent_std_AtExit.DAT", "rb");
  PixeEvent ppp;
  int proj = -1;
  while (fread(&ppp, 7, 1, input2)) {
    if (ppp.projectionIndex != proj) {
      printf("__ProjectionIndex=%d\n", ppp.projectionIndex);
      proj = ppp.projectionIndex;
    }
    // if(proj<20) printf("__ProjectionIndex=%d, SliceIndex=%d, PixelIndex=%d, Energy_10eV=%d\n",
    // ppp.projectionIndex, ppp.sliceIndex, ppp.pixelIndex, ppp.energy_10eV);
  }
  fclose(input2);
}
