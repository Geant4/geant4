//***********************************************************************************************************
// LocateInterruption_ProtonAtExit.C
// Root command file
// Type: root LocateInterruption_ProtonAtExit.C
//
// It is used by reading ProtonAtExit.dat file to locate at which projection the interruption is
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

struct ParticleInfo
{
  float energy_keV;
  float mx;
  float my;
  float mz;
};

// struct ParticleInfo
// {
// float energy_keV;
// float mx;
// float my;
// float mz;
// float x;
// float y;
// float z;
// };

struct RunInfo
{
  // uint_16t
  uint8_t projectionIndex;  // 1 byte
  uint16_t sliceIndex;  //
  uint16_t pixelIndex;
  uint32_t nbParticle;  // 4 bytes int
};

void LocateInterruption_ProtonAtExit()
{
  FILE* input = fopen("../build/ProtonAtExit_1.dat", "rb");
  if (input == NULL) {
    printf("error for opening the input file\n");
    return;
  }

  RunInfo runInfo;
  int projection = 0;  // the projection when interruption occurs

  //***********************************************************************
  //**************************Detection parameters (begin)*****************
  //***********************************************************************

  const int nbProjection = 10;
  const int nbSlice = 128;
  const int nbPixel = 20;

  //***********************************************************************
  //**************************Detection parameters (end)*******************
  //***********************************************************************

  int runID = -1;
  while (fread(&runInfo, sizeof(RunInfo), 1, input)) {
    runID++;

    runInfo.projectionIndex = runID / (nbSlice * nbPixel);
    int remain = runID % (nbSlice * nbPixel);
    runInfo.sliceIndex = remain / nbPixel;
    runInfo.pixelIndex = remain % nbPixel;

    int nbParticle = runInfo.nbParticle;
    std::vector<ParticleInfo> protonAtExit(nbParticle);
    fread(&protonAtExit[0], sizeof(ParticleInfo), nbParticle, input);

    printf("---------ProjectionIndex=%d, SliceIndex=%d, PixelIndex=%d, nbParticle = %d\n",
           runInfo.projectionIndex, runInfo.sliceIndex, runInfo.pixelIndex, nbParticle);

    projection = runInfo.projectionIndex;
  }

  printf("-----------------------It is interrupted at ProjectionIndex = %d--------------------\n",
         projection);
  fclose(input);
}