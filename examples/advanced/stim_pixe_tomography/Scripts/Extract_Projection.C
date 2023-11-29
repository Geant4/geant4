//***********************************************************************************************************
// Extract_Slice.C
// Root command file
// Use it by typing in the command line of Root terminal: root Extract_Slice.C
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

// to extract a certain slice or slices

void Extract_Projection()
{
  // FILE *in =fopen("PixeEvent_std_AtCreation.DAT","rb");
  FILE* in = fopen("../build/PixeEvent_std_AtExit_Detector135_Aperture70.DAT", "rb");

  // FILE* out = fopen("PixeEvent_std_AtCreation_50Projections.DAT","wb");
  FILE* out = fopen("../build/PixeEvent_std_AtExit_Detector135_Aperture70_50Projections.DAT", "wb");

  if (in == NULL) {
    printf("error for opening the intput file\n");
    return;
  }

  PixeEvent p;
  PixeEvent pp;
  vector<int> valid_projections;
  for (int i = 0; i < 50; ++i) {
    int p = 2 * i;
    valid_projections.push_back(p);
  }

  while (fread(&p, 7, 1, in)) {
    int key = p.projectionIndex;

    if (std::find(valid_projections.begin(), valid_projections.end(), key)
        != valid_projections.end()) {
      pp.energy_10eV = p.energy_10eV;
      pp.projectionIndex = p.projectionIndex / 2;
      pp.sliceIndex = p.sliceIndex;  // index of slices should be reset, starting from 0
      pp.pixelIndex = p.pixelIndex;
      pp.pixelIndex = p.pixelIndex;
      // printf("__ProjectionIndex=%d, SliceIndex=%d, PixelIndex=%d, Energy_10eV=%d\n",
      // pp.projectionIndex, pp.sliceIndex, pp.pixelIndex, pp.energy_10eV);
      fwrite(&pp, 7, 1, out);
    }
  }
  fclose(in);
  fclose(out);
}
