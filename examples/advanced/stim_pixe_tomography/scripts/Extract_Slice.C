//***********************************************************************************************************
//     Extract_Slice.C
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

void Extract_Slice()
{
  int start_slice = 0;  // start_slice: the first slice you would like to select
  int end_slice = 0;  // end_slice: the last slice you would like to select

  FILE* in = fopen("../build/PixeEvent_std_AtCreation.DAT", "rb");
  // FILE *in =fopen("PixeEvent_std_AtExit.DAT.DAT","rb");

  FILE* out = fopen("../build/PixeEvent_std_AtCreation_slice.DAT", "wb");
  // FILE* out = fopen("PixeEvent_std_AtExit_slice.DAT","wb");

  if (in == NULL) {
    printf("error for opening the intput file\n");
    return;
  }

  PixeEvent p;
  PixeEvent pp;

  while (fread(&p, 7, 1, in)) {
    if (p.sliceIndex >= start_slice && p.sliceIndex <= end_slice) {
      pp.energy_10eV = p.energy_10eV;
      pp.projectionIndex = p.projectionIndex;
      pp.sliceIndex =
        p.sliceIndex - start_slice;  // index of slices should be reset, starting from 0
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
