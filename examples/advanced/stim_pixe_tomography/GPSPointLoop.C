//***********************************************************************************************************
// GPSPointLoop.C
// Root command file
// Type: root GPSPointLoop.C
//
// It generates a macro file to run the simulation.
//
// More information is available in UserGuide
// Created by Z.LI LP2i Bordeaux 2022
//***********************************************************************************************************

// include <stdio.h>
// include <string.h>
// include <stdint.h>
// include <vector>
// include <math.h>
// using namespace std;

void GPSPointLoop()
{
  gSystem->CopyFile("pixe3d_initial.mac", "pixe3d.mac", true);
  FILE* pfile = fopen("pixe3d.mac", "a+");

  //    gSystem->CopyFile("pixe3d_initial.mac", "pixe3d_stim.mac", true);
  //    FILE* pfile = fopen("pixe3d_stim.mac", "a+");

  //***********************************************
  //***(begin)** Define scan parameters************
  //***********************************************
  //***********************************************
  int NumberOfProjections = 10;  // Define the number of Projections from zero to TotalAngleSpan
                                 // (last value "TotalAngleSpan" is excluded)
  int NumberOfSlices = 1;  // Define the number of Slices
  int NumberOfPixels = 20;  // Define the number of Pixels for square YZ Scan
  double TotalAngleSpan = 180;  // scan angular range in degrees
  double ScanSize = 40 * 1.8;  // unit um, scan size for cube of 40 um
  //  double ScanSize   = 42.48*1.8;       // unit um, scan size for C.elegans
  //  double ScanSize   = 500;       // unit um, scan size for GDP
  double ScanHeight = ScanSize;  // Height of the scan, it depends on the need
  // double ScanHeight = 201.127;   //Height of the scan for STIM-T simulatio of C. elegans, for 128
  // slices
  int NbParticles = 1000000;
  double energy = 1.5;  // MeV
  char typeParticle[10] = "proton";

  double PixelWidth = 1. * ScanSize / NumberOfPixels;  // Width of each pixel
  double SliceHeight = 1. * ScanHeight / NumberOfSlices;  // Height of each
                                                          // slice
  double AngleStep =
    1. * TotalAngleSpan / NumberOfProjections;  // angular increment (in degrees)
                                                // between two consecutive projections
  //
  // The beam position is at the center of each pixel
  // Starting position of the beam = StartScan + 0.5 x PixelWidth
  // The scan starts from the bottom left of the square
  //
  double StartScanXY = -0.5 * ScanSize;
  double StartScanZ = -0.5 * ScanHeight;
  // double StartScanZ = 0;

  bool isInterrupted = false;
  int P_interrupt = 0;  // the start of projection index to resume a simulation

  //***********************************************
  //***(end)** Define scan parameters**************
  //***********************************************

  //************************************
  //***(begin)** SCAN IMPLEMENTATION ***
  //************************************
  fprintf(pfile, "/tomography/run/scanParameters %d %d %d\n", NumberOfProjections, NumberOfSlices,
          NumberOfPixels);
  fprintf(pfile, "#\n");

  if (isInterrupted) {
    fprintf(pfile, "/tomography/run/resumeSimulation true\n");
    fprintf(pfile, "/tomography/run/resumeProjectionIndex %d\n", P_interrupt);
    fprintf(pfile, "#\n");
  }
  fprintf(pfile, "/run/initialize\n");
  fprintf(pfile, "#\n");
  fprintf(pfile, "/run/printProgress 500000\n");
  fprintf(pfile, "#\n");
  fprintf(pfile, "# Source definition : energy, type\n");
  fprintf(pfile, "#\n");
  fprintf(pfile, "/gps/energy %.2f MeV\n", energy);
  fprintf(pfile, "/gps/particle %s\n", typeParticle);
  fprintf(pfile, "#\n");
  fprintf(pfile, "# SOURCE POSITION AND DIRECTION\n");
  fprintf(pfile, "#\n");
  for (int projectionIndex = 0; projectionIndex < NumberOfProjections;
       ++projectionIndex)  // projections
  {
    if (isInterrupted) {
      if (projectionIndex < P_interrupt) continue;
    }

    for (int sliceIndex = 0; sliceIndex < NumberOfSlices; ++sliceIndex)  // slices
    {
      // if(sliceIndex<15) continue;
      for (int pixelIndex = 0; pixelIndex < NumberOfPixels; ++pixelIndex)  // pixels
      {
        double px = cos(projectionIndex * AngleStep * TMath::DegToRad());  // beam direction
        double py = sin(projectionIndex * AngleStep * TMath::DegToRad());
        double pz = 0.0;
        double x =
          StartScanXY * px - (StartScanXY + (pixelIndex + 0.5) * PixelWidth) * py;  // beam position
        double y = StartScanXY * py + (StartScanXY + (pixelIndex + 0.5) * PixelWidth) * px;
        double z = StartScanZ + (sliceIndex + 0.5) * SliceHeight;
        // z = 18.07;  //if z is fixed
        //  z = z + 1.953125;
        //  z = z + 3.90625;
        fprintf(pfile, "/gps/direction %f %f %f\n", px, py, pz);
        // fprintf(pfile, "/gps/pos/centre %.6f %.6f %.6f um\n",x, y, z );
        fprintf(pfile, "/gps/pos/centre %f %f %f um\n", x, y, z);
        fprintf(pfile, "/run/beamOn %d\n", NbParticles);
        fprintf(pfile, "#\n");
      }
    }
  }
  fclose(pfile);

  //************************************
  //***(end)** SCAN IMPLEMENTATION ***
  //************************************
}
