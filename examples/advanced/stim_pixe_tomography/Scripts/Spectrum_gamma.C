//***********************************************************************************************************
// Spectrum_gamma.C
// Root command file
// Type: root Spectrum_gamma.C
//
// It visualizes the spectrum of X-rays and plots a histogram by reading
// simulation result GammaAtCreation.dat or GammaAtExit.dat
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

// Define a structure to read and write each event in the required binary format
struct RunInfo
{
  // uint_16t
  uint8_t projectionIndex;  // 1 byte
  uint16_t sliceIndex;  //
  uint16_t pixelIndex;
  uint32_t nbParticle;  // 4 bytes int
};

struct ParticleInfo
{
  float energy_keV;
  float mx;
  float my;
  float mz;
};

// struct ParticleInfo
//{
//  float energy_keV;
//  float mx;
//  float my;
//  float mz;
//  float x;
//  float y;
//  float z;
//};

void Plot(vector<double>& energies, int bin, double eMin, double eMax)
{
  auto mycanvas = new TCanvas("canvas", "canvas", 800, 50, 600, 600);
  gPad->SetLeftMargin(0.15);

  // unit is in keV
  auto hist = new TH1D("hist (keV)", "Spectrum of photons", bin, eMin, eMax);

  for (int i = 0; i < energies.size(); ++i) {
    hist->Fill(energies[i]);
  }

  hist->Draw();
  hist->GetXaxis()->SetTitle("Energy (keV)");
  hist->GetYaxis()->SetTitle("Counts");
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();

  mycanvas->Print("spectrum_gamma.png");
}

void Spectrum_gamma()
{
  FILE* input = fopen("../build/GammaAtExit.dat", "rb");

  if (input == NULL) {
    printf("error for opening the input file\n");
    return;
  }

  //***********************************************************************
  //**************************Selection parameters (begin)*****************
  //***********************************************************************

  const int nbProjection = 10;
  const int nbSlice = 1;
  const int nbPixel = 20;

  int projection_index_begin = 0;  // starter of the projection selected
  int projection_index_end = 0;  // end of the projection selected

  int slice_index_begin = 0;  // starter of the slice selected
  int slice_index_end = 0;  // end of the slice selected

  //********************Parameters for spectrum***************************
  int bin = 100;
  double eMin = 0;  // keV
  double eMax = 0;  // keV

  //***********************************************************************
  //**************************Selection parameters (end)*******************
  //***********************************************************************

  RunInfo runInfo;
  vector<double> energies;
  int runID = -1;  // index of simulations, namely runID, starting from 0
  // while(!feof(input)) //if not the end, read
  while (fread(&runInfo, sizeof(RunInfo), 1, input)) {
    runID++;
    int nbParticle = runInfo.nbParticle;

    // ***********the following codes are used
    // if**************************************(begin)
    // ***********the index of projection, slice and pixel is not correctly
    // configured in the simulation
    runInfo.projectionIndex = runID / (nbSlice * nbPixel);
    int remain = runID % (nbSlice * nbPixel);
    runInfo.sliceIndex = remain / nbPixel;
    runInfo.pixelIndex = remain % nbPixel;
    //************************************************************************(end)

    if (!nbParticle) continue;
    std::vector<ParticleInfo> particles(nbParticle);
    fread(&particles[0], sizeof(ParticleInfo), nbParticle, input);

    if (runInfo.projectionIndex >= projection_index_begin
        && runInfo.projectionIndex <= projection_index_end)
    {
      if (runInfo.sliceIndex >= slice_index_begin && runInfo.sliceIndex <= slice_index_end) {
        for (int i = 0; i < nbParticle; ++i) {
          // printf("--%d, %.9e\n", i, particles[i].energy_keV);

          energies.push_back(particles[i].energy_keV);
          if (particles[i].energy_keV > eMax) eMax = particles[i].energy_keV;
        }
      }
    }
    else
      break;
  }

  fclose(input);
  Plot(energies, bin, eMin, eMax + 10);
}
