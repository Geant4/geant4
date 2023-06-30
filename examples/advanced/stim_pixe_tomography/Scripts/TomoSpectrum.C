//***********************************************************************************************************
// TomoSpectrum.C
// Root command file
// Type: root TomoSpectrum.C
//
// It visualizes the spectrum of X-rays and plots a graph by reading PixeEvent data.
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

struct PixeEvent
{
  uint16_t energy_10eV;
  uint16_t pixelIndex;
  uint16_t sliceIndex;
  uint8_t projectionIndex;
};

void Plot(int nbChannels, vector<int>& X, vector<int>& Y)
{
  gROOT->Reset();

  auto mycanvas = new TCanvas("canvas", "canvas", 800, 50, 600, 600);
  mycanvas->ToggleEventStatus();

  gPad->SetLeftMargin(0.15);

  auto graph = new TGraph(nbChannels, X.data(), Y.data());
  graph->SetLineColor(8);

  graph->Draw("AL");

  graph->SetLineColor(8);
  graph->SetTitle("TOMO Energy Spectrum");
  graph->GetXaxis()->SetTitle("ADC channels");
  graph->GetYaxis()->SetTitle("Nb events");
  graph->GetXaxis()->CenterTitle();
  graph->GetYaxis()->CenterTitle();

  mycanvas->Print("TomoSpectrum.png");
}

void TomoSpectrum()
{
  FILE* input = fopen("../build/PixeEvent_std_AtExit_Detector135_Aperture70.DAT", "rb");

  if (input == NULL) {
    printf("----------error for opening the input file--------------\n");
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

  //***********************************************************************
  //**************************Selection parameters (end)*******************
  //***********************************************************************

  int nbChannels = 4096;
  vector<int> X(nbChannels);  // save channels 1-4096, index X: 0-4095
  vector<int> Y(nbChannels);  // save event counts for channel 1-4096, index Y: 0-4095
  PixeEvent p;

  while (fread(&p, 7, 1, input)) {
    if (p.projectionIndex >= projection_index_begin && p.projectionIndex <= projection_index_end) {
      if (p.sliceIndex >= slice_index_begin && p.sliceIndex <= slice_index_end) {
        // printf("%d %d %d\n",p.projectionIndex, p.sliceIndex, p.energy_10eV);
        Y[p.energy_10eV - 1] = Y[p.energy_10eV - 1] + 1;
      }
    }
  }
  fclose(input);
  for (int i = 0; i < nbChannels; ++i) {
    X[i] = 1 + i;
  }

  FILE* out = fopen("Spectrum.txt", "wb");
  for (int i = 0; i < nbChannels; ++i) {
    fprintf(out, "%d\t%d\n", X[i], Y[i]);
  }

  fclose(out);
  Plot(nbChannels, X, Y);
}
