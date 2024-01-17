//***********************************************************************************************************
// TomoSpectrum_HIST_proton.C
// Root command file
// Type: root TomoSpectrum_HIST_proton.C
//
// It visualizes the spectrum of protons and plots a histogram by reading StimEvent data
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

struct StimEvent
{
  uint16_t energy_keV;  // different from Pixe Event, it is in keV
  uint16_t pixelIndex;
  uint16_t sliceIndex;
  uint8_t projectionIndex;
};

void Plot(vector<double>& energies, int bin, double eMin, double eMax)
{
  gROOT->Reset();

  auto mycanvas = new TCanvas("canvas", "canvas", 800, 50, 600, 600);
  mycanvas->ToggleEventStatus();

  gPad->SetLeftMargin(0.15);

  auto hist = new TH1D("HIST", "Spectrum", bin, eMin, eMax);

  for (int i = 0; i < energies.size(); ++i) {
    hist->Fill(energies[i]);
  }

  hist->Draw();
  hist->SetTitle("TOMO Energy Spectrum");
  hist->GetXaxis()->SetTitle("ADC channels");
  hist->GetYaxis()->SetTitle("Nb events");

  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();

  // hist->GetYaxis()->SetTitleOffset(2);

  mycanvas->Print("TomoSpectrum_hist_proton.png");
}

void TomoSpectrum_HIST_proton()
{
  FILE* input = fopen("../build/StimEvent_std_Detector0_Aperture10.2.DAT", "rb");

  if (input == NULL) {
    printf("----------error for opening the input file--------------\n");
    return;
  }

  //***********************************************************************
  //**************************Selection parameters (begin)*****************
  //***********************************************************************
  const int nbProjection = 10;
  const int nbSlice = 128;
  const int nbPixel = 20;

  int projection_index_begin = 0;  // starter of the projection selected
  int projection_index_end = 0;  // end of the projection selected

  int slice_index_begin = 64;  // starter of the slice selected
  int slice_index_end = 64;  // end of the slice selected

  //********************Parameters for spectrum***************************
  int nbChannels = 4096;
  double eMin = 0;  // initialization
  double eMax = 0;  //

  //***********************************************************************
  //**************************Selection parameters (end)*******************
  //***********************************************************************

  vector<double> energies;
  StimEvent s;
  while (fread(&s, 7, 1, input)) {
    if (s.projectionIndex >= projection_index_begin && s.projectionIndex <= projection_index_end) {
      if (s.sliceIndex >= slice_index_begin && s.sliceIndex <= slice_index_end) {
        energies.push_back(s.energy_keV);
        if (s.energy_keV > eMax) eMax = s.energy_keV;
      }
    }
  }

  fclose(input);

  if (eMax > 4096) printf("---error in data----\n");
  long int size = energies.size();
  vector<int> X(nbChannels);  // save channels 1-4096
  vector<int> Y(nbChannels);

  for (long int i = 0; i < size; ++i) {
    int energy = energies[i];
    Y[energy - 1] = Y[energy - 1] + 1;
  }
  for (int i = 0; i < nbChannels; ++i) {
    X[i] = 1 + i;
  }

  FILE* out = fopen("Spectrum_hist_proton.txt", "wb");

  for (int i = 0; i < nbChannels; ++i) {
    fprintf(out, "%d\t%d\n", X[i], Y[i]);
  }

  fclose(out);

  Plot(energies, nbChannels, 0, nbChannels);
}
