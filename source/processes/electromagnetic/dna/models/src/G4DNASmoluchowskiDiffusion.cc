//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/*
 * G4DNASmoluchowskiDiffusion.cc
 *
 *  Created on: 2 févr. 2015
 *      Author: matkara
 */

//#define DNADEV_TEST

#ifdef DNADEV_TEST
#include "../include/G4DNASmoluchowskiDiffusion.hh"
#else
#include "G4DNASmoluchowskiDiffusion.hh"
#endif

//#if __cplusplus >= 201103L
#ifdef DNADEV_TEST
#include "TRint.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TMath.h"
#endif

G4DNASmoluchowskiDiffusion::G4DNASmoluchowskiDiffusion(double epsilon) :  fEpsilon(epsilon)
{
  fNbins = (int) trunc(1/fEpsilon);
  // std::cout << "fNbins: " << fNbins << std::endl;
#ifdef DNADEV
  assert(fNbins > 0);
#endif
  fInverse.resize(fNbins+2); // trunc sous-estime + borne max a rajouter ==> 2

  // std::cout << "fInverse.capacity(): "<< fInverse.capacity() << std::endl;
}

G4DNASmoluchowskiDiffusion::~G4DNASmoluchowskiDiffusion()
{
}
//#endif

// --> G4DNASmoluchowskiDiffusion -- DEVELOPMENT TEST
#ifdef DNADEV_TEST

static G4DNASmoluchowskiDiffusion gDiff;

double time_test = 1e-6 /*s*/;
double D = 4.9e-9 /*m2/s*/;
double test_distance = 1e-9; // m

double Plot(double* x, double* )
{
  double diff = gDiff.GetDensityProbability(x[0], time_test, D);
  return diff;
}

static double InvErfc(double x)
{
  return TMath::ErfcInverse(x);
}

Axis_t* BinLogX(Int_t bins, Axis_t from, Axis_t to) // en puissance de 10
{
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);
//     std::cout << new_bins[i] << std::endl;
   }
   return new_bins;
}

int main(int argc, char **argv)
{
  gDiff.InitialiseInverseProbability();
//  srand (time(NULL));
  TRint* root = new TRint("G4DNASmoluchowskiDiffusion",&argc, argv);
  double interval = 1e-5;
  G4DNASmoluchowskiDiffusion* diff = new G4DNASmoluchowskiDiffusion(interval);
  diff->InitialiseInverseProbability();

//  for(size_t i = 0 ; i < diff->fInverse.size() ; ++i)
//  {
//    std::cout << i*interval << " "<< diff->fInverse[i] << std::endl;
//  }

  std::cout << diff->fInverse.size() << std::endl;

  TCanvas* canvas = new TCanvas();
  //canvas->SetLogx();
  //canvas->SetLogy();
//
//  TF1 * f = new TF1("f",diff,&G4DNASmoluchowskiDiffusion::PlotInverse,0,10,0,"G4DNASmoluchowskiDiffusion","Plot");   // create TF1 class.
//  f->SetNpx(100000);
//  f->Draw();
//  canvas->Draw();
//
//  canvas = new TCanvas();
  TH1D* h1 = new TH1D("h1", "h1", 100, 0., 1e-6);
  double distance = -1;

  int N = 100000;

  for(size_t i = 0 ; i < N ; ++i)
  {
    distance = diff->GetRandomDistance(time_test,D);
    h1->Fill(distance);
    //std::cout << distance << std::endl;
  }

  double scalf;

  {
  int integral_h1 = h1->Integral();
  h1->Scale(1./integral_h1);
  scalf=h1->GetBinWidth ( 1 ) ;
  h1->Scale(1./scalf);
  h1->GetXaxis()->SetTitle("distance");
  }

  TH1D* h2 = new TH1D("h2", "h2", 100, 0., 1e-6);
  TH1D* h_irt_distance = new TH1D("h2", "h2", 100, 0., 1e-6);

  for(size_t i = 0 ; i < N ; ++i)
  {
    double x = std::sqrt(2*D*time_test)*root_random.Gaus();
    double y = std::sqrt(2*D*time_test)*root_random.Gaus();
    double z = std::sqrt(2*D*time_test)*root_random.Gaus();

    distance = std::sqrt(x*x+y*y+z*z);
    h2->Fill(distance);
    //std::cout << distance << std::endl;

    double proba = root_random.Rndm();
    double irt_distance = InvErfc(proba)*2*std::sqrt(D*time_test);
    h_irt_distance->Fill(irt_distance);
  }

  {
  int integral_h2 = h2->Integral();
  h2->Scale(1./integral_h2);
  scalf=h2->GetBinWidth ( 1 ) ;
  h2->Scale(1./scalf);
  }

  {
  int integral_h_irt_distance = h_irt_distance->Integral();
  h_irt_distance->Scale(1./integral_h_irt_distance);
  scalf = h_irt_distance->GetBinWidth ( 1 ) ;
  h_irt_distance->Scale(1./scalf);
  h_irt_distance->GetXaxis()->SetTitle("distance");
  }


  TF1 * f2 = new TF1("f2",&Plot,0,1e-6,0,"Plot");   // create TF1 class.
  //f2->SetNpx(1000);
  h1->Draw();
  // h1->DrawNormalized();
  f2->Draw("SAME");
  h2->Draw("SAME");
  h_irt_distance->Draw("SAME");
  double integral = f2->Integral(0., 1e-6);
  std::cout << "integral = " << integral << std::endl;
  std::cout << "integral h1 = " << h1->Integral() << std::endl;
  canvas->Draw();

  std::vector<double> rdm(3);
  int nbins = 100;
  Axis_t* bins = BinLogX(nbins, -12, -1);

  TH1D* h3 = new TH1D("h3", "h3", 100, bins);
  TH1D* h4 = new TH1D("h4", "h4", 100, bins);
  TH1D* h_irt = new TH1D("h_irt", "h_irt", 100, bins);

  for(size_t i = 0 ; i < N ; ++i)
  {
    for(size_t j = 0 ; j < 3 ; ++j)
      rdm[j] = root_random.Gaus();

    double denum = 1./(rdm[0]*rdm[0] + rdm[1]*rdm[1] + rdm[2]*rdm[2]);

    double t = ((test_distance*test_distance)*denum)*1./(2*D);
    h3->Fill(t);

    double t_h4 =  diff->GetRandomTime(test_distance,D);
    h4->Fill(t_h4);
//    std::cout << t  << " " << t_h4 << std::endl;

    double proba = root_random.Rndm();
    double t_irt =  1./(4*D)*std::pow((test_distance)/InvErfc(proba),2);
    h_irt ->Fill(t_irt);
  }

  {
    TCanvas* c1 = new TCanvas();
    c1->SetLogx();
    int integral_h3 = h3->Integral();
    h3->Scale(1./integral_h3);
    scalf=h3->GetBinWidth ( 1 ) ;
    h3->Scale(1./scalf);
    h3->SetLineColor(1);
    h3->GetXaxis()->SetTitle("time");;
    h3->Draw();
  }

  {
//    TCanvas* c1 = new TCanvas();
//    c1->SetLogx();
    int integral_h4 = h4->Integral();
    h4->Scale(1./integral_h4);
    scalf=h4->GetBinWidth ( 1 ) ;
    h4->Scale(1./scalf);
    h4->SetLineColor(6);
    h4->Draw("SAME");
  //  h4->Draw("SAME");
  }

  {
//    TCanvas* c1 = new TCanvas();
//    c1->SetLogx();
    int integral_h_irt = h_irt->Integral();
    h_irt->Scale(1./integral_h_irt);
    scalf=h_irt->GetBinWidth ( 1 ) ;
    h_irt->Scale(1./scalf);
    h_irt->SetLineColor(4);
    h_irt->Draw("SAME");
  //  h4->Draw("SAME");
  }
  root->Run();
  return 0;
}
#endif
