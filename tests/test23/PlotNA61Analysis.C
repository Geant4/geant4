#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <list>

#include <math.h>
#include <vector>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"

/*
int    NPoints    = 0;

float* SecMom     = 0;

float* SecSigma   = 0; 
float* SecESigma  = 0;

float* K2PiRatio  = 0;
float* K2PiERatio = 0;
*/

const int NModels = 3;
std::string ModelName[3] = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };
// int         ColorModel[2]    = { 14, 7 }; // 14 = grey, 7 = light "sky"-blue
int         ColorModel[3] = { kGreen, 7, kRed }; // 14 = grey, 7 = light "sky"-blue

const int NVersions = 4;
std::string Versions[4] = { "geant4-09-06-ref07", "geant4-10-00-b01", "geant4-09-06-ref03", "geant4-09-06-p01" };
// std::string Versions[2] = { "geant4-10-00-b01", "geant4-09-06-ref03" };
std::string CurrentVersion = "geant4-10-00-b01";
int ColorVersion[5] = { kRed, kGreen, 7, kBlack, 14 };

static int isNA61UtilLoaded = 0;


void plotSecondarySum8( std::string secondary )
{

   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA61Data.C");
      gROOT->LoadMacro("./shared-root-macros/DrawNA61Spectra.C");
      isNA61UtilLoaded = 1;
   }

   TCanvas *myc1 = new TCanvas("myc1","",900,900);
   myc1->Divide(2,2);
   
   myc1->cd(1);
   gPad->SetLogy(); //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "0", "20" );

   myc1->cd(2);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "20", "40" );

   myc1->cd(3);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "40", "60" );

   myc1->cd(4);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "60", "100" );

   TCanvas *myc2 = new TCanvas("myc2","",900,900);
   myc2->Divide(2,2);
   
   myc2->cd(1);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "100", "140" );

   myc2->cd(2);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "140", "180" );

   myc2->cd(3);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "180", "240" );
   
   if ( secondary == "proton" ) return;

   myc2->cd(4);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "240", "300" );

   return;

}

void plotSecondarySum8Regre( std::string secondary, std::string model )
{

   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA61Data.C");
      gROOT->LoadMacro("./shared-root-macros/DrawNA61Spectra.C");
      isNA61UtilLoaded = 1;
   }

   TCanvas *myc1 = new TCanvas("myc1","",900,900);
   myc1->Divide(2,2);
   
   myc1->cd(1);
   gPad->SetLogy(); //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "0", "20", model );

   myc1->cd(2);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "20", "40", model );

   myc1->cd(3);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "40", "60", model );

   myc1->cd(4);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "60", "100", model );

   TCanvas *myc2 = new TCanvas("myc2","",900,900);
   myc2->Divide(2,2);
   
   myc2->cd(1);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "100", "140", model );

   myc2->cd(2);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "140", "180", model );

   myc2->cd(3);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "180", "240", model );
   
   if ( secondary == "proton" ) return;

   myc2->cd(4);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "240", "300", model );

   return;

}

