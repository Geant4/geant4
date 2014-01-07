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

// for all types of secondaries except pions
//

const int NModels = 3;
// std::string ModelName[4]  = { "pythia8180", "ftfp-with-decays", "qgsp-with-decays", "qgsp-g4lund-str-fragm-with-decays" };
// std::string ModelName[3]  = { "ftfp-with-decays", "qgsp-with-decays", "qgsp-g4lund-str-fragm-with-decays" };
std::string ModelName[3]  = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };
int         ColorModel[5] = { kBlack, kGreen, kRed, 14, 7 }; // 14 = grey, 7 = light "sky"-blue
//
//int         ColorModel[4] = { kGreen, 7, kRed, kBlack }; // 14 = grey, 7 = light "sky"-blue
// int         ColorModel[4] = { kMagenta, 7, kRed, kBlack }; // 14 = grey, 7 = light "sky"-blue
int         SymbModel[4]     = { 8, 21, 23, 25 };

const int NVersions = 4;
// ---> std::string Versions[3] = { "geant4-10-00-b01", "geant4-09-06-ref03", "geant4-09-06-p01" };
std::string Versions[4] = { "geant4-09-06-ref07", "geant4-10-00-b01", "geant4-09-06-ref03", "geant4-09-06-p01" };
// std::string Versions[2] = { "geant4-10-00-b01", "geant4-09-06-ref03" };
std::string CurrentVersion = "geant4-10-00-b01";
int ColorVersion[5] = { kRed, kGreen, 7, kBlack, 14 };


static int isNA49UtilLoaded = 0;


void plot_dNdxF( std::string beam, std::string target )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   TCanvas *myc = new TCanvas("myc","",1200,800);
   myc->Divide(3,2);

   myc->cd(1);
   drawIntegratedSpectrum( beam, target, "proton", "dNdxF" );
   
   myc->cd(2);
   drawIntegratedSpectrum( beam, target, "antiproton", "dNdxF" );

   myc->cd(3);
   drawIntegratedSpectrum( beam, target, "neutron", "dNdxF" );

   myc->cd(4);
   drawIntegratedSpectrum( beam, target, "piplus", "dNdxF" );

   myc->cd(5);
   drawIntegratedSpectrum( beam, target, "piminus", "dNdxF" );


   return;

}

void plot_dNdxF_Regre( std::string beam, std::string target, std::string model )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   TCanvas *myc = new TCanvas("myc","",1200,800);
   myc->Divide(3,2);
   
   myc->cd(1);
   drawIntergratedSpectrumRegre( beam, target, "proton", "dNdxF", model );

   myc->cd(2);
   drawIntergratedSpectrumRegre( beam, target, "antiproton", "dNdxF", model );

   myc->cd(3);
   drawIntergratedSpectrumRegre( beam, target, "neutron", "dNdxF", model );

   myc->cd(4);
   drawIntergratedSpectrumRegre( beam, target, "piplus", "dNdxF", model );

   myc->cd(5);
   drawIntergratedSpectrumRegre( beam, target, "piminus", "dNdxF", model );
   
   return;

}

void plot_dNdxF_MC2Data( std::string beam, std::string target )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }
   
   TCanvas* myc = new TCanvas("myc","",1200,800);
   myc->Divide(3,2);
   
   myc->cd(1);
   drawIntSpectrumMC2Data( beam, target, "proton", "dNdxF" );
   
   myc->cd(2);
   drawIntSpectrumMC2Data( beam, target, "antiproton", "dNdxF" );
   
   myc->cd(3);
   drawIntSpectrumMC2Data( beam, target, "neutron", "dNdxF" );

   myc->cd(4);
   drawIntSpectrumMC2Data( beam, target, "piplus", "dNdxF" );

   myc->cd(5);
   drawIntSpectrumMC2Data( beam, target, "piminus", "dNdxF" );

   return ;

} 

void plot_dNdxF_pions_combined( std::string beam, std::string target )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   TCanvas* myc   = new TCanvas("myc","",1000,600);
   
/*
   TPad*    pad11 = new TPad( "pad11", "", 0.01, 0.71, 0.33, 0.99 );
   TPad*    pad12 = new TPad( "pad12", "", 0.01, 0.51, 0.33, 0.709 );
   pad11->Draw();
   pad12->Draw();
   pad11->cd();
   drawIntegratedSpectrum( beam, target, "piplus", "dNdxF" );
   pad12->cd();
   drawIntSpectrumMC2Data( beam, target, "piplus", "dNdxF" );
*/

   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.01, 0.49, 0.99 );
   
   pad1->Draw();
   pad1->Divide(1.,2.,0.,0.);
   pad1->cd(1); gPad->SetRightMargin(0.025);
   drawIntegratedSpectrum( beam, target, "piplus", "dNdxF" );
   pad1->cd(2); gPad->SetRightMargin(0.025);
   drawIntSpectrumMC2Data( beam, target, "piplus", "dNdxF" );
   
   myc->cd();
      
   TPad* pad2 = new TPad( "pad2", "", 0.51, 0.01, 0.99, 0.99 );

   pad2->Draw();
   pad2->Divide(1.,2.,0.,0.);
   pad2->cd(1); gPad->SetRightMargin(0.025);
   drawIntegratedSpectrum( beam, target, "piminus", "dNdxF" );
   pad2->cd(2); gPad->SetRightMargin(0.025);
   drawIntSpectrumMC2Data( beam, target, "piminus", "dNdxF" );
   
   return;

} 

void plot_std_diff_graph( std::string beam, std::string target )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }
   
   TCanvas* myc1 = new TCanvas("myc1","",1000,600);
   myc1->Divide(2,1);
   
   myc1->cd(1); // gPad->SetLeftMargin(0.15);
   drawStdDiffGraph( beam, target, "piplus", "dNdxF" );
   
   myc1->cd(2); // gPad->SetLeftMargin(0.15);
   drawStdDiffGraph( beam, target, "piminus", "dNdxF" );

   TCanvas* myc2 = new TCanvas("myc2","",1000,600);
   myc2->Divide(2,1);
   
   myc2->cd(1); // gPad->SetLeftMargin(0.15);
   drawStdDiffGraph( beam, target, "piplus", "pT" );
   
   myc2->cd(2); // gPad->SetLeftMargin(0.15);
   drawStdDiffGraph( beam, target, "piminus", "pT" );
   
   return;

}

void plot_std_diff_histo( std::string beam, std::string target )
{


   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }
   
   TCanvas* myc1 = new TCanvas("myc1","",1000,600);
   myc1->Divide(2,1);

   myc1->cd(1); // gPad->SetLeftMargin(0.15);
   drawStdDiffHisto( beam, target, "piplus", "dNdxF" );
   
   myc1->cd(2); // gPad->SetLeftMargin(0.15);
   drawStdDiffHisto( beam, target, "piminus", "dNdxF" );

   TCanvas* myc2 = new TCanvas("myc2","",1000,600);
   myc2->Divide(2,1);

   myc2->cd(1); // gPad->SetLeftMargin(0.15);
   drawStdDiffHisto( beam, target, "piplus", "pT" );
   
   myc2->cd(2); // gPad->SetLeftMargin(0.15);
   drawStdDiffHisto( beam, target, "piminus", "pT" );
   
   return;

}

void plot_pT( std::string beam, std::string target )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   TCanvas *myc = new TCanvas("myc","",800,800);
   myc->Divide(2,2);

   myc->cd(1);
   drawIntegratedSpectrum( beam, target, "proton", "pT" );
   
   myc->cd(2);
   drawIntegratedSpectrum( beam, target, "antiproton", "pT" );

   myc->cd(3);
   drawIntegratedSpectrum( beam, target, "piplus", "pT" );

   myc->cd(4);
   drawIntegratedSpectrum( beam, target, "piminus", "pT" );


   return;

}

void plot_pT_Regre( std::string beam, std::string target, std::string model )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   TCanvas *myc = new TCanvas("myc","",800,800);
   myc->Divide(2,2);
   
   myc->cd(1);
   drawIntergratedSpectrumRegre( beam, target, "proton", "pT", model );

   myc->cd(2);
   drawIntergratedSpectrumRegre( beam, target, "antiproton", "pT", model );

   myc->cd(3);
   drawIntergratedSpectrumRegre( beam, target, "piplus", "pT", model );

   myc->cd(4);
   drawIntergratedSpectrumRegre( beam, target, "piminus", "pT", model );
   
   return;

} 


