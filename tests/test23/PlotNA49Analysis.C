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

const int NModels = 2;
// std::string ModelName[4]  = { "qgsp_ftfp_bert", "ftfp_bert", "ftfp", "qgsp" };  
std::string ModelName[2]  = { "NuBeam", "ftfp_bert" };  
//std::string ModelName[5]  = { "NuBeam", "qgsp_bert", "ftfp_bert", "NuBeam-with-res-decays", "qgsp-g4lund-str-fragm"};  
// std::string ModelName[4]  = { "ftfp", "qgsp", "ftfp_bert", "qgsp_ftfp_bert" };  
// const int NModels = 3;
// std::string ModelName[3]  = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };
int         ColorModel[5] = { kMagenta, 7, kRed, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue

static int isNA49UtilLoaded = 0;

void plot_pions_for_numi_talk()
{
     
   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }
      
   TCanvas* myc1 = new TCanvas("myc1", "", 1000, 600 );
   myc1->Divide(2,1);
   
   myc1->cd(1); 
   drawIntegratedSpectrum( "proton", "C", "piplus", "dNdxF" );
   myc1->cd(2); 
   drawIntegratedSpectrum( "proton", "C", "piminus", "dNdxF" );

   TCanvas* myc2 = new TCanvas("myc2", "", 1000, 600 );
   myc2->Divide(2,1);
   
   myc2->cd(1); 
   drawIntegratedSpectrum( "proton", "C", "piplus", "pT" );
   myc2->cd(2); 
   drawIntegratedSpectrum( "proton", "C", "piminus", "pT" );

   return;

}

void plot_protons_for_numi_talk()
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   TCanvas* myc = new TCanvas("myc", "", 1000, 600 );
   myc->Divide(2,1);
   
   myc->cd(1); 
   drawIntegratedSpectrum( "proton", "C", "proton", "dNdxF" );
   myc->cd(2); 
   drawIntegratedSpectrum( "proton", "C", "proton", "pT" );
   
   return;

}


void plot_dNdxF( std::string beam, std::string target )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
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

void plot_pT( std::string beam, std::string target )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
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

