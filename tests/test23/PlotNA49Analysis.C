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

const int NModels = 3;
// std::string ModelName[4]  = { "qgsp_ftfp_bert", "ftfp_bert", "ftfp", "qgsp" };  
std::string ModelName[2]  = { "NuBeam-with-decays", "ftfp_bert-with-decays", "qgsp_bert-with-decays" };  
//std::string ModelName[5]  = { "NuBeam", "qgsp_bert", "ftfp_bert", "NuBeam-with-res-decays", "qgsp-g4lund-str-fragm"};  
// std::string ModelName[4]  = { "ftfp", "qgsp", "ftfp_bert", "qgsp_ftfp_bert" };  
// const int NModels = 3;
// std::string ModelName[3]  = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };
int         ColorModel[5] = { kMagenta, 7, kRed, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue
int         SymbModel[4]     = { 8, 21, 23, 25 };

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

void plot_dNdxF_pT( std::string beam, std::string target, std::string secondary )
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
   drawIntegratedSpectrum( beam, target, secondary, "dNdxF" );
   pad1->cd(2); gPad->SetRightMargin(0.025);
   drawIntSpectrumMC2Data( beam, target, secondary, "dNdxF" );
   
   myc->cd();
      
   TPad* pad2 = new TPad( "pad2", "", 0.51, 0.01, 0.99, 0.99 );

   pad2->Draw();
   pad2->Divide(1.,2.,0.,0.);
   pad2->cd(1); gPad->SetRightMargin(0.025);
   drawIntegratedSpectrum( beam, target, secondary, "pT" );
   pad2->cd(2); gPad->SetRightMargin(0.025);
   drawIntSpectrumMC2Data( beam, target, secondary, "pT" );
   
   return;

} 

void plotDDiffXSec( std::string beam, std::string target, std::string secondary, int start=0, int end=22 )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   readDDiffSpectra( beam, target, secondary );

   int N1 = std::min(0,start);
   int N2 = std::max(end,NSubSets); 
   int NN = N2 - N1 ;
   int NN1 = 0;

   
   TCanvas** cnv = 0;
   if ( NN%2 == 0 )
   {
      NN1 = NN/2;
   }
   else
   {
      NN1 = NN/2 + 1;
   }

   cnv = new TCanvas*[NN1];

   for ( int i=0; i<NN1; ++i )
   {
      std::ostringstream cnt;
      cnt << i;
      std::string name = "cnv" + cnt.str();
      cnv[i] = new TCanvas( name.c_str(), "", 800, 500 );
      cnv[i]->Divide( 2, 1 );
   }
   
   for ( int i=0; i<NN; ++i )
   {
      cnv[i/2]->cd((i%2)+1);
      draw1DDiffXSec( beam, target, secondary, i );
   }

   return;

}

