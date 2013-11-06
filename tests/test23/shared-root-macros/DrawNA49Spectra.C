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

//const int NModels = 2;
// std::string ModelName[4]  = { "qgsp_ftfp_bert", "ftfp_bert", "ftfp", "qgsp" };  
//std::string ModelName[2]  = { "NuBeam", "ftfp_bert" };  
// std::string ModelName[5]  = { "NuBeam", "qgsp_bert", "ftfp_bert", "NuBeam-with-res-decays", "qgsp-g4lund-str-fragm"};  
// std::string ModelName[4]  = { "ftfp", "qgsp", "ftfp_bert", "qgsp_ftfp_bert" };  
// const int NModels = 3;
// std::string ModelName[3]  = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };
//int         ColorModel[5] = { kMagenta, 7, kRed, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue

// static int isReadNA49Loaded = 0;


void drawIntegratedSpectrum( std::string beam, std::string target,  
                             std::string secondary, std::string histo )
{

//   std::cout << "About to load Reader" << std::endl;
   
//   if ( isReadNA49Loaded <= 0 )
//   {
//      gROOT->LoadMacro("ReadNA49Data.C");
//      isReadNA49Loaded = 1;
//   }
   
//   std::cout << "Reader loaded" << std::endl;
 
   readIntegratedSpectra( beam, target, secondary );

   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;

   float* Value = new float[NPoints];
   float* Error = new float[NPoints];
   
   for ( int i=0; i<NPoints; i++ )
   {
      if ( histo == "dNdxF" )
      {
         Value[i] = dNdxF[i];
	 Error[i] = err_dNdxF[i];
      }
      else if ( histo == "pT" )
      {
         Value[i] = pT[i];
	 Error[i] = err_pT[i];
      }
      else if ( histo == "pT2" )
      {
         Value[i] = pT2[i];
	 Error[i] = err_pT2[i];
      }
      if ( Value[i]+Error[i] > ymax ) ymax = Value[i] + Error[i];
      if ( Value[i]-Error[i] < ymin ) ymin = Value[i] - Error[i];
      if ( ymin < 0. ) ymin = 0.; 
   }
 
   TH1F* hi[NModels];
   std::string YTitle;
   if ( histo == "dNdxF" )
   {
      YTitle = "dN/dxF";
   }
   else if ( histo == "pT" )
   {
      YTitle = "dpT/dxF, GeV/c";
   }
   else if ( histo == "pT2" )
   {
      YTitle = "dpT^2/dxF, (GeV/c)^2";
   }
   
   for ( int m=0; m<NModels; m++ )
   {

      std::string histofile = "./na49-histo/" + beam + target + "158.0GeV" + ModelName[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_" + histo ;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
/*
      if ( histo == "pT" && m<3 )
      {
         hi[m]->Scale(32.);
      }
*/
      
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      if ( m == 0 ) hi[m]->SetLineWidth(3.5);
      
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) 
      {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
         hi[m]->Draw();
	 hi[m]->GetXaxis()->SetTitle("xF");
	 hi[m]->GetYaxis()->SetTitle( YTitle.c_str() );
	 hi[m]->GetYaxis()->SetTitleOffset(1.5);
      }
      else hi[m]->Draw("same");     

   }
   
   TLegend* leg = new TLegend(0.55, 0.65, 0.9, 0.9);
   leg->SetTextSize(0.025);
   
   for ( int m=0; m<NModels; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelName[m].c_str(), "L" );
   }
     
   TGraph* gr = new TGraphErrors( NPoints, xF, Value, 0, Error );

   // gr->GetYaxis()->SetRangeUser( 0., 2.5 );
   // gr->GetXaxis()->SetRangeUser( -0.3, 0.4 );
   // gr->SetRangeUser( 0., 2.5 );
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
    
   gr->Draw("p");
      
   leg->AddEntry( gr, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);   
   
   return;

}

