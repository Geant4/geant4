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
int    NPoints         = 0;

float* xF              = 0;

float* dXSecdxF     = 0;
float* err_dXSecdxF = 0;
float* dNdxF           = 0;
float* err_dNdxF       = 0;
float* pT              = 0;
float* err_pT          = 0;
float* pT2             = 0;
float* err_pT2         = 0;

float* XSec            = 0;
float* err_XSec        = 0;

// for pions (pi+ & pi-)
//
float* Y =0;
float* dNdY[2]     = { 0, 0 };
float* err_dNdY[2] = { 0, 0 };


const int NModels = 3;
// std::string ModelName[4]  = { "qgsp_ftfp_bert", "ftfp_bert", "ftfp", "qgsp" };  
std::string ModelName[5]  = { "NuBeam", "qgsp_bert", "ftfp_bert", "NuBeam-with-res-decays", "qgsp-g4lund-str-fragm"};  
// std::string ModelName[4]  = { "ftfp", "qgsp", "ftfp_bert", "qgsp_ftfp_bert" };  
// const int NModels = 3;
// std::string ModelName[3]  = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };
int         ColorModel[5] = { kMagenta, 7, kRed, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue

void readIntegratedSpectra( std::string beam, std::string target, std::string secondary )
{

   std::string dirname = "./na49-exp-data/";
   
   std::string filename = beam + "_" + target + "_" + secondary;
   
   std::string filenam1;
   
   if ( secondary == "pion" )
   { 
      filenam1 = filename + "_dndy.dat";
   }
   else
   {
      filename1 = filename + "_integr.dat";
   }
   
   std::string file1 = dirname + filename1;
   
   ifstream infile1;
   infile1.open( file1.c_str() );

   infile1 >> NPoints;
   
   if ( secondary == "pion" )
   {
      if (Y) delete Y;
      Y = new float[NPoints];
      if ( dNdY[0] = 0 ) delete dNdY[0];
      if ( dNdY[1] = 0 ) delete dNdY[1];
      dNdY[0] = new float[NPoints]; // pi+
      dNdY[1] = new float[NPoints]; // pi-
      if ( err_dNdY[0] = 0 ) delete err_dNdY[0];
      if ( err_dNdY[1] = 0 ) delete err_dNdY[1];      
      err_dNdY[0] = new float[NPoints];
      err_dNdY[1] = new float[NPoints];
      for ( int i=0; i<NPoints; i++ )
      {
         infile1 >> Y[i] >> dNdY[0][i] >> dNdY[1][i];
	 err_dNdY[0][i] = dNdY[0][i]*0.02; // errors are said to be 2%
	 err_dNdY[1][i] = dNdY[1][i]*0.02;
      } 
      return;     
   }
   

   if (xF) delete [] xF;
   xF = new float[NPoints];
   
  
   if (dNdxF) delete dNdxF;
   dNdxF = new float[NPoints];
   
   if (err_dNdxF) delete err_dNdxF;
   err_dNdxF = new float[NPoints];
   
   if ( secondary != "neutron" )
   {
     if  (dXSecdxF) delete [] dXSecdxF;
     dXSecdxF = new float[NPoints];
     if  (err_dXSecdxF) delete [] err_dXSecdxF;
     err_dXSecdxF = new float[NPoints];
     if (pT) delete [] pT;
     pT = new float[NPoints];
     if ( err_pT ) delete [] err_pT;
     err_pT = new float[NPoints];
     if (pT2) delete [] pT2;
     pT2 = new float[NPoints];
     if ( err_pT2 ) delete [] err_pT2;
     err_pT2 = new float[NPoints];
   }

   for ( int i=0; i<NPoints; i++ )
   {
      if ( secondary == "neutron" )
      {
         infile1 >> xF[i] >> dNdxF[i] >> err_dNdxF[i];
	 err_dNdxF[i] *= 0.01;
      }
      else
      {
         infile1 >> xF[i] >> dXSecdxF[i] >> err_dXSecdxF[i]  >> dNdxF[i] >> err_dNdxF[i]  
	         >> pT[i] >> err_pT[i] >> pT2[i] >> err_pT2[i];
         err_dXSecdxF[i]    *= 0.01;
	 err_dNdxF[i]       *= 0.01;
	 err_pT[i]          *= 0.01;
	 err_pT2[i]         *= 0.01;
      }
   }
   
   return;

}

void plot_pions_for_numi_talk()
{

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


void drawIntegratedSpectrum( std::string beam, std::string target,  
                             std::string secondary, std::string histo )
{

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

      std::string histofile = "./na49-histo-no-res-decays/" + beam + target + "158.0GeV" + ModelName[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_" + histo ;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
      if ( histo == "pT" && m<3 )
      {
         hi[m]->Scale(32.);
      }
      
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

