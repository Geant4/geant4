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
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TRefArray.h"
#include "TStyle.h"
#include "TGraph.h"

const int NModels=2;
const int NModelsPiMinus=3;

//
// std::string Models[2]  = { "CHIPS", "stopping" };
std::string Models[3] = { "CHIPS", "stopping", "stopping-alt" };

int         ColorModel[3]    = {1, 2, 3};
int         SymbModel[3]     = {24, 29, 25};

const int   NPointsMadey = 30;
const int   NTargetsMadey = 7; 
float KENeut[30];  
float Value[7][30], Error[7][30];
std::string TargetsMadey[7] = { "C", "N", "O", "Al", "Cu", "Ta", "Pb" };

const int NPointsAntiProtonMom  = 44;
const int NPointsAntiProtonMult = 6;
float MomX[44], MomValue[44], MomError[44];
float MultX[6], MultValue[6], MultError[6];

 
void readMadey()
{

   std::string fname = "./piminus/madey_spectra.dat";
   std::cout << "Reads data from file " << fname << "\n";
   ifstream infile;
   infile.open(fname.c_str());
      
   for ( int i=0; i<NPointsMadey; i++ )
   {

      infile >> KENeut[i];
      for ( int j=0; j<NTargetsMadey; j++ )
      {
         infile >> Value[j][i] >> Error[j][i];
	 // rescale as numbers are per 100 pi-'s
	 Value[j][i] /= 100. ;
	 Error[j][i] /= 100. ;
      }
   }
      
   return;

}

void readAntiProton()
{

   std::string fname1 = "./antiproton/data_apc_h.dat";
   std::cout << "Reads data from file " << fname1 << "\n";
   ifstream infile;
   infile.open(fname1.c_str());
      
   for ( int i=0; i<NPointsAntiProtonMom; i++ )
   {
      // NOTE: the numbers in the input file came from arXiv:hep-ph/9504362v1,
      //       by the use of DigitizeIt to extract data from the graph 
      //       (the "real original" probably comes from a conf.talk in 1975...)
      infile >> MomX[i] >> MomValue[i] >> MomError[i];
   }
      
   return;
   
}

void plotPiMinus( std::string target )
{
   
   int TargetID = -1;
   
   for ( int i=0; i<NTargetsMadey; i++ )
   {
      if ( TargetsMadey[i] == target ) 
      {
         TargetID = i;
	 break;
      }
   }
      
   if ( TargetID == -1 || TargetID >= NTargetsMadey )
   {
      std::cout << " Invalid Target: " << target << std::endl;
      return;
   }
   
   TH1F* hi[NModelsPiMinus];
   
   TCanvas *myc = new TCanvas("myc","",800,600);
   myc->SetLogx(1);
   myc->SetLogy(1);
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int m=0; m<NModelsPiMinus; m++ )
   {
      std::string histofile = "piminus" + TargetsMadey[TargetID] + Models[m];
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      hi[m] = (TH1F*)f->Get("NvsT");
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      hi[m]->GetXaxis()->SetTitle("Kinetic energy of secondary neutron (MeV)");
      hi[m]->GetYaxis()->SetTitle("Number of neutrons per MeV");
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) hi[m]->Draw();
      else hi[m]->Draw("same");
   }
   
   // std::cout << " ymin= " << ymin << " ymax= " << ymax << std::endl;
   
   for ( int m=0; m<NModelsPiMinus; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*5.); // hi[m]->SetTitle("");
   }

   readMadey();
   
   TGraph*  gr1 = new TGraphErrors(NPointsMadey,KENeut,Value[TargetID],0,Error[TargetID]);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   gr1->Draw("p");

   return;

}

void plotAntiProton( std::string target )
{

   TH1F* hi[NModels];
   
   TCanvas *myc = new TCanvas("myc","",800,600);
   //myc->SetLogx(1);
   //myc->SetLogy(1);
   
   double ymin = 100000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int m=0; m<NModels; m++ )
   {
      std::string histofile = "antiproton" + target + Models[m];
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      // f->ls();
      hi[m] = (TH1F*)f->Get("ChargedPionMomentum");
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      hi[m]->GetXaxis()->SetTitle("Charged Pion Momentum (GeV/c)");
      hi[m]->GetYaxis()->SetTitle("dN/dP (GeV/c)^{-1}");
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) hi[m]->Draw();
      else hi[m]->Draw("same");
   }

   std::cout << " ymin= " << ymin << " ymax= " << ymax << std::endl;
   
   for ( int m=0; m<NModels; m++ )
   {
      //hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*10.);
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.1);
   }
   
   readAntiProton();
   
   TGraph*  gr1 = new TGraphErrors(NPointsAntiProtonMom,MomX,MomValue,0,MomError);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   gr1->Draw("p");   

   return;

}

void setStyle() 
{

  gStyle->SetCanvasBorderMode(0); gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);    gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);   gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);   gStyle->SetFrameLineWidth(1);
  gStyle->SetTitleOffset(1.6,"Y");  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(1);
  
  return;

}

