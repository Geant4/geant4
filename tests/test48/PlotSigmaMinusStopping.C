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
#include "TLegend.h"


//const int NModelsSigma=2;
//std::string ModelsSigma[2] = { "CHIPS", "Bertini" };
const int NModelsSigma=1;
std::string ModelsSigma[1] = { "Bertini" };

int         ColorModel[2]    = { 6, 3 };
// int         ColorModel[3]    = { 1, 6, 3 };
// not needed, in fact...
// int         SymbModel[4]     = { 20, 29, 21, 8 };

const int NDModes = 2;
// NOTE: Error on teh Lambda mode is given at 0.06 (0.48+/-0.06)
// which corresponds to 10.7%
// Systematic error is said to be of the same level, i.e. 10-11%
//
float SigmaHyBR[2] = { 0.58, 0.42 };
float SigmaEHyBR[2] = { 0.06, 0.045 };
int ModeID[2] =  { 7, 9 };
std::string Labels[2] = { "#Lambda^{0}", "#Sigma^{-}" }

/*
const int NVersions = 4;
int ColorVersion[5] = { kRed, kGreen, 7, kBlack, 14 };
std::string Versions[4] = { "geant4-09-06-p03", "geant4-10-00-p01", "geant4-10-00-p02", "geant4-10-01-b01" };
*/

#include "../test23/shared-root-macros/REGRESSION_TEST.h" 

void plotSigmaMinus( std::string target )
{

   TCanvas* myc = new TCanvas("myc1","",800,600);
   // myc->Divide(2,1);
   
   // myc->cd(1);
   gPad->SetLeftMargin(0.15);
   gPad->SetLogy();
   drawSigmaMinus( target, "PartTypeMult" );
      
   myc->cd();
   
   myc->Print("sigmamin-H-models.gif");

   return;

}

void plotSigmaMinusRegre( std::string target, std::string model )
{

   TCanvas* myc = new TCanvas("myc1","",800,600);
   // myc->Divide(2,1);
   
   // myc->cd(1);
   gPad->SetLeftMargin(0.15);
   gPad->SetLogy();
   drawSigmaMinusRegre( target, "PartTypeMult", model );
      
   myc->cd();
   myc->Print( "sigmamin-H-Bert-regre.gif" );

   return;

}



void drawSigmaMinus( std::string target, std::string histo, 
                     bool useStatEr=true, bool useSysEr=true )
{

   TH1F** hi = new TH1F*[NModelsSigma];
      
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
      
   for ( int m=0; m<NModelsSigma; m++ )
   {
      std::string histofile = "sigma-" + target + ModelsSigma[m];
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      hi[m] = (TH1F*)f->Get( histo.c_str() );
      // turn off the stats pad; use this space to draw color codes, etc.
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      if ( histo == "PartTypeMult" )
      {
         hi[m]->GetXaxis()->SetTitle("secondary particle type" );
         hi[m]->GetYaxis()->SetTitle("Production rate");
	 hi[m]->GetYaxis()->SetTitleOffset(1.5);
      }
      else
      {
         hi[m]->GetXaxis()->SetTitle( histo.c_str() );
      }
//      int nx = hi[m]->GetNbinsX();
//      for (int k=1; k <= nx; k++) {
//	double yy = hi[m]->GetBinContent(k);
//	if ( yy > ymax ) ymax = yy;
//	if ( yy < ymin && yy > 0. ) ymin = yy;
//      }
      if ( m == 0 ) hi[m]->Draw("histo");
      else hi[m]->Draw("histosame");
   }
   
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
   
   for ( int m=0; m<NModelsSigma; m++ )
   {
//      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelsSigma[m].c_str(), "L" );
   }
      
   if ( histo == "PartTypeMult" && target == "H" )
   {
      TH1F* h = new TH1F("h", " ", 15, 0, 15. );

      for ( int i=0; i<NDModes; i++ )
      {
	 //float val = (float)ModeID[i] + 0.5;
	 //h->Fill( val, SigmaHyBR[i] );
	 h->SetBinContent( ModeID[i]+1, SigmaHyBR[i] );
	 float totalerr2 = 0.;
	 if ( useStatEr ) totalerr2 += SigmaEHyBR[i]*SigmaEHyBR[i];
	 if ( useSysEr )  totalerr2 *= 2; // since sys err is said to be 
	                                  // the same as stat we can assume
					  // 50% (factor of 2 ) increase 
	 h->SetBinError( ModeID[i]+1, sqrt(totalerr2) );
	 h->GetXaxis()->SetBinLabel( ModeID[i], Labels[i].c_str() );
      }
      h->SetStats(0);
      h->SetMarkerStyle(22);
      h->SetMarkerColor(4);
      h->SetMarkerSize(1.6);
      h->Draw("psame");
      leg->AddEntry( h, "exp.data", "p" );

//      
//      TGraph* gr1 = new TGraph( NDModes, ModeID, KaonicHyBR );
//      gr1->SetLineColor(4);
//      gr1->SetMarkerStyle(22);
//      gr1->SetMarkerSize(1.6);
//      gr1->Draw("p");
//      leg->AddEntry( gr1, "exp.data", "p" ); 
//
   }
   
      
   leg->Draw();
   leg->SetFillColor(kWhite);

   return ;

}

void drawSigmaMinusRegre( std::string target, std::string histo, std::string model )
{

   TH1F** hi = new TH1F*[NVersions];
      
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
      
   for ( int m=0; m<NVersions; m++ )
   {
//      std::string histofile = Versions[m] + "/sigma-" + target + model;
//      histofile += ".root";

      std::string location = "";
      if ( Versions[m] == CurrentVersion || Versions[m] == "." )
      {
         location = "";
      }
      else
      {
         location = regre_test_dir + "/test48/" + Versions[m] + "/";
      }
      std::string histofile = location + "sigma-" + target + model + ".root"; 

      TFile* f = new TFile( histofile.c_str() );
      
      if ( f ) std::cout << " histofile " << histofile << " is open" << std::endl;
      
      hi[m] = (TH1F*)f->Get( histo.c_str() );
      // turn off the stats pad; use this space to draw color codes, etc.
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorVersion[m]);
      hi[m]->SetLineWidth(6-m);
      if ( histo == "PartTypeMult" )
      {
         hi[m]->GetXaxis()->SetTitle("secondary particle type" );
         hi[m]->GetYaxis()->SetTitle("Production rate");
	 hi[m]->GetYaxis()->SetTitleOffset(1.5);
      }
      else
      {
         hi[m]->GetXaxis()->SetTitle( histo.c_str() );
      }
//      int nx = hi[m]->GetNbinsX();
//      for (int k=1; k <= nx; k++) {
//	double yy = hi[m]->GetBinContent(k);
//	if ( yy > ymax ) ymax = yy;
//	if ( yy < ymin && yy > 0. ) ymin = yy;
//      }
      if ( m == 0 ) hi[m]->Draw();
      else hi[m]->Draw("same");
   }
   
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
   
   for ( int m=0; m<NVersions; m++ )
   {
//      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], Versions[m].c_str(), "L" );
   }
      
   if ( histo == "PartTypeMult" && target == "H" )
   {
      TH1F* h = new TH1F("h", " ", 15, 0, 15. );

      for ( int i=0; i<NDModes; i++ )
      {
	 float val = (float)ModeID[i] + 0.5;
	 h->Fill( val, SigmaHyBR[i] );
	 h->GetXaxis()->SetBinLabel( ModeID[i], Labels[i].c_str() );
      }
      h->SetStats(0);
      h->SetMarkerStyle(22);
      h->SetMarkerColor(4);
      h->SetMarkerSize(1.6);
      h->Draw("psame");
      leg->AddEntry( h, "exp.data", "p" );

//      
//      TGraph* gr1 = new TGraph( NDModes, ModeID, KaonicHyBR );
//      gr1->SetLineColor(4);
//      gr1->SetMarkerStyle(22);
//      gr1->SetMarkerSize(1.6);
//      gr1->Draw("p");
//      leg->AddEntry( gr1, "exp.data", "p" ); 
//
   }   
      
   leg->Draw();
   leg->SetFillColor(kWhite);

   return ;

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
