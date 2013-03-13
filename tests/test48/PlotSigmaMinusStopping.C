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
float SigmaHyBR[2] = { 0.58, 0.42 };
int ModeID[2] =  { 7, 9 };
std::string Labels[2] = { "#Lambda^{0}", "#Sigma^{-}" }


const int NVersions = 1;
int ColorVersion[3] = { kGreen, kRed, kBlack };
//std::string Versions[2] = { "geant4-09-04-ref10", "geant4-09-05-ref01" };
std::string Versions[1] = { "geant4-09-05-ref02" };

void plotSigmaMinus( std::string target )
{

   TCanvas* myc = new TCanvas("myc1","",800,600);
   // myc->Divide(2,1);
   
   // myc->cd(1);
   gPad->SetLeftMargin(0.15);
   gPad->SetLogy();
   drawSigmaMinus( target, "PartTypeMult" );
      
   myc->cd();

   return;

}

void drawSigmaMinus( std::string target, std::string histo )
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
      if ( m == 0 ) hi[m]->Draw();
      else hi[m]->Draw("same");
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
