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

const int NModels=2;

const int NModelsPiMinus=2;
// const int NModelsPiMinus=3;

//
std::string Models[2]  = { "CHIPS", "stopping" };
// std::string Models[3] = { "CHIPS", "stopping", "stopping-alt" };

int         ColorModel[3]    = {1, 2, 3};
int         SymbModel[3]     = {24, 29, 25};

const int   NPointsMadey = 30;
const int   NTargetsMadey = 7; 
float KENeut[30];  
float Value[7][30], Error[7][30];
std::string TargetsMadey[7] = { "C", "N", "O", "Al", "Cu", "Ta", "Pb" };

const int NPointsPbar_PionMom  = 44;
const int NPointsPbar_PionMult = 6;
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

   std::string fname1 = "./antiproton/pbarH_charged_pions_mom.dat";
   std::cout << "Reads data from file " << fname1 << "\n";
   
   ifstream infile;
   
   infile.open(fname1.c_str());
      
   for ( int i=0; i<NPointsPbar_PionMom; i++ )
   {
      // NOTE: the numbers in the input file came from arXiv:hep-ph/9504362v1,
      //       by the use of DigitizeIt to extract data from the graph 
      //       (the "real original" probably comes from a conf.talk in 1975...)
      infile >> MomX[i] >> MomValue[i] >> MomError[i];
   }
   infile.close();
   
   std::string fname2 = "./antiproton/pbarH_pions_mult.dat";   
   std::cout << "Reads data from file " << fname2 << "\n";
   
   infile.open(fname2.c_str());
   
   for ( int i=0; i<NPointsPbar_PionMult; i++ )
   {
      infile >> MultX[i] >> MultValue[i] >> MultError[i];
      // numbers in % - scale to ratio's
      MultValue[i] /= 100.; 
      MultError[i] /= 100.;
      // std::cout << MultX[i] << MultValue[i] << MultError[i] << std::endl;
   }
   infile.close();
      
   return;
   
}

void plotPiMinusAll()
{

   readMadey();

   TCanvas *myc = new TCanvas("myc","",1000,600);
   myc->Divide(4,2);

   for ( int i=0; i<NTargetsMadey; i++ )
   {      
      myc->cd(i+1); gPad->SetLogx(1); gPad->SetLogy(1);
      drawPiMinus( TargetsMadey[i] ); 
   }

   return ;
}

void plotPiMinus( std::string target )
{

   readMadey();

   TCanvas *myc = new TCanvas("myc","",800,600);
   myc->SetLogx(1);
   myc->SetLogy(1);
   
   drawPiMinus( target );
   
   return;   

}

void plotKMinus( std::string target )
{

   TCanvas *myc = new TCanvas("myc","",800,600);
   //myc->SetLogx(1);
   //myc->SetLogy(1);

   myc->Divide(3,2);
   
   myc->cd(1); // gPad->SetLogx(1); gPad->SetLogy(1);
   drawKMinus( target, "NSecondaries" );
   myc->cd(2); // gPad->SetLogx(1); gPad->SetLogy(1);
   drawKMinus( target, "NChargedSecondaries" );
   myc->cd(3);
   drawKMinus( target, "NNeutrons" );
   myc->cd(4);
   drawKMinus( target, "ChargeOfSecondary" );
   myc->cd(5);
   drawKMinus( target, "ChargedSecondaryMomentum" );

   return;

}

void drawPiMinus( std::string target )
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
      
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int m=0; m<NModelsPiMinus; m++ )
   {
      std::string histofile = "piminus" + TargetsMadey[TargetID] + Models[m];
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      hi[m] = (TH1F*)f->Get("NvsT");
      // turn off the stats pad; use this space to draw color codes, etc.
      hi[m]->SetStats(0);
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
   
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
   
   for ( int m=0; m<NModelsPiMinus; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*5.); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], Models[m].c_str(), "L" );
   }
      
   TGraph*  gr1 = new TGraphErrors(NPointsMadey,KENeut,Value[TargetID],0,Error[TargetID]);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   gr1->Draw("p");
   
   leg->AddEntry( gr1, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);

   return;

}

void drawKMinus( std::string target, std::string histo )
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
      
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int m=0; m<NModelsPiMinus; m++ )
   {
      std::string histofile = "kminus" + TargetsMadey[TargetID] + Models[m];
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      hi[m] = (TH1F*)f->Get( histo.c_str() );
      // turn off the stats pad; use this space to draw color codes, etc.
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      hi[m]->GetXaxis()->SetTitle( histo.c_str() );
      // hi[m]->GetYaxis()->SetTitle("Number of neutrons per MeV");
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
   
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
   
   for ( int m=0; m<NModelsPiMinus; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], Models[m].c_str(), "L" );
   }
      
   leg->Draw();
   leg->SetFillColor(kWhite);

   return ;

}


void plotAntiProton( std::string target )
{

   TH1F* hi_mom[NModels];
   TH1F* hi_mult[NModels];
   
   TCanvas *myc = new TCanvas("myc","",800,600);
   myc->Divide(2,1);
   myc->cd(1); gPad->SetLeftMargin(0.15);
   //gPad->SetLogx(1); gPad->SetLogy(1);
   myc->cd(2);gPad->SetLeftMargin(0.15); 
   //gPad->SetLogx(1); gPad->SetLogy(1);
   
   double ymin_mom = 100000., ymin_mult = 100000. ; // something big... don't know if I can use FLT_MAX
   double ymax_mom = -1., ymax_mult = -1. ;
   for ( int m=0; m<NModels; m++ )
   {
      std::string histofile = "antiproton" + target + Models[m];
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      // f->ls();
      hi_mom[m] = (TH1F*)f->Get("ChargedPionMomentum");
      hi_mult[m] = (TH1F*)f->Get("NPions");
      hi_mom[m]->SetLineColor(ColorModel[m]);
      hi_mult[m]->SetLineColor(ColorModel[m]);
      hi_mom[m]->SetLineWidth(2);
      hi_mult[m]->SetLineWidth(2);
      hi_mom[m]->GetXaxis()->SetTitle("Charged Pion Momentum (GeV/c)");
      hi_mom[m]->GetYaxis()->SetTitle("dN/dP (GeV/c)^{-1} / NEvents");
      hi_mult[m]->GetXaxis()->SetTitle("Number of pions (#pi^{+} + #pi^{-} + #pi^{0})");
      hi_mult[m]->GetYaxis()->SetTitle("Fraction [%] of events" );
      hi_mom[m]->GetYaxis()->SetTitleOffset(1.5);
      hi_mult[m]->GetYaxis()->SetTitleOffset(1.5);
      int nx = hi_mom[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi_mom[m]->GetBinContent(k);
	if ( yy > ymax_mom ) ymax_mom = yy;
	if ( yy < ymin_mom && yy > 0. ) ymin_mom = yy;
      }
      nx = hi_mult[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi_mult[m]->GetBinContent(k);
	if ( yy > ymax_mult ) ymax_mult = yy;
	if ( yy < ymin_mult && yy > 0. ) ymin_mult = yy;
      }      
      myc->cd(1);
      if ( m == 0 ) hi_mom[m]->Draw();
      else hi_mom[m]->Draw("same");
      myc->cd(2);
      if ( m == 0 ) hi_mult[m]->Draw();
      else hi_mult[m]->Draw("same");      
   }
   
   TLegend* leg1 = new TLegend(0.6, 0.70, 0.9, 0.9);
   TLegend* leg2 = new TLegend(0.6, 0.70, 0.9, 0.9);

   for ( int m=0; m<NModels; m++ )
   {
      //hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*10.);
      hi_mom[m]->GetYaxis()->SetRangeUser(ymin_mom,ymax_mom*1.4);
      hi_mom[m]->SetStats(0);
      leg1->AddEntry( hi_mom[m], Models[m].c_str(), "L" );
      hi_mult[m]->GetYaxis()->SetRangeUser(ymin_mult,ymax_mult*1.1);
      hi_mult[m]->SetStats(0);
      leg2->AddEntry( hi_mult[m], Models[m].c_str(), "L" );
   }
   
   readAntiProton();
   
   TGraph*  gr1 = new TGraphErrors(NPointsPbar_PionMom,MomX,MomValue,0,MomError);
   TGraph*  gr2 = new TGraphErrors(NPointsPbar_PionMult,MultX,MultValue,0,MultError);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   
   gr2->SetMarkerColor(4);  gr2->SetMarkerStyle(22);
   gr2->SetMarkerSize(1.6);
   
   myc->cd(1);
   gr1->Draw("p"); 
   leg1->AddEntry(gr1, "exp.data", "p");
   leg1->Draw();
   leg1->SetFillColor(kWhite);
     
   myc->cd(2);
   gr2->Draw("p");
   leg2->AddEntry(gr2, "exp.data", "p");
   leg2->Draw();
   leg2->SetFillColor(kWhite);
   
   myc->cd();

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

