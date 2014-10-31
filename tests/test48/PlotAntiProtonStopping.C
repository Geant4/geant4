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

#include "../test23/shared-root-macros/REGRESSION_TEST.h"
#include "../test23/shared-root-macros/Chi2Calc.C"

// provision for near-future use of FTF for baryons
//
/*
const int NModelsBaryons=3;
std::string ModelsBaryons[3] = { "stopping", "CHIPS", "FTF" };
int         ColorModel[4]    = { 2, 1, 6, 3 };
int         SymbModel[4]     = { 29, 20, 21, 8 };
*/

#ifndef G4VAL_PLOTANTIPROTONSTOPPING_C
#define G4VAL_PLOTANTIPROTONSTOPPING_C

const int NModelsBaryons=1;
std::string ModelsBaryons[1] = { "FTF" };
int         ColorModel[1]    = { 6 };


// pbar beam business
const int NPointsPbar_PionMom  = 44;
const int NPointsPbar_PionMult = 6;
float MomX[44], MomValue[44], MomError[44];
float MultX[6], MultValue[6], MultError[6];


/* 
void PlotAntiProtonStopping()
{

  plotAntiProton("H");
  return; 

}
*/
/*
void PlotAntiProtonStoppingRegre()
{

   drawAntiProtonRegression( "H", "FTF" );
   return;
}
*/

void readAntiProton()
{

   std::string fname1 = "./antiproton/pbarH_charged_pions_mom.dat";
   std::cout << "Reads data from file " << fname1 << "\n";
   
   ifstream infile;
   
   infile.open(fname1.c_str());
      
   for ( int i=0; i<NPointsPbar_PionMom; i++ )
   {      
      MomX[i] = 0.;
      MomValue[i] = 0.;
      MomError[i] = 0.;
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
      MultX[i] = 0.;
      MultValue[i] = 0.;
      MultError[i] = 0.;
      infile >> MultX[i] >> MultValue[i] >> MultError[i];
      // numbers in % - scale to ratio's
      // MultX[i] += 0.5;
      MultValue[i] /= 100.; 
      MultError[i] /= 100.;
      // std::cout << MultX[i] << MultValue[i] << MultError[i] << std::endl;
   }
   infile.close();
      
   return;
   
}

TGraphErrors* getPionMultAsGraph()
{

   readAntiProton();

   TGraphErrors*  gr1 = new TGraphErrors(NPointsPbar_PionMult,MultX,MultValue,0,MultError);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   gr1->GetXaxis()->SetTitle("Charged Pion Momentum (GeV/c)");
   gr1->GetYaxis()->SetTitle("MC/Data (dN/dP (GeV/c)^{-1} / NEvents)");
   
   return gr1;

}

double calcChi2PionMult( std::string target="H", std::string model="FTF", int NDF )
{
   double chi2 = 0.;
//   int NDF = 0;

   std::string histofile = "antiproton" + target + model;
   histofile += ".root";
   TFile* f = new TFile( histofile.c_str() );
   TH1F* hi_mult = (TH1F*)f->Get("NPions");
   
   TGraphErrors* gdata = getPionMultAsGraph();
   
   chi2 = Chi2( gdata, hi_mult, NDF );
   
   std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << (chi2/NDF) << std::endl;
   
   return chi2;
   
}

TGraphErrors* getChPiMomAsGraph()
{

   readAntiProton();

   TGraphErrors*  gr1 = new TGraphErrors(NPointsPbar_PionMom,MomX,MomValue,0,MomError);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   gr1->GetXaxis()->SetTitle("Charged Pion Momentum (GeV/c)");
   gr1->GetYaxis()->SetTitle("MC/Data (dN/dP (GeV/c)^{-1} / NEvents)");
   
   return gr1;

}

double calcChi2ChPiMom( std::string target="H", std::string model="FTF", int& NDF )
{
   double chi2 = 0.;
//   int NDF = 0;

   std::string histofile = "antiproton" + target + model;
   histofile += ".root";
   TFile* f = new TFile( histofile.c_str() );
   TH1F* hi_mom = (TH1F*)f->Get("ChargedPionMomentum");
   
   TGraphErrors* gdata = getChPiMomAsGraph();
   
   chi2 = Chi2( gdata, hi_mom, NDF );
   
   std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << (chi2/NDF) << std::endl;
   
   return chi2;
   
}

void drawAntiProtonMC2DataMom( std::string target="H" )
{

   readAntiProton();
   
   float YY[NPointsPbar_PionMom], DYY[NPointsPbar_PionMom];

   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;


   for ( int ip=0; ip<NPointsPbar_PionMom; ip++)
   {
            
      YY[ip] = 1.0;
      DYY[ip] = MomError[ip] / MomValue[ip] ;
      if ( (YY[ip]+DYY[ip]) > ymax ) ymax = YY[ip]+DYY[ip];
      if ( (YY[ip]-DYY[ip]) < ymin ) ymin = YY[ip]-DYY[ip];
   }

   TGraph*  gr1 = new TGraphErrors(NPointsPbar_PionMom,MomX,YY,0,DYY);

   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   gr1->GetXaxis()->SetTitle("Charged Pion Momentum (GeV/c)");
   gr1->GetYaxis()->SetTitle("MC/Data (dN/dP (GeV/c)^{-1} / NEvents)");
  
   float MC2DataX[NPointsPbar_PionMom];
   float MC2DataY[NPointsPbar_PionMom];
   float DX[NPointsPbar_PionMom], DY[NPointsPbar_PionMom];
   int np =0;
   
   TGraph* gr[NModelsBaryons-1];
   
   TH1F* hi_mult[NModelsBaryons-1];

   for ( int m=0; m<NModelsBaryons; m++ )
   {
//      std::string histofile = "antiproton" + target + ModelsBaryons[m+1]; // shift by 1 to skip "stopping"
      std::string histofile = "antiproton" + target + ModelsBaryons[m]; // "stopping" has been decommisioned
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      hi_mult[m] = (TH1F*)f->Get("ChargedPionMomentum");
      // turn off the stats pad; use this space to draw color codes, etc.
      int nx = hi_mult[m]->GetNbinsX();
      
      // std::cout << " nx = " << nx << std::endl;
      
      // yes... but nx != NPointsMadey !!!
      
      np=0;
      for (int k=1; k <= nx; k++) {
        double xx1 = hi_mult[m]->GetBinLowEdge(k);
	double xx2 = hi_mult[m]->GetBinWidth(k);
	// std::cout << " xx1= " << xx1 << " xx2= " << xx2 << std::endl;
	for (int kk=0; kk<NPointsPbar_PionMom; kk++ )
	{
	   if ( xx1 < MomX[kk] && xx1+xx2 > MomX[kk] )
	   {
	      double yy = hi_mult[m]->GetBinContent(k);
	      MC2DataX[np] = MomX[kk];
	      DX[np] = 0.;
	      MC2DataY[np] = yy / MomValue[kk];
	      // also error calc here !...
	      DY[np]=0.;
	      if ( (MC2DataY[np]+DY[np]) > ymax ) ymax = MC2DataY[np]+DY[np];
	      if ( (MC2DataY[np]-DY[np]) < ymin ) ymin = MC2DataY[np]-DY[np];
	      np++;
	      break;
	   }
	}
      }
      gr[m] = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
      gr[m]->SetTitle(hi_mult[m]->GetTitle());
      gr[m]->GetXaxis()->SetTitle("Charged Pion Momentum (GeV/c)");
      gr[m]->GetYaxis()->SetTitle("MC/Data (dN/dP (GeV/c)^{-1} / NEvents)");
      gr[m]->SetMarkerColor(ColorModel[m+1]);  
      gr[m]->SetMarkerStyle(SymbModel[m+1]);
      gr[m]->SetMarkerSize(1.6);
      
      if ( m==0 ) gr1->SetTitle(hi_mult[m]->GetTitle());
      
   }
   
   gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->Draw("apl");    
   
   for ( int m=0; m<NModelsBaryons-1; m++ )
   {
      gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[m]->Draw("lpsame");
   }
  
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9); 
     
   for ( int m=0; m<NModelsBaryons-1; m++ )
   {
      leg->AddEntry( gr[m], ModelsBaryons[m+1].c_str(), "p" );
   }

   leg->AddEntry( gr1, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);

   return;

} 

void plotAntiProton( std::string target="H" )
{

   readAntiProton();
   
   TH1F* hi_mom[NModelsBaryons];
   TH1F* hi_mult[NModelsBaryons];
   
   TCanvas *myc = new TCanvas("myc","pbar annihilation on H",1100,500);

// This part is for plotting model-vs data only.
// It needs to be commented out if both model-vs-data and regression are wanted (see below).
// However, leave the TCanvas instanciation here !

/*
   myc->Divide(2,1);
   myc->cd(1); gPad->SetLeftMargin(0.15);
   //gPad->SetLogx(1); gPad->SetLogy(1);
   myc->cd(2);gPad->SetLeftMargin(0.15); gPad->SetLogy();
   //gPad->SetLogx(1); gPad->SetLogy(1);
*/   
   // This part is for plotting both model-vs-data & regression, all in one canvas. 
   // Make sure the TCanvas is instantiated !
   //

/*
   TPad* pad1 = new TPad("pad1","pbar annihilation on H",0.01, 0.01,0.99,0.99);
   pad1->Draw();
   pad1->Divide(1,2);
   pad1->cd(1); gPad->SetLeftMargin(0.15);
   pad1->cd(2); gPad->SetLeftMargin(0.15);
*/
   TText* txt = new TText( 0.35, 0.93, "pbar annihilation on H" );
   txt->SetTextSize(0.05);
   txt->Draw();
   
   TPad* pad1 = new TPad("pad1","",0.01, 0.01,0.49,0.92);
   pad1->Draw();
   TPad* pad2 = new TPad("pad2","",0.51, 0.01,0.99,0.92);
   pad2->Draw();
   // pad2->SetLogy();
      
   double ymin_mom = 100000., ymin_mult = 100000. ; // something big... don't know if I can use FLT_MAX
   double ymax_mom = -1., ymax_mult = -1. ;
   
   // now run over exp.data to determine min/max
   //
   for ( int ipp=0; ipp<NPointsPbar_PionMom; ipp++ )
   {
      if ( MomValue[ipp]+MomError[ipp] > ymax_mom ) ymax_mom = MomValue[ipp]+MomError[ipp];
      if ( MomValue[ipp]-MomError[ipp] < ymin_mom )
      {
         ymin_mom = MomValue[ipp]-MomError[ipp];
	 if ( ymin_mom < 0. ) ymin_mom = 0.;
      }
   }
   
   
   for ( int m=0; m<NModelsBaryons; m++ )
   {
      std::string histofile = "antiproton" + target + ModelsBaryons[m];
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      // f->ls();
      hi_mom[m] = (TH1F*)f->Get("ChargedPionMomentum");
      hi_mom[m]->SetTitle("");
      hi_mult[m] = (TH1F*)f->Get("NPions");
      hi_mult[m]->SetTitle("");
      hi_mom[m]->SetLineColor(ColorModel[m]);
      hi_mult[m]->SetLineColor(ColorModel[m]);
      hi_mom[m]->SetLineWidth(2);
      hi_mult[m]->SetLineWidth(2);
      hi_mom[m]->GetXaxis()->SetTitle("Charged Pion Momentum (GeV/c)");
      hi_mom[m]->GetYaxis()->SetTitle("dN/dP (GeV/c)^{-1} / NEvents");
      hi_mult[m]->GetXaxis()->SetTitle("Number of pions (#pi^{+} + #pi^{-} + #pi^{0})");
//      hi_mult[m]->GetYaxis()->SetTitle("Fraction [%] of events" );
      hi_mult[m]->GetYaxis()->SetTitle("Fraction of events" );
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
//       myc->cd(1);
      pad1->cd();
      if ( m == 0 ) hi_mom[m]->Draw("histo");
      else hi_mom[m]->Draw("histosame");
//      myc->cd(2);
      pad2->cd();
      if ( m == 0 ) hi_mult[m]->Draw("histo");
      else hi_mult[m]->Draw("histosame");      
   }
   
   TLegend* leg1 = new TLegend(0.6, 0.70, 0.9, 0.9);
   TLegend* leg2 = new TLegend(0.6, 0.70, 0.9, 0.9);

   for ( int m=0; m<NModelsBaryons; m++ )
   {
      //hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*10.);
      hi_mom[m]->GetYaxis()->SetRangeUser(ymin_mom,ymax_mom*1.4);
      hi_mom[m]->SetStats(0);
      leg1->AddEntry( hi_mom[m], ModelsBaryons[m].c_str(), "L" );
      // hi_mult[m]->GetYaxis()->SetRangeUser(ymin_mult,ymax_mult*1.1);
      hi_mult[m]->GetYaxis()->SetRangeUser(ymin_mult,1.);
      hi_mult[m]->SetStats(0);
      leg2->AddEntry( hi_mult[m], ModelsBaryons[m].c_str(), "L" );
   }
   
   // readAntiProton();
   
   TGraph*  gr1 = new TGraphErrors(NPointsPbar_PionMom,MomX,MomValue,0,MomError);
   TGraph*  gr2 = new TGraphErrors(NPointsPbar_PionMult,MultX,MultValue,0,MultError);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   
   gr2->SetMarkerColor(4);  gr2->SetMarkerStyle(22);
   gr2->SetMarkerSize(1.6);
   
//    myc->cd(1);
   pad1->cd();
   gr1->Draw("p"); 
   leg1->AddEntry(gr1, "exp.data", "p");
   leg1->Draw();
   leg1->SetFillColor(kWhite);
     
//   myc->cd(2);
   pad2->cd();
   gr2->Draw("p");
   leg2->AddEntry(gr2, "exp.data", "p");
   leg2->Draw();
   leg2->SetFillColor(kWhite);
   
   myc->cd();
   myc->Print( "pbar-H-models.gif");

   return;

}


void drawAntiProtonRegression( std::string target="H", std::string model="FTF" )
{

   TH1F* hi_mom[NVersions];
   TH1F* hi_mult[NVersions];
   
// This part is for plotting model-vs data only.
// It needs to be commented out if both model-vs-data and regression are wanted (see below).
// TCanvas isn't needed here, because it gets created in the earlier method ! 
//
   TCanvas* myc = new TCanvas("myc", "", 800, 600);
   myc->Divide(2,1);
   myc->cd(1); gPad->SetLeftMargin(0.15);
   myc->cd(2); gPad->SetLeftMargin(0.15);

   // this part is for plotting both models-vs-data & regression, all in one canvas
   //
//   TPad* pad2 = new TPad("pad2","",0.51, 0.01,0.99,0.99);
//   pad2->Draw();
//   pad2->Divide(1,2);
//   pad2->cd(1); gPad->SetLeftMargin(0.15);
//   pad2->cd(2); gPad->SetLeftMargin(0.15);

   double ymin_mom = 100000., ymin_mult = 100000. ; // something big... don't know if I can use FLT_MAX
   double ymax_mom = -1., ymax_mult = -1. ;

   readAntiProton();

   for ( int ip=0; ip<NPointsPbar_PionMom; ip++)
   {
            
      if ( (MomValue[ip]+MomError[ip]) > ymax_mom ) ymax_mom = MomValue[ip]+MomError[ip];
      if ( (MomValue[ip]-MomError[ip]) < ymin_mom ) ymin_mom = MomValue[ip]-MomError[ip];
   }
   for ( int ip=0; ip<NPointsPbar_PionMult; ip++)
   {
            
      if ( (MultValue[ip]+MultError[ip]) > ymax_mult ) ymax_mult = MultValue[ip]+MultError[ip];
      if ( (MultValue[ip]-MultError[ip]) < ymin_mult ) ymin_mult = MultValue[ip]-MultError[ip];
   }
   
   for ( int iv=0; iv<NVersions; iv++ )
   {
//      std::string histofile = Versions[iv] + "/" + "antiproton" + target + model;
//      histofile += ".root";

      std::string location = "";
      if ( Versions[iv] == CurrentVersion || Versions[iv] == "." )
      {
         location = "";
      }
      else
      {
         location = regre_test_dir + "/test48/" + Versions[iv] + "/";
      }
      std::string histofile = location + "antiproton" + target + model + ".root"; 

      TFile* f = new TFile( histofile.c_str() );
      // f->ls();
      hi_mom[iv] = (TH1F*)f->Get("ChargedPionMomentum");
      hi_mult[iv] = (TH1F*)f->Get("NPions");
      hi_mom[iv]->SetLineColor(ColorVersion[iv]);
      hi_mult[iv]->SetLineColor(ColorVersion[iv]);
      hi_mom[iv]->SetLineWidth(6-iv);
//      hi_mom[iv]->SetLineWidth(2);
      hi_mult[iv]->SetLineWidth(6-iv);
//      hi_mult[iv]->SetLineWidth(2);
      hi_mom[iv]->GetXaxis()->SetTitle("Charged Pion Momentum (GeV/c)");
      hi_mom[iv]->GetXaxis()->SetTitleSize(0.04);
      hi_mom[iv]->GetXaxis()->CenterTitle();
      hi_mom[iv]->GetYaxis()->SetTitle("dN/dP (GeV/c)^{-1} / NEvents");
      hi_mult[iv]->GetXaxis()->SetTitle("Number of pions (#pi^{+} + #pi^{-} + #pi^{0})");
      hi_mult[iv]->GetXaxis()->SetTitleSize(0.04);
      hi_mult[iv]->GetXaxis()->CenterTitle();
//      hi_mult[iv]->GetYaxis()->SetTitle("Fraction [%] of events" );
      hi_mult[iv]->GetYaxis()->SetTitle("Fraction of events" );
      hi_mom[iv]->GetYaxis()->SetTitleOffset(1.5);
      hi_mom[iv]->GetYaxis()->SetTitleSize(0.04);
      hi_mult[iv]->GetYaxis()->SetTitleOffset(1.5);
      hi_mult[iv]->GetYaxis()->SetTitleSize(0.04);
      int nx = hi_mom[iv]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi_mom[iv]->GetBinContent(k);
	if ( yy > ymax_mom ) ymax_mom = yy;
	if ( yy < ymin_mom && yy > 0. ) ymin_mom = yy;
      }
      nx = hi_mult[iv]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi_mult[iv]->GetBinContent(k);
	if ( yy > ymax_mult ) ymax_mult = yy;
	if ( yy < ymin_mult && yy > 0. ) ymin_mult = yy;
      }      
      myc->cd(1);
//      pad2->cd(1);
      if ( iv == 0 ) hi_mom[iv]->Draw();
      else hi_mom[iv]->Draw("same");
      myc->cd(2);
//      pad2->cd(2);
      if ( iv == 0 ) hi_mult[iv]->Draw();
      else hi_mult[iv]->Draw("same");      
   
   }

   TLegend* leg1 = new TLegend(0.45, 0.70, 0.95, 0.9);
   leg1->SetTextSize(0.03);
   TLegend* leg2 = new TLegend(0.45, 0.70, 0.95, 0.9);
   leg2->SetTextSize(0.03);
   
   for ( int iv=0; iv<NVersions; iv++ )
   {
      hi_mom[iv]->GetYaxis()->SetRangeUser(ymin_mom,ymax_mom*1.4); // hi[m]->SetTitle("");
      hi_mom[iv]->SetStats(0);
      std::string entry1 = model + ", " + Versions[iv];
      leg1->AddEntry( hi_mom[iv], entry1.c_str(), "L" );
      hi_mult[iv]->GetYaxis()->SetRangeUser(ymin_mult,1.); // hi[m]->SetTitle("");
      hi_mult[iv]->SetStats(0);
      std::string entry2 = model + ", " + Versions[iv];
      leg2->AddEntry( hi_mult[iv], entry2.c_str(), "L" );
   }
   
   TGraph*  gr1 = new TGraphErrors(NPointsPbar_PionMom,MomX,MomValue,0,MomError);
   TGraph*  gr2 = new TGraphErrors(NPointsPbar_PionMult,MultX,MultValue,0,MultError);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   
   gr2->SetMarkerColor(4);  gr2->SetMarkerStyle(22);
   gr2->SetMarkerSize(1.6);
   gr2->GetYaxis()->SetRangeUser(0., 1.);
   
   myc->cd(1);
//   pad2->cd(1); 
   gr1->Draw("p"); 
   leg1->AddEntry(gr1, "exp.data", "p");
   leg1->Draw();
   leg1->SetFillColor(kWhite);
   myc->cd(2);
//   pad2->cd(2); 
   gr2->Draw("p"); 
   leg2->AddEntry(gr2, "exp.data", "p");
   leg2->Draw();
   leg2->SetFillColor(kWhite);
   
   myc->cd();
   myc->Print("pbar-H-FTF-regre.gif");
     
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

#endif
