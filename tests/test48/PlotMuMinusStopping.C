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


const int NModelsMuons=2;
std::string ModelsMuons[2] = { "stopping", "captureUpdate" };

// int         ColorModel[4]    = { 1, 6, 2, 3 };
int         ColorModel[4]    = { 1, 3, 2, 6 };
int         SymbModel[4]     = { 20, 29, 21, 8 };

// mu- beam business
const int   NPointsNeutMultSinger = 5;
const int   NTargetsSingerExp = 8; 
const int   NTargetsSingerTheo = 4;
float NeutMult[5];  
float NeutMultTheo[5];
float ValueExp[8][5], ErrorExp[8][5];
float ValueTheo[8][5];
std::string TargetsSingerExp[8] = { "Al", "Si", "Ca", "Fe", "Ag", "I", "Au", "Pb" };
std::string TargetsSingerTheo[4] = { "Ag", "I", "Au", "Pb" };

const int NPointsSundelin = 7;
const int NTargetsSundelin = 3;
std::string TargetsSundelin[3] = { "Si", "S", "Ca" };
float EKinMin[7], EKinMax[7], EKin[7];
float NNeut[3][7], ErNeut[3][7]; 
 

const int NVersions = 3;
int ColorVersion[4] = { kBlack, kRed, kGreen, kMagenta };

std::string Versions[3] = { "geant4-09-04-ref10", "geant4-09-05", "geant4-09-05-ref01" };
// std::string Versions[2] = { "geant4-09-04-ref10", "geant4-09-05-ref01" };

 
void readNeutMultSinger()
{
 
/*
   read data poins extracted from paper by P.Singer
*/   
   
   ifstream infile;

   std::string fname = "./muminus/neutron_mult_singer_exp.dat"; 
   std::cout << "Reads data from file " << fname << "\n";
   infile.open(fname.c_str());
      
   for ( int i=0; i<NPointsNeutMultSinger; i++ )
   {

      infile >> NeutMult[i];
      // shift to the middle of bin
      NeutMult[i] += 0.5;
      for ( int j=0; j<NTargetsSingerExp; j++ )
      {
         infile >> ValueExp[j][i] >> ErrorExp[j][i];
      }
   }
   
   infile.close();
   
   std::string fname1 = "./muminus/neutron_mult_singer_theo.dat";
   std::cout << "Reads data from file " << fname1 << "\n";
   
   infile.open(fname1.c_str());
   
   for ( int i=0; i<NPointsNeutMultSinger; i++ )
   {

      infile >> NeutMultTheo[i];
      // shift to the middle of bin
      NeutMultTheo[i] += 0.5;
      for ( int j=0; j<NTargetsSingerTheo; j++ )
      {
         infile >> ValueTheo[j][i];
      }
   }
   
   infile.close();
      
   return;

}

void readNNeutEKin()
{

/*
   Read data points from R.M.Sundelin et al., 
   "Spectrum of Neutrons from Muon Capture in Silicon, Sulfur, and Calcium"
   Phys.Rev>Lett, Vol.20, Number 21, 1198 (1968)
*/

   int ip=0;
   int it = 0;
   
   for ( ip=0; ip<NPointsSundelin; ip++ )
   {
      EKinMin[ip] = 0;
      EKinMax[ip] = 0;
      EKin[ip] = 0;
      for ( it=0; it<NTargetsSundelin; it++ )
      {
         NNeut[it][ip] = 0;
	 ErNeut[it][ip] = 0;
      }
   }

   ifstream infile;

   std::string fname = "./muminus/num_neut_vs_ekin_exp.dat"; 
   std::cout << "Reads data from file " << fname << "\n";
   infile.open(fname.c_str());
   
   for ( ip=0; ip<NPointsSundelin; ip++ )
   {
      infile >> EKinMin[ip] >> EKinMax[ip];
      EKin[ip] = 0.5 * ( EKinMin[ip] + EKinMin[ip] );
      for ( it=0; it<NTargetsSundelin; it++ )
      {
         infile >> NNeut[it][ip] >> ErNeut[it][ip];
      }
   }
   
   infile.close();
   
   return;

}

void plotMuMinusNeutMult( std::string target )
{

   readNeutMultSinger();
   
   TCanvas* myc = new TCanvas("myc","",800,600);
   //myc->SetLogx(1);
   //myc->SetLogy(1);
   
   drawMuMinusNeutMult( target );
   
   return;   
}

void plotMuMinusNNeutEKin( std::string target )
{

   readNNeutEKin();
   
   TCanvas* myc = new TCanvas("myc","",800,600);
   myc->SetLogy();
   
   drawMuMinusNNeutEKin( target );
     
   return;

}

void plotMuMinNMultForTalk()
{

   readNeutMultSinger();
   
   TCanvas* myc1 = new TCanvas("myc1","",800,800);
   myc1->Divide(2,2);

   myc1->cd(1);
   //gPad->SetLogy();
   drawMuMinNMultForTalk("Al");

   myc1->cd(2);
   // gPad->SetLogy();
   drawMuMinNMultForTalk("Si");
   // drawMuMinusNeutMult("Si");

   myc1->cd(3);
   //gPad->SetLogy();
   drawMuMinNMultForTalk("Ca");
   // drawMuMinusNeutMult("Ca");

   myc1->cd(4);
   //gPad->SetLogy();
   drawMuMinNMultForTalk("Fe");
   // drawMuMinusNeutMult("Fe");

   TCanvas* myc2 = new TCanvas("myc2","",800,800);
   myc2->Divide(2,2);

   myc2->cd(1);
   //gPad->SetLogy();
   drawMuMinNMultForTalk("Ag");
   // drawMuMinusNeutMult("Ag");

   myc2->cd(2);
   //gPad->SetLogy();
   drawMuMinNMultForTalk("I");
   // drawMuMinusNeutMult("I");

   myc2->cd(3);
   //gPad->SetLogy();
   drawMuMinNMultForTalk("Au");
   // drawMuMinusNeutMult("Au");

   myc2->cd(4);
   //gPad->SetLogy();
   drawMuMinNMultForTalk("Pb");
   // drawMuMinusNeutMult("Pb");

   return;

}

void plotMuMinNRateEKinForTalk()
{

   readNNeutEKin();
   
   TCanvas* myc = new TCanvas("myc","",1200,600);
   myc->Divide(3,1);
   
   myc->cd(1);
   gPad->SetLogy();
   drawMuMinNRateEKinForTalk( "Si" );

   myc->cd(2);
   gPad->SetLogy();
   drawMuMinNRateEKinForTalk( "S" );

   myc->cd(3);
   gPad->SetLogy();
   drawMuMinNRateEKinForTalk( "Ca" );
     
   return;

}

void drawMuMinNMultForTalk( std::string target )
{

   int TargetID = findTargetSingerExp( target );
         
   if ( TargetID == -1 || TargetID >= NTargetsSingerExp )
   {
      std::cout << " Invalid Target: " << target << std::endl;
      return;
   }
   
   int TargetIDTheo = findTargetSingerTheo( target ); // if exists...
   
   TH1F* hi[3];

   //double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   //double ymax = -1. ;

   for ( int m=0; m<3; m++ )
   {
      
      std::string histofile ="";
      if ( m == 0 || m == 1 )
      {
         histofile = "9.6.p02/";
      }
      else if ( m == 2 )
      {
         histofile = "9.6.ref08/";
      }
      
      histofile += "muminus" +  TargetsSingerExp[TargetID];
      
      if ( m == 0 )
      {
         histofile += "stopping";
      }
      else
      {
         histofile += "captureUpdate";
      }
      
      histofile += ".root";
      
      TFile* f = new TFile( histofile.c_str() );
      hi[m] = (TH1F*)f->Get("NNeutrons");
      // turn off the stats pad; use this space to draw color codes, etc.
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      if ( m == 1 ) hi[m]->SetLineWidth(5);
      hi[m]->GetXaxis()->Set( 10, 0., 10.);
      hi[m]->GetYaxis()->SetRangeUser( 0., 1. );
      hi[m]->GetXaxis()->SetTitle("Number of secondary neutrons per mu- capture");
      hi[m]->GetYaxis()->SetTitle("Normalized yield");
      hi[m]->GetYaxis()->SetTitleOffset(1.5);
/*
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
*/
      if ( m == 0 ) hi[m]->Draw();
      else hi[m]->Draw("same");
   }
      
   TLegend* leg = new TLegend(0.3, 0.5, 0.9, 0.9);

   leg->AddEntry( "", "CaptureAtRest:", "" );
   leg->AddEntry( hi[0], "9.5.p02 - equivalent", "L" );
   leg->AddEntry( "", "CaptureAtRest - Restructured:", "" );
   leg->AddEntry( hi[1], "9.6.p02", "L" );
   leg->AddEntry( hi[2], "9.6.ref08", "L" );

   TGraph*  gr1 = new TGraphErrors(NPointsNeutMultSinger,NeutMult,ValueExp[TargetID],0,ErrorExp[TargetID]);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.8);
   gr1->Draw("psame");
   
   leg->AddEntry( gr1, "exp.data (P.Singer)", "p");
   
   if ( TargetIDTheo != -1 )
   {
      TGraph*  gr2 = new TGraphErrors(NPointsNeutMultSinger,NeutMultTheo,ValueTheo[TargetIDTheo],0,0);
      gr2->SetMarkerColor(7); // 7 = light blue (turquoise) 
      gr2->SetMarkerStyle(21);
      gr2->SetMarkerSize(1.8);
      gr2->Draw("psame");     
      leg->AddEntry( gr2, "theory (P.Singer)", "p");
 
   }
   
   leg->Draw();
   leg->SetFillColor(kWhite);

   return;

} 



void drawMuMinNRateEKinForTalk( std::string target )
{

   int TargetID = findTargetSundelin( target );
   
   if ( TargetID == -1 ||  TargetID >= NTargetsSundelin ) return;
      
   TH1F* hi[3];

   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;

   for ( int m=0; m<3; m++ )
   {
      
      std::string histofile ="";
      if ( m == 0 || m == 1 )
      {
         histofile = "9.6.p02/";
      }
      else if ( m == 2 )
      {
         histofile = "9.6.ref08/";
      }
      
      histofile += "muminus" +  TargetsSundelin[TargetID];
      
      if ( m == 0 )
      {
         histofile += "stopping";
      }
      else
      {
         histofile += "captureUpdate";
      }
      
      histofile += ".root";
      
      TFile* f = new TFile( histofile.c_str() );

      hi[m] = (TH1F*)f->Get("NeutronKineticEnergy");
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      if ( m == 1 ) hi[m]->SetLineWidth(5);
      // -> hi[m]->GetXaxis()->Set( 10, 0., 10.);
      hi[m]->GetYaxis()->SetRangeUser( 0.0000001, 1. );
      hi[m]->GetXaxis()->SetTitle("Ekin of secondary neutron (MeV)");
      hi[m]->GetYaxis()->SetTitle("Number of neutrons per capture per MeV");

      hi[m]->GetYaxis()->SetTitleOffset(1.5);
/*      
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
*/
      if ( m == 0 ) hi[m]->Draw();
      else hi[m]->Draw("same");
   }
      
   TLegend* leg = new TLegend(0.4, 0.6, 0.9, 0.9);
   leg->SetTextSize(0.035);

   leg->AddEntry( "", "CaptureAtRest", "" );
   leg->AddEntry( hi[0], "9.5.p02 - equivalent", "L" );
   leg->AddEntry( "", "CaptureAtRest-Restruct.", "" );
   leg->AddEntry( hi[1], "9.6.p02", "L" );
   leg->AddEntry( hi[2], "9.6.ref08", "L" );

   TGraph*  gr1 = new TGraphErrors(NPointsSundelin,EKin,NNeut[TargetID],0,ErNeut[TargetID]);
   // TGraph*  gr1 = new TGraphErrors(NPointsNeutMultSinger,NeutMult,ValueExp[TargetID],0,ErrorExp[TargetID]);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.8);
   gr1->Draw("psame");
   
   leg->AddEntry( gr1, "exp.data (R.M.Sundelin)", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);

   return;

} 

void drawMuMinusNeutMult( std::string target )
{
      
   int TargetID = findTargetSingerExp( target );
         
   if ( TargetID == -1 || TargetID >= NTargetsSingerExp )
   {
      std::cout << " Invalid Target: " << target << std::endl;
      return;
   }
   
   int TargetIDTheo = findTargetSingerTheo( target ); // if exists...
   
   TH1F* hi[NModelsMuons];
      
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int m=0; m<NModelsMuons; m++ )
   {
      std::string histofile = "muminus" + TargetsSingerExp[TargetID] + ModelsMuons[m];
      histofile += ".root";

      std::cout << "About to open file: " << histofile << std::endl;

      TFile* f = new TFile( histofile.c_str() );
      hi[m] = (TH1F*)f->Get("NNeutrons");
      // turn off the stats pad; use this space to draw color codes, etc.
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      hi[m]->GetXaxis()->SetLimits(0.,10.);
      hi[m]->GetXaxis()->SetTitle("Number of secondary neutrons per mu- capture");
      hi[m]->GetYaxis()->SetTitle("Normalized yield");
      hi[m]->GetYaxis()->SetTitleOffset(1.5);
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

   for ( int m=0; m<NModelsMuons; m++ )
   {
      // hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*5.); // hi[m]->SetTitle("");
      double delta = 0.5;
      if ( delta > ymax/2. ) delta = ymax/2.;
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax+delta); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelsMuons[m].c_str(), "L" );
   }
      
      // draw a frame to define the range
//   TH1F *hr = myc->DrawFrame(-1.,-0.1,10.,1.);
//   hr->SetXTitle("X title");
//   hr->SetYTitle("Y title");
//   myc->GetFrame()->SetFillColor(21);
//   myc->GetFrame()->SetBorderSize(12);
   
   TGraph*  gr1 = new TGraphErrors(NPointsNeutMultSinger,NeutMult,ValueExp[TargetID],0,ErrorExp[TargetID]);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.8);
   gr1->Draw("psame");
   
   leg->AddEntry( gr1, "exp.data", "p");
   
   if ( TargetIDTheo != -1 )
   {
      TGraph*  gr2 = new TGraphErrors(NPointsNeutMultSinger,NeutMultTheo,ValueTheo[TargetIDTheo],0,0);
      gr2->SetMarkerColor(7); // 7 = light blue (turquoise) 
      gr2->SetMarkerStyle(21);
      gr2->SetMarkerSize(1.8);
      gr2->Draw("psame");     
      leg->AddEntry( gr2, "theory", "p");
 
   }

   leg->Draw();
   leg->SetFillColor(kWhite);

   return;

}

void drawMuMinusNNeutEKin( std::string target )
{

   int TargetID = findTargetSundelin( target );
   
   if ( TargetID == -1 ||  TargetID >= NTargetsSundelin ) return;
   
   TH1F* hi[NModelsMuons] ;

   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int m=0; m<NModelsMuons; m++ )
   {
      std::string histofile = "muminus" + TargetsSundelin[TargetID] + ModelsMuons[m] ;
      histofile += ".root";
      
      std::cout << "About to open file: " << histofile << std::endl;
      
      TFile* f = new TFile( histofile.c_str() );
      hi[m] = (TH1F*)f->Get("NeutronKineticEnergy");
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      hi[m]->GetXaxis()->SetTitle("Ekin of secondary neutron (MeV)");
      hi[m]->GetYaxis()->SetTitle("Number of neutrons per capture per MeV");
      hi[m]->GetYaxis()->SetTitleOffset(1.5);
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) hi[m]->Draw();
      else hi[m]->Draw("same");   
   }   

   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);

   for ( int m=0; m<NModelsMuons; m++ )
   {
      // hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*5.); // hi[m]->SetTitle("");
      double delta = 0.5;
      if ( delta > ymax/2. ) delta = ymax/2.;
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax+delta); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelsMuons[m].c_str(), "L" );
   }

   TGraph*  gr1 = new TGraphErrors(NPointsSundelin,EKin,NNeut[TargetID],0,ErNeut[TargetID]);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.8);
   gr1->Draw("psame");
   
   leg->AddEntry( gr1, "exp.data", "p");
   
   leg->Draw();
   leg->SetFillColor(kWhite);
   
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

int findTargetSingerExp( std::string target )
{

   int TargetID = -1;
   
   for ( int i=0; i<NTargetsSingerExp; i++ )
   {
      if ( TargetsSingerExp[i] == target ) 
      {
         TargetID = i;
	 break;
      }
   }
      
   if ( TargetID == -1 || TargetID >= NTargetsSingerExp )
   {
      std::cout << " Invalid Target: " << target << std::endl;
      return -1;
   }
   
   return TargetID;

}

int findTargetSingerTheo( std::string target )
{

   int TargetID = -1;
   
   for ( int i=0; i<NTargetsSingerTheo; i++ )
   {
      if ( TargetsSingerTheo[i] == target ) 
      {
         TargetID = i;
	 break;
      }
   }
      
   
   return TargetID;

}

int findTargetSundelin( std::string target )
{
   
   int TargetID = -1;
   
   for ( int i=0; i<NTargetsSundelin; i++ )
   {
      if ( TargetsSundelin[i] == target )
      {
         TargetID = i;
	 break;
      }
   
   }

   if ( TargetID == -1 || TargetID >= NTargetsSundelin )
   {
      std::cout << " Invalid Target: " << target << std::endl;
      return -1;
   }
     
   return TargetID;

}
