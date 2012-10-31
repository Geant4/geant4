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


const int NModelsMuons=3;
std::string ModelsMuons[3] = { "stopping", "captureUpdate", "CHIPS" };

// int         ColorModel[4]    = { 1, 6, 2, 3 };
int         ColorModel[4]    = { 1, 3, 2, 6 };
int         SymbModel[4]     = { 20, 29, 21, 8 };

// pi- beam business
const int   NPointsNeutMultSinger = 5;
const int   NTargetsSingerExp = 8; 
const int   NTargetsSingerTheo = 4;
float NeutMult[5];  
float NeutMultTheo[5];
float ValueExp[8][5], ErrorExp[8][5];
float ValueTheo[8][5];
std::string TargetsSingerExp[8] = { "Al", "Si", "Ca", "Fe", "Ag", "I", "Au", "Pb" };
std::string TargetsSingerTheo[4] = { "Ag", "I", "Au", "Pb" };

const int NVersions = 3;
int ColorVersion[4] = { kBlack, kRed, kGreen, kMagenta };

std::string Versions[3] = { "geant4-09-04-ref10", "geant4-09-05", "geant4-09-05-ref01" };
// std::string Versions[2] = { "geant4-09-04-ref10", "geant4-09-05-ref01" };

 
void readNeutMultSinger()
{
 
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
      
   return;

}

void plotMuMinus( std::string target )
{

   readNeutMultSinger();

   TCanvas *myc = new TCanvas("myc","",800,600);
   //myc->SetLogx(1);
   //myc->SetLogy(1);
   
   drawMuMinus( target );
   // drawPiMinusMC2Data( target );
   
   return;   
}

void drawMuMinus( std::string target )
{
      
   int TargetID = findTargetExp( target );
         
   if ( TargetID == -1 || TargetID >= NTargetsSingerExp )
   {
      std::cout << " Invalid Target: " << target << std::endl;
      return;
   }
   
   int TargetIDTheo = findTargetTheo( target ); // if exists...
   
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
      hi[m]->GetXaxis()->SetTitle("Number of secondary neutrons per mu- capture");
      hi[m]->GetYaxis()->SetTitle("Normalized yield");
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
   gr1->Draw("p");
   
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

/*
void plotPiMinusDataAndRegression( std::string target )
{

   TCanvas *myc = new TCanvas("myc","",800,600);
   myc->Divide(2,1);
   
   readMadey();
   
   myc->cd(1);
   gPad->SetLogx();
   gPad->SetLogy();
   drawPiMinus( target );
   
   myc->cd(2);
   gPad->SetLogx();
   gPad->SetLogy();
   drawPiMinusRegression( target, "BertiniPreCo" );   
   
   return;

}

void plotPiMinusAll()
{

   readMadey();

   TCanvas *myc = new TCanvas("myc","",1000,600);
   myc->Divide(4,2);

   for ( int i=0; i<NTargetsMadey; i++ )
   {      
      myc->cd(i+1);  gPad->SetLogy(1);
      drawPiMinus( TargetsMadey[i] ); 
   }

   return ;
}


void plotPiMinusRegressionSummary2( std::string target, std::string model )
{

   
   readMadey();
   
   TCanvas* myc = new TCanvas("myc","",800,600);
   
   myc->Divide(2,1);
   
   myc->cd(1);
   gPad->SetLogx();
   gPad->SetLogy();
   drawPiMinusRegression( target, model );
   
   myc->cd(2);
   gPad->SetLogx();
   drawPiMinusMC2DataRegression( target, model );
   
   myc->cd();
   
   return;

}

void plotPiMinusSummary2( std::string target )
{

   readMadey();

   TCanvas *myc = new TCanvas("myc","",800,600);
   
   myc->Divide(2,1);
   
   myc->cd(1);
   gPad->SetLogx(1);
   gPad->SetLogy(1);
   drawPiMinus( target );

   myc->cd(2);
   gPad->SetLogx(1);
   drawPiMinusMC2Data( target );
   
   myc->cd();
   
   return;   

}

void plotPiMinusSummary4( std::string target )
{

   readMadey();
   
   TCanvas* myc = new TCanvas("myc", "", 1200, 700 );
   
   myc->Divide(2,2);
   
   myc->cd(1); gPad->SetLogx(); 
   drawPiMinus( target );
   
   myc->cd(2); gPad->SetLogy();
   drawPiMinus( target );
   
   myc->cd(3); gPad->SetLogx(); gPad->SetLogy();
   drawPiMinus( target );
   
   myc->cd(4); myc->SetLogx();
   drawPiMinusMC2Data( target );
   
   myc->cd();

   return;

}

void drawPiMinusMC2Data( std::string target )
{
      
   int TargetID = findTarget( target );
         
   if ( TargetID == -1 || TargetID >= NTargetsMadey )
   {
      std::cout << " Invalid Target: " << target << std::endl;
      return;
   }
   
   float YY[NPointsMadey], DYY[NPointsMadey];
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;

   for ( int ip=0; ip<NPointsMadey; ip++)
   {
            
      YY[ip] = 1.0;
      DYY[ip] = Error[TargetID][ip] / Value[TargetID][ip] ;
      if ( (YY[ip]+DYY[ip]) > ymax ) ymax = YY[ip]+DYY[ip];
      if ( (YY[ip]-DYY[ip]) < ymin ) ymin = YY[ip]-DYY[ip];
   }

   // TGraph*  gr1 = new TGraphErrors(NPointsMadey,KENeut,Value[TargetID],0,Error[TargetID]);
   TGraph*  gr1 = new TGraphErrors(NPointsMadey,KENeut,YY,0,DYY);
   // gr1->SetTitle("");
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   gr1->GetXaxis()->SetTitle("Kinetic energy of secondary neutron (MeV)");
   gr1->GetYaxis()->SetTitle("MC/Data (Number of neutrons per MeV)");


   TH1F* hi[NModelsMesons];
         
   float MC2DataX[NPointsMadey];
   float MC2DataY[NPointsMadey];
   float DX[NPointsMadey], DY[NPointsMadey];
   int np =0;
   
   TGraph* gr[NModelsMesons];
   
   for ( int m=0; m<NModelsMesons; m++ )
   {
      std::string histofile = "piminus" + TargetsMadey[TargetID] + ModelsMesons[m];
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      hi[m] = (TH1F*)f->Get("NvsT");
      // turn off the stats pad; use this space to draw color codes, etc.
      int nx = hi[m]->GetNbinsX();
      
      // std::cout << " nx = " << nx << std::endl;
      
      // yes... but nx != NPointsMadey !!!
      
      np=0;
      for (int k=1; k <= nx; k++) {
        double xx1 = hi[m]->GetBinLowEdge(k);
	double xx2 = hi[m]->GetBinWidth(k);
	// std::cout << " xx1= " << xx1 << " xx2= " << xx2 << std::endl;
	for (int kk=0; kk<NPointsMadey; kk++ )
	{
	   if ( xx1 < KENeut[kk] && xx1+xx2 > KENeut[kk] )
	   {
	      double yy = hi[m]->GetBinContent(k);
	      MC2DataX[np] = KENeut[kk];
	      DX[np] = 0.;
	      MC2DataY[np] = yy / Value[TargetID][kk];
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
      gr[m]->SetTitle(hi[m]->GetTitle());
      gr[m]->GetXaxis()->SetTitle("Kinetic energy of secondary neutron (MeV)");
      gr[m]->GetYaxis()->SetTitle("MC/Data (Number of neutrons per MeV)");
      gr[m]->SetMarkerColor(ColorModel[m]);  
      gr[m]->SetMarkerStyle(SymbModel[m]);
      gr[m]->SetMarkerSize(1.6);
      
      if ( m==0 ) gr1->SetTitle(hi[m]->GetTitle());
      
   }

   gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->Draw("apl");    
   
   for ( int m=0; m<NModelsMesons; m++ )
   {
      gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[m]->Draw("lpsame");
   }
  
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);   
   
   for ( int m=0; m<NModelsMesons; m++ )
   {
      leg->AddEntry( gr[m], ModelsMesons[m].c_str(), "p" );
   }

   leg->AddEntry( gr1, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);

   return;

}

void drawPiMinusMC2DataRegression( std::string target, std::string model )
{
   
   // readMadey();
   
   int TargetID = findTarget( target );
         
   if ( TargetID == -1 || TargetID >= NTargetsMadey )
   {
      std::cout << " Invalid Target: " << target << std::endl;
      return;
   }
   
   float YY[NPointsMadey], DYY[NPointsMadey];
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;

   for ( int ip=0; ip<NPointsMadey; ip++)
   {
            
      YY[ip] = 1.0;
      DYY[ip] = Error[TargetID][ip] / Value[TargetID][ip] ;
      if ( (YY[ip]+DYY[ip]) > ymax ) ymax = YY[ip]+DYY[ip];
      if ( (YY[ip]-DYY[ip]) < ymin ) ymin = YY[ip]-DYY[ip];
   }

   // TGraph*  gr1 = new TGraphErrors(NPointsMadey,KENeut,Value[TargetID],0,Error[TargetID]);
   TGraph*  gr1 = new TGraphErrors(NPointsMadey,KENeut,YY,0,DYY);
   // gr1->SetTitle("");
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   gr1->GetXaxis()->SetTitle("Kinetic energy of secondary neutron (MeV)");
   gr1->GetYaxis()->SetTitle("MC/Data (Number of neutrons per MeV)");


   TH1F* hi[NVersions];
         
   float MC2DataX[NPointsMadey];
   float MC2DataY[NPointsMadey];
   float DX[NPointsMadey], DY[NPointsMadey];
   int np =0;
   
   TGraph* gr[NVersions];

   for ( int iver=0; iver<NVersions; iver++ )
   {
      std::string histofile = Versions[iver] + "/" + "piminus" + target + model;
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      hi[iver] = (TH1F*)f->Get("NvsT");
      // turn off the stats pad; use this space to draw color codes, etc.
      int nx = hi[iver]->GetNbinsX();
      
      // std::cout << " nx = " << nx << std::endl;
      
      // yes... but nx != NPointsMadey !!!
      
      np=0;
      for (int k=1; k <= nx; k++) {
        double xx1 = hi[iver]->GetBinLowEdge(k);
	double xx2 = hi[iver]->GetBinWidth(k);
	// std::cout << " xx1= " << xx1 << " xx2= " << xx2 << std::endl;
	for (int kk=0; kk<NPointsMadey; kk++ )
	{
	   if ( xx1 < KENeut[kk] && xx1+xx2 > KENeut[kk] )
	   {
	      double yy = hi[iver]->GetBinContent(k);
	      MC2DataX[np] = KENeut[kk];
	      DX[np] = 0.;
	      MC2DataY[np] = yy / Value[TargetID][kk];
	      // also error calc here !...
	      DY[np]=0.;
	      if ( (MC2DataY[np]+DY[np]) > ymax ) ymax = MC2DataY[np]+DY[np];
	      if ( (MC2DataY[np]-DY[np]) < ymin ) ymin = MC2DataY[np]-DY[np];
	      np++;
	      break;
	   }
	}
      }
      gr[iver] = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
      gr[iver]->SetTitle(hi[iver]->GetTitle());
      gr[iver]->GetXaxis()->SetTitle("Kinetic energy of secondary neutron (MeV)");
      gr[iver]->GetYaxis()->SetTitle("MC/Data (Number of neutrons per MeV)");
      gr[iver]->SetMarkerColor(ColorVersion[iver]);  
      gr[iver]->SetMarkerStyle(20);
      gr[iver]->SetMarkerSize(1.6);
      
      if ( iver==0 ) 
      {
	 std::string title = "pi- on " + target + ", " + model;
	 gr1->SetTitle( title.c_str() );
	 //gr1->SetTitle( hi[iver]->GetTitle() );
      }
      
   }
   gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->Draw("apl");    
   
   for ( int iver=0; iver<NVersions; iver++ )
   {
      gr[iver]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[iver]->Draw("lpsame");
   }
  
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);   
   
   for ( int iver=0; iver<NVersions; iver++ )
   {
      leg->AddEntry( gr[iver], Versions[iver].c_str(), "p" );
   }

   leg->AddEntry( gr1, "exp.data", "p");
   
   leg->Draw();
   leg->SetFillColor(kWhite);

   return;
}

void plotPiMinusRegression( std::string target="C", std::string model="CHIPS" )
{

   readMadey();
   
   TCanvas *myc = new TCanvas("myc","",800,600);
   myc->SetLogx(1);
   myc->SetLogy(1);
   
   drawPiMinusRegression( target, model );
   
   return;

}

void plotPiMinusRegressionAllTargets( std::string model="CHIPS" )
{

   
   readMadey();
   
   TCanvas *myc = new TCanvas("myc","",1200,900);
   myc->Divide(4,2);

   for ( int i=0; i<NTargetsMadey; i++ )
   {      
      myc->cd(i+1); gPad->SetLogx(1); gPad->SetLogy(1);
      drawPiMinusRegression( TargetsMadey[i], model ); 
   }
   
   myc->cd();

   return;

}


void drawPiMinusRegression( std::string target="C", std::string model="CHIPS" )
{

   int TargetID = findTarget( target );
         
   if ( TargetID == -1 || TargetID >= NTargetsMadey )
   {
      std::cout << " Invalid Target: " << target << std::endl;
      return;
   }

   TH1F* hi[NVersions];
   
   double ymin = 10000.;
   double ymax = -1.;
   for ( int iver=0; iver<NVersions; iver++ )
   {
      std::string histofile = Versions[iver] + "/" + "piminus" + target + model;
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      hi[iver] = (TH1F*)f->Get("NvsT");
      hi[iver]->SetStats(0);
      hi[iver]->SetLineColor(ColorVersion[iver]);
      hi[iver]->SetLineWidth(2);
      hi[iver]->GetXaxis()->SetTitle("Kinetic energy of secondary neutron (MeV)");
      hi[iver]->GetYaxis()->SetTitle("Number of neutrons per MeV");
      int nx = hi[iver]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi[iver]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( iver == 0 ) hi[iver]->Draw();
      else hi[iver]->Draw("same");
      
      
   }
   // std::cout << " ymin= " << ymin << " ymax= " << ymax << std::endl;
   
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
   leg->SetTextSize(0.017);
   
   for ( int iver=0; iver<NVersions; iver++ )
   {
      hi[iver]->GetYaxis()->SetRangeUser(ymin,ymax*10.); // hi[m]->SetTitle("");
      std::string htitle = "pi- on " + target + ", " + model;
      hi[iver]->SetTitle( htitle.c_str() );
      // std::string entry = model + ", " + Versions[iver];
      // leg->AddEntry( hi[iver], entry.c_str(), "L" );
      leg->AddEntry( hi[iver], Versions[iver].c_str(), "L" );
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
*/

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

int findTargetExp(std::string target )
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

int findTargetTheo(std::string target )
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
      
   if ( TargetID == -1 || TargetID >= NTargetsSingerTheo )
   {
      std::cout << " Invalid Target: " << target << std::endl;
      return -1;
   }
   
   return TargetID;

}
