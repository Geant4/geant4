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

// provision for near-future use of FTF for baryons
//
const int NModelsBaryons=3;
std::string ModelsBaryons[3] = { "CHIPS", "stopping", "FTF" };

const int NModelsMesons=4;
std::string ModelsMesons[4] = { "CHIPS", "stopping", "Bertini", "BertiniPreCo" };

int         ColorModel[4]    = { 1, 2, 6, 3 };
int         SymbModel[4]     = { 20, 29, 21, 8 };

// pi- beam business
const int   NPointsMadey = 30;
const int   NTargetsMadey = 7; 
float KENeut[30];  
float Value[7][30], Error[7][30];
std::string TargetsMadey[7] = { "C", "N", "O", "Al", "Cu", "Ta", "Pb" };

// pbar beam business
const int NPointsPbar_PionMom  = 44;
const int NPointsPbar_PionMult = 6;
float MomX[44], MomValue[44], MomError[44];
float MultX[6], MultValue[6], MultError[6];

// K- beam business
const int NPointsKMinus_Pi0Energy = 82;
float Pi0Energy[82], Pi0EnergyStat[82];


//const int NVersions = 2;
//int ColorVersion[2] = { 1, 2};
//std::string Versions[2] = { "geant4-09-02-patch01", "geant4-09-04-beta01" };
// const int NVersions = 3;
// int ColorVersion[3] = { 1, 2, 3 };
// int ColorVersion[3] = { kBlack, kRed, kGreen };
// std::string Versions[3] = { "geant4-09-02-patch01", "geant4-09-04-beta01", "geant4-09-04-ref05" };
const int NVersions = 2;
int ColorVersion[4] = { kBlack, kRed, kGreen, kMagenta };

//std::string Versions[4] = { "geant4-09-02-patch01", "geant4-09-04-beta01", "geant4-09-04-ref05", "geant4-09-05-beta01" };
// std::string Versions[1] = { "." };
// std::string Versions[3] = { "geant4-09-04-ref10", "geant4-09-05", "geant4-09-05-ref01" };
std::string Versions[2] = { "geant4-09-04-ref10", "geant4-09-05-ref01" };

 
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
      MultX[i] += 0.5;
      MultValue[i] /= 100.; 
      MultError[i] /= 100.;
      // std::cout << MultX[i] << MultValue[i] << MultError[i] << std::endl;
   }
   infile.close();
      
   return;
   
}

void readKMinusPi0Energy()
{

   std::string fname = "./kminus/pi0_energy_larson.dat";
   std::cout << "Reads data from file " << fname << "\n";
   ifstream infile;
   infile.open(fname.c_str());
      
   for ( int i=0; i<NPointsKMinus_Pi0Energy; i++ )
   {
   
      infile >> Pi0Energy[i] >> Pi0EnergyStat[i] ;
      // rescale for (approx.) stat
      Pi0EnergyStat[i] /= 1028571;
      // std::cout << Pi0Energy[i] << " " << Pi0EnergyStat[i] << std::endl;
   
   }
   
   infile.close();
   
   return;

}

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
      myc->cd(i+1);  /* gPad->SetLogx(1); */ gPad->SetLogy(1);
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
   // drawPiMinusMC2Data( target );
   
   return;   
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

void plotKMinus( std::string target )
{
   
   TCanvas* myc1 = new TCanvas("myc1","",1200,600);
   // myc1->Divide(3,2);
   TCanvas* myc2 = new TCanvas("myc2","",1200,600);

   if ( target =="H" || target =="He" )
   {
      //myc1->cd();
      //myc1->Divide( 3, 1 );    
      myc1->Divide(3,2);
      myc1->cd(1); gPad->SetLogy(1);
      drawKMinus( target, "NSecondaries" );
      myc1->cd(2); gPad->SetLogy(1);
      drawKMinus( target, "Topology" );
      myc1->cd(3); gPad->SetLogy(1);
      drawKMinus( target, "PartTypeMult" );  
      //myc2->cd(); 
      //myc2->Divide( 3, 1 );   
      // myc2->cd(1); 
      myc1->cd(4); 
      gPad->SetLogy(1);
      drawKMinus( target, "EnergyPi0" );
      // myc2->cd(2); 
      myc1->cd(5); 
      gPad->SetLogy(1);
      drawKMinus( target, "EnergyPhoton" );
      // myc2->cd(3); 
      myc1->cd(6); 
      gPad->SetLogy(1);
      drawKMinus( target, "EnergyChargedPions" );
   }
   else
   {
      myc1->cd();
      myc1->Divide( 2, 1 );
      myc1->cd(1); gPad->SetLogy(1);
      drawKMinus( target, "NSecondaries" );
      myc1->cd(2); gPad->SetLogy(1);
      drawKMinus( target, "PartTypeMult" );  
      myc2->cd(); 
      myc2->Divide( 2, 1 );   
      myc2->cd(1); gPad->SetLogy(1);
      drawKMinus( target, "EnergyPi0" );
      myc2->cd(2); gPad->SetLogy(1);
      drawKMinus( target, "EnergyChargedPions" );
   }
   
   return;

}

void drawPiMinus( std::string target )
{
      
   int TargetID = findTarget( target );
         
   if ( TargetID == -1 || TargetID >= NTargetsMadey )
   {
      std::cout << " Invalid Target: " << target << std::endl;
      return;
   }
   
   TH1F* hi[NModelsMesons];
      
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int m=0; m<NModelsMesons; m++ )
   {
      std::string histofile = "piminus" + TargetsMadey[TargetID] + ModelsMesons[m];
      histofile += ".root";

      std::cout << "About to open file: " << histofile << std::endl;

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
   
   for ( int m=0; m<NModelsMesons; m++ )
   {
      // hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*5.); // hi[m]->SetTitle("");
      double delta = 0.5;
      if ( delta > ymax/2. ) delta = ymax/2.;
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax+delta); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelsMesons[m].c_str(), "L" );
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

void drawKMinus( std::string target, std::string histo )
{

   TH1F** hi;
      
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   int NCounts;
   if ( target == "H" )
   {
      NCounts = NModelsMesons-1;
      hi = new TH1F*[NModelsMesons-1];
   }
   else
   {
      NCounts = NModelsMesons;
      hi = new TH1F*[NModelsMesons];
   }
   for ( int m=0; m<NCounts; m++ )
   {
      std::string histofile = "kminus" + target + ModelsMesons[m];
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      hi[m] = (TH1F*)f->Get( histo.c_str() );
      // turn off the stats pad; use this space to draw color codes, etc.
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      if ( histo == "EnergyPi0" )
      {
         hi[m]->GetXaxis()->SetTitle("Total energy of secondary pi0 (MeV)" );
         hi[m]->GetYaxis()->SetTitle("dN / dE (MeV)^{-1} / NEvents");
      }
      else
      {
         hi[m]->GetXaxis()->SetTitle( histo.c_str() );
         // hi[m]->GetYaxis()->SetTitle("Number of neutrons per MeV");
      }
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
   
   for ( int m=0; m<NCounts; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelsMesons[m].c_str(), "L" );
   }
   
   if ( histo == "EnergyPi0" && target == "H" )
   {
      readKMinusPi0Energy();
      TGraph* gr = new TGraph( NPointsKMinus_Pi0Energy, Pi0Energy, Pi0EnergyStat );
      gr->SetMarkerColor(4); // blue, alt. would be kBlue
      gr->SetMarkerStyle(22); // solid triangle
      gr->SetMarkerSize(1.6);
      gr->Draw("p");
      leg->AddEntry( gr, "exp.data", "p" );
   }
   
      
   leg->Draw();
   leg->SetFillColor(kWhite);

   return ;

}


void plotAntiProton( std::string target="H" )
{

   TH1F* hi_mom[NModelsBaryons];
   TH1F* hi_mult[NModelsBaryons];
   
   TCanvas *myc = new TCanvas("myc","",800,600);

// This part is for plotting model-vs data only.
// It needs to be commented out if both model-vs-data abd regression are wanted (see below).
// However, leave the TCanvas instanciation here !

   myc->Divide(2,1);
   myc->cd(1); gPad->SetLeftMargin(0.15);
   //gPad->SetLogx(1); gPad->SetLogy(1);
   myc->cd(2);gPad->SetLeftMargin(0.15); 
   //gPad->SetLogx(1); gPad->SetLogy(1);
   
   // This part is for plotting both model-vs-data & regression, all in one canvas. 
   // Make sure the TCanvas is instanciated !
   //
/*
   TPad* pad1 = new TPad("pad1","",0.01, 0.01,0.49,0.99);
   pad1->Draw();
   pad1->Divide(1,2);
   pad1->cd(1); gPad->SetLeftMargin(0.15);
   pad1->cd(2); gPad->SetLeftMargin(0.15);
*/      
   double ymin_mom = 100000., ymin_mult = 100000. ; // something big... don't know if I can use FLT_MAX
   double ymax_mom = -1., ymax_mult = -1. ;
   for ( int m=0; m<NModelsBaryons; m++ )
   {
      std::string histofile = "antiproton" + target + ModelsBaryons[m];
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
//      pad1->cd(1);
      if ( m == 0 ) hi_mom[m]->Draw();
      else hi_mom[m]->Draw("same");
      myc->cd(2);
//      pad1->cd(2);
      if ( m == 0 ) hi_mult[m]->Draw();
      else hi_mult[m]->Draw("same");      
   }
   
   TLegend* leg1 = new TLegend(0.6, 0.70, 0.9, 0.9);
   TLegend* leg2 = new TLegend(0.6, 0.70, 0.9, 0.9);

   for ( int m=0; m<NModelsBaryons; m++ )
   {
      //hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*10.);
      hi_mom[m]->GetYaxis()->SetRangeUser(ymin_mom,ymax_mom*1.4);
      hi_mom[m]->SetStats(0);
      leg1->AddEntry( hi_mom[m], ModelsBaryons[m].c_str(), "L" );
      hi_mult[m]->GetYaxis()->SetRangeUser(ymin_mult,ymax_mult*1.1);
      hi_mult[m]->SetStats(0);
      leg2->AddEntry( hi_mult[m], ModelsBaryons[m].c_str(), "L" );
   }
   
   readAntiProton();
   
   TGraph*  gr1 = new TGraphErrors(NPointsPbar_PionMom,MomX,MomValue,0,MomError);
   TGraph*  gr2 = new TGraphErrors(NPointsPbar_PionMult,MultX,MultValue,0,MultError);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   
   gr2->SetMarkerColor(4);  gr2->SetMarkerStyle(22);
   gr2->SetMarkerSize(1.6);
   
   myc->cd(1);
//   pad1->cd(1);
   gr1->Draw("p"); 
   leg1->AddEntry(gr1, "exp.data", "p");
   leg1->Draw();
   leg1->SetFillColor(kWhite);
     
   myc->cd(2);
//   pad1->cd(2);
   gr2->Draw("p");
   leg2->AddEntry(gr2, "exp.data", "p");
   leg2->Draw();
   leg2->SetFillColor(kWhite);
   
   myc->cd();

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

void drawAntiProtonRegression( std::string target="H", std::string model="CHIPS" )
{

   TH1F* hi_mom[NVersions];
   TH1F* hi_mult[NVersions];
   
// This part is for plotting model-vs data only.
// It needs to be commented out if both model-vs-data abd regression are wanted (see below).
// TCanvas isn't needed here, because it gets created in the earlier method ! 
//
//   TCanvas* myc = new TCanvas("myc", "", 800, 600);
//   myc->Divide(2,1);
//   myc->cd(1); gPad->SetLeftMargin(0.15);
//   myc->cd(2); gPad->SetLeftMargin(0.15);

   // this part is for plotting both models-vs-data & regression, all in one canvas
   //
   TPad* pad2 = new TPad("pad2","",0.51, 0.01,0.99,0.99);
   pad2->Draw();
   pad2->Divide(1,2);
   pad2->cd(1); gPad->SetLeftMargin(0.15);
   pad2->cd(2); gPad->SetLeftMargin(0.15);

   double ymin_mom = 100000., ymin_mult = 100000. ; // something big... don't know if I can use FLT_MAX
   double ymax_mom = -1., ymax_mult = -1. ;

   for ( int iv=0; iv<NVersions; iv++ )
   {
      std::string histofile = Versions[iv] + "/" + "antiproton" + target + model;
      histofile += ".root";
      TFile* f = new TFile( histofile.c_str() );
      // f->ls();
      hi_mom[iv] = (TH1F*)f->Get("ChargedPionMomentum");
      hi_mult[iv] = (TH1F*)f->Get("NPions");
      hi_mom[iv]->SetLineColor(ColorVersion[iv]);
      hi_mult[iv]->SetLineColor(ColorVersion[iv]);
      hi_mom[iv]->SetLineWidth(2);
      hi_mult[iv]->SetLineWidth(2);
      hi_mom[iv]->GetXaxis()->SetTitle("Charged Pion Momentum (GeV/c)");
      hi_mom[iv]->GetYaxis()->SetTitle("dN/dP (GeV/c)^{-1} / NEvents");
      hi_mult[iv]->GetXaxis()->SetTitle("Number of pions (#pi^{+} + #pi^{-} + #pi^{0})");
      hi_mult[iv]->GetYaxis()->SetTitle("Fraction [%] of events" );
      hi_mom[iv]->GetYaxis()->SetTitleOffset(1.5);
      hi_mult[iv]->GetYaxis()->SetTitleOffset(1.5);
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
//      myc->cd(1);
      pad2->cd(1);
      if ( iv == 0 ) hi_mom[iv]->Draw();
      else hi_mom[iv]->Draw("same");
//      myc->cd(2);
      pad2->cd(2);
      if ( iv == 0 ) hi_mult[iv]->Draw();
      else hi_mult[iv]->Draw("same");      
   
   }

   TLegend* leg1 = new TLegend(0.6, 0.70, 0.9, 0.9);
   leg1->SetTextSize(0.017);
   TLegend* leg2 = new TLegend(0.6, 0.70, 0.9, 0.9);
   leg2->SetTextSize(0.017);
   
   for ( int iv=0; iv<NVersions; iv++ )
   {
      hi_mom[iv]->GetYaxis()->SetRangeUser(ymin_mom,ymax_mom*1.4); // hi[m]->SetTitle("");
      hi_mom[iv]->SetStats(0);
      std::string entry1 = model + ", " + Versions[iv];
      leg1->AddEntry( hi_mom[iv], entry1.c_str(), "L" );
      hi_mult[iv]->GetYaxis()->SetRangeUser(ymin_mult,ymax_mult*1.1); // hi[m]->SetTitle("");
      hi_mult[iv]->SetStats(0);
      std::string entry2 = model + ", " + Versions[iv];
      leg2->AddEntry( hi_mult[iv], entry2.c_str(), "L" );
   }
   
   readAntiProton();

   TGraph*  gr1 = new TGraphErrors(NPointsPbar_PionMom,MomX,MomValue,0,MomError);
   TGraph*  gr2 = new TGraphErrors(NPointsPbar_PionMult,MultX,MultValue,0,MultError);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   
   gr2->SetMarkerColor(4);  gr2->SetMarkerStyle(22);
   gr2->SetMarkerSize(1.6);
   
//   myc->cd(1);
   pad2->cd(1); 
   gr1->Draw("p"); 
   leg1->AddEntry(gr1, "exp.data", "p");
   leg1->Draw();
   leg1->SetFillColor(kWhite);
//   myc->cd(2);
   pad2->cd(2); 
   gr2->Draw("p"); 
   leg2->AddEntry(gr2, "exp.data", "p");
   leg2->Draw();
   leg2->SetFillColor(kWhite);
   
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

int findTarget (std::string target )
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
      return -1;
   }
   
   return TargetID;

}
