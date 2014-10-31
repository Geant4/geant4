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


//const int NModelsMesons=3;
//std::string ModelsMesons[4] = { "CHIPS", "stopping", "BertiniPreCo", "Bertini" };
const int NModelsMesons=1;
std::string ModelsMesons[2] = { "BertiniPreCo", "Bertini" };

//int         ColorModel[4]    = { 1, 2, 3, 6 };
//int         SymbModel[4]     = { 20, 29, 8, 21 };
int         ColorModel[2]    = { 3, 6 };
int         SymbModel[2]     = { 8, 21 };



/*
const int NVersions = 4;
//std::string Versions[4] = { "geant4-09-06-p02", "geant4-09-06-p03", "geant4-10-00", "geant4-10-00-p01" };
//int ColorVersion[5] = { kRed, kGreen, 7, kBlack, 14 };
// std::string Versions[4] = { "geant4-10-00-p01", "geant4-10-00-ref03", "geant4-10-00-ref04", "geant4-10-00-ref05" };
std::string Versions[4] = { "geant4-09-06-p03", "geant4-10-00-p01", "geant4-10-00-p02", "geant4-10-01-b01" };
// int ColorVersion[3] = { kBlack, 14, kYellow };
int ColorVersion[5] = { kRed, kGreen, 7, kBlack, 14 };
// int ColorVersion[5] = { kRed, kMagenta, 7, kBlack, 14 };
// int ColorVersion[6] = { kGreen, kRed, kMagenta, 7, kBlack, 14 };
*/

#include "../test23/shared-root-macros/REGRESSION_TEST.h"


// static int isExpReadLoaded = 0;

#include "../test48/ReadExpData.C"
#include "../test23/shared-root-macros/Chi2Calc.C"

#ifndef G4VAL_PLOTPIMINUSSTOPPING_C
#define G4VAL_PLOTPIMINUSSTOPPING_C
/*
void PlotPiMinusStopping()
{

   plotPiMinusSummary2("C");
   plotPiMinusSummary2("N");
   plotPiMinusSummary2("O");
   plotPiMinusSummary2("Al");
   plotPiMinusSummary2("Cu");
   plotPiMinusSummary2("Ta");
   plotPiMinusSummary2("Pb");
   
   return;

}
*/
/*
void PlotPiMinusStoppingRegre()
{

   plotPiMinusRegressionSummary2( "C", "BertiniPreCo" );
   plotPiMinusRegressionSummary2( "N", "BertiniPreCo" );
   plotPiMinusRegressionSummary2( "O", "BertiniPreCo" );
   plotPiMinusRegressionSummary2( "Al", "BertiniPreCo" );
   plotPiMinusRegressionSummary2( "Cu", "BertiniPreCo" );
   plotPiMinusRegressionSummary2( "Ta", "BertiniPreCo" );
   plotPiMinusRegressionSummary2( "Pb", "BertiniPreCo" );

   return;

}
*/
double calcChi2Madey( std::string target, std::string model, int& NDF )
{

   double chi2 = 0.;
//   int NDF = 0;

   std::string histofile = "piminus" + target + model;
   histofile += ".root";

   std::cout << "About to open file: " << histofile << std::endl;

   TFile* f = new TFile( histofile.c_str() );
   TH1F* hi = (TH1F*)f->Get("NvsT");
   
   TGraphErrors* gdata = getMadeyAsGraph(target);

   chi2 = Chi2( gdata, hi, NDF );

   std::cout << " chi2/NDF = " << (chi2/NDF) << " for model " << model << std::endl;


}

void plotPiMinusDataAndRegression( std::string target )
{

   TCanvas *myc = new TCanvas("myc","",800,600);
   myc->Divide(2,1);
   
   //   if ( isExpReadLoaded <= 0 )
   //   {
   //      gROOT->LoadMacro("ReadExpData.C");
   //      isExpReadLoaded = 1;
   //   }

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

  //   if ( isExpReadLoaded <= 0 )
  //   {
  //      gROOT->LoadMacro("ReadExpData.C");
  //      isExpReadLoaded = 1;
  //   }

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

  //   if ( isExpReadLoaded <= 0 )
  //   {
  //      gROOT->LoadMacro("ReadExpData.C");
  //      isExpReadLoaded = 1;
  //   }

   readMadey();

   TCanvas *myc = new TCanvas("myc","",800,600);
   myc->SetLogx(1);
   myc->SetLogy(1);
   
   drawPiMinus( target );
   // drawPiMinusMC2Data( target );
   
   return;   
}

void plotPiMiBertiniRegression( std::string target )
{

  //   if ( isExpReadLoaded <= 0 )
  //   {
  //      gROOT->LoadMacro("ReadExpData.C");
  //      isExpReadLoaded = 1;
  //   }

   readMadey();
   
   TCanvas* myc = new TCanvas("myc","",800,600);
   

   myc->Divide(2,1);
   
   myc->cd(1);
   gPad->SetLogx();
//   gPad->SetLogy();
   drawPiMinusMC2DataRegression( target, "BertiniPreCo" );

   myc->cd(2);
   gPad->SetLogx();
   drawPiMinusMC2DataRegression( target, "BertiniPreCo" );
   
   myc->cd();

   return;

} 

void plotPiMinusRegressionSummary2( std::string target, std::string model )
{
   
   
  //   if ( isExpReadLoaded <= 0 )
  //   {
  //      gROOT->LoadMacro("ReadExpData.C");
  //      isExpReadLoaded = 1;
  //   }
   
   readMadey();
   
   TCanvas* myc = new TCanvas("myc","",800,600);
   
   myc->Divide(2,1);
   
   myc->cd(1);
   gPad->SetLogx();
   gPad->SetLogy();
   drawPiMinusRegression( target, model );
   
   myc->cd(2);
   gPad->SetLogx();
   gPad->SetLogy();
   drawPiMinusMC2DataRegression( target, model );
   
   std::string output = "pim-" + target + "-Bert-regre.gif";
   myc->Print( output.c_str() );

   myc->cd();
   
   return;

}

void plotPiMinusSummary2( std::string target )
{

  //   if ( isExpReadLoaded <= 0 )
  //   {
  //      gROOT->LoadMacro("ReadExpData.C");
  //      isExpReadLoaded = 1;
  //   }

   readMadey();

   TCanvas *myc = new TCanvas("myc","",800,600);
   
   myc->Divide(2,1);
   
   myc->cd(1);
   gPad->SetLogx(1);
   gPad->SetLogy(1);
   drawPiMinus( target );

   myc->cd(2);
   gPad->SetLogx(1);
   gPad->SetLogy();
   drawPiMinusMC2Data( target );
   
   myc->cd();

   for ( int i=0; i<NModelsMesons; ++i )
     {
       double chi2 = 0.;
       int NDF = 0;
       chi2 = calcChi2Madey( target, ModelsMesons[i], NDF );
     }
   
   
   std::string output = "pim-" + target + "-models.gif";
   myc->Print( output.c_str() );
   
   return;   

}

void plotPiMinusSummary4( std::string target )
{

  //   if ( isExpReadLoaded <= 0 )
  //   {
  //      gROOT->LoadMacro("ReadExpData.C");
  //      isExpReadLoaded = 1;
  //   }

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

void drawPiMinus( std::string target )
{
      
 
   TGraphErrors* gr1 = getMadeyAsGraph(target);
   if ( !gr1 ) 
   {
      return;
   }
    
  /*
   int TargetID = findTarget( target );
         
   if ( TargetID == -1 || TargetID >= NTargetsMadey )
   {
      std::cout << " Invalid Target: " << target << std::endl;
      return;
   }
   */
   
   TH1F* hi[NModelsMesons];
      
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int m=0; m<NModelsMesons; m++ )
   {
     //      std::string histofile = "piminus" + TargetsMadey[TargetID] + ModelsMesons[m];
      std::string histofile = "piminus" + target + ModelsMesons[m];
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
   
//   TLegend* leg = new TLegend(0.55, 0.70, 0.9, 0.9);
   TLegend* leg = new TLegend(0.1, 0.20, 0.6, 0.4);
   leg->SetTextSize(0.04);
   
   for ( int m=0; m<NModelsMesons; m++ )
   {
      // hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*5.); // hi[m]->SetTitle("");
      double delta = 0.5;
      if ( delta > ymax/2. ) delta = ymax/2.;
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax+delta); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelsMesons[m].c_str(), "L" );
   }
      
   //   TGraphErrors*  gr1 = new TGraphErrors(NPointsMadey,KENeut,Value[TargetID],0,Error[TargetID]);
   //gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   //gr1->SetMarkerSize(1.6);
   //
   //    TGraphError* gr1 = getMadeyAsGraph( target );
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
   gr1->GetYaxis()->SetRangeUser(0.2,5.);


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

//   gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->GetYaxis()->SetRangeUser(0.2, 5.);  
   gr1->Draw("apl");    
   
   for ( int m=0; m<NModelsMesons; m++ )
   {
//      gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[m]->GetYaxis()->SetRangeUser( 0.2, 5. );
      gr[m]->Draw("lpsame");
   }
  
   TLegend* leg = new TLegend(0.45, 0.70, 0.95, 0.9);  
//   TLegend* leg = new TLegend(0.1, 0.20, 0.6, 0.4);  
   leg->SetTextSize(0.04); 
   
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
   gr1->GetYaxis()->SetRangeUser(0.2,5.);


   TH1F* hi[NVersions];
         
   float MC2DataX[NPointsMadey];
   float MC2DataY[NPointsMadey];
   float DX[NPointsMadey], DY[NPointsMadey];
   int np =0;
   
   TGraph* gr[NVersions];

   for ( int iver=0; iver<NVersions; iver++ )
   {
//      std::string histofile = "";
// //      if ( iver == 0 )
// //      {
// //         histofile = "piminus" + target + model;
// //      }
// //      else
// //      {
//         histofile = Versions[iver] + "/" + "piminus" + target + model;
// //      }
//      histofile += ".root";

      std::string location = "";
      if ( Versions[iver] == CurrentVersion || Versions[iver] == "." )
      {
         location = "";
      }
      else
      {
         location = regre_test_dir + "/test48/" + Versions[iver] + "/";
      }
      std::string histofile = location + "piminus" + target + model + ".root"; 
      
      std::cout << "MC2Data Regression, histofile: " << histofile << std::endl;
      
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
//   gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->GetYaxis()->SetRangeUser(0.2, 5.);  
   gr1->Draw("apl");    
   
   for ( int iver=0; iver<NVersions; iver++ )
   {
//      gr[iver]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[iver]->GetYaxis()->SetRangeUser( 0.2, 5. );
      gr[iver]->Draw("lpsame");
   }
  
   TLegend* leg = new TLegend(0.45, 0.70, 0.95, 0.9);  
//   TLegend* leg = new TLegend(0.1, 0.20, 0.6, 0.4);  
   leg->SetTextSize(0.04); 
   
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

  //   if ( isExpReadLoaded <= 0 )
  //   {
  //      gROOT->LoadMacro("ReadExpData.C");
  //      isExpReadLoaded = 1;
  //   }

   readMadey();
   
   TCanvas *myc = new TCanvas("myc","",800,600);
   myc->SetLogx(1);
   myc->SetLogy(1);
   
   drawPiMinusRegression( target, model );
   
   return;

}

void plotPiMinusRegressionAllTargets( std::string model="CHIPS" )
{

   
  //   if ( isExpReadLoaded <= 0 )
  //   {
  //      gROOT->LoadMacro("ReadExpData.C");
  //      isExpReadLoaded = 1;
  //   }

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

//      std::string histofile = "";
// //      if ( iver == 0 )
// //      {
// //         histofile = "piminus" + target + model;
// //      }
// //      else
// //      {
//         histofile = Versions[iver] + "/" + "piminus" + target + model;
// //      }
//      histofile += ".root";
      
      std::string location = "";
      if ( Versions[iver] == CurrentVersion || Versions[iver] == "." )
      {
         location = "";
      }
      else
      {
         location = regre_test_dir + "/test48/" + Versions[iver] + "/";
      }
      std::string histofile = location + "piminus" + target + model + ".root"; 

      std::cout << " Regression, histofile: " << histofile << std::endl;
      
      TFile* f = new TFile( histofile.c_str() );
      hi[iver] = (TH1F*)f->Get("NvsT");
      hi[iver]->SetStats(0);
      hi[iver]->SetLineColor(ColorVersion[iver]);
//      hi[iver]->SetLineWidth(2);
      hi[iver]->SetLineWidth(6-iver);
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
   
//   TLegend* leg = new TLegend(0.5, 0.70, 0.99, 0.9);
   TLegend* leg = new TLegend(0.1, 0.20, 0.6, 0.4);
   leg->SetTextSize(0.04);
   
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

/*
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
*/

#endif
