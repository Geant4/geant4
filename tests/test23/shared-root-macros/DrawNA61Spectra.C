
#include <iostream>
#include <sstream>
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

void drawMomSpectrum( std::string beam, std::string target, 
                      std::string secondary, 
		      std::string sec_angle_min, std::string sec_angle_max)
{

/*
   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA61Data.C");
      // gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA61UtilLoaded = 1;
   }
*/   
   readMomSpectrum( beam, target, secondary, sec_angle_min, sec_angle_max );
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
   for ( int ip=0; ip<NPoints; ip++ )
   {
      if ( (SecSigma[ip]+SecESigma[ip]) > ymax ) ymax = SecSigma[ip]+SecESigma[ip];
      if ( (SecSigma[ip]-SecESigma[ip]) < ymin ) ymin = SecSigma[ip]-SecESigma[ip];
      if ( ymin < 0. ) ymin = 0.;
   }
   
   TH1F* hi[NModels];
  
   TString YTitle = "d#sigma/dp, [mb/(GeV/c)]";
   
   for ( int m=0; m<NModels; m++ )
   {
   
      std::string histofile = "./na61-histo/" + beam + target + "31.0GeV" + ModelName[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );
      
      std::string histoname = secondary + "Mult_" + sec_angle_min + "_" + sec_angle_max;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      // FIXME !!!
      // hi[m]->GetXaxis()->SetTitle("Kinetic energy of secondary neutron (MeV)");
      // hi[m]->GetYaxis()->SetTitle("Number of neutrons per MeV");
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
         hi[m]->Draw();
	 hi[m]->GetXaxis()->SetTitle("p, GeV/c");
	 hi[m]->GetYaxis()->SetTitle(YTitle);
	 hi[m]->GetYaxis()->SetTitleOffset(1.5);
      }
      else hi[m]->Draw("same");     
   }
   
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
   
   for ( int m=0; m<NModels; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelName[m].c_str(), "L" );
   }
      
   TGraph* gr = new TGraphErrors( NPoints, SecMom, SecSigma, 0, SecESigma );
   // gr->GetYaxis()->SetRangeUser( 0., 0.02 );
   // gr->SetRangeUser( ymin, ymax*1.5 );
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
   
   gr->Draw("p");
   
   leg->AddEntry( gr, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);

   return ;

}

void drawKPlus2PiPlusRatio( std::string beam, std::string target, 
		            std::string sec_angle_min, std::string sec_angle_max)
{

/*   
   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA61Data.C");
      // gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA61UtilLoaded = 1;
   }
*/

   readKPlus2PiPlusRatio( beam, target, sec_angle_min, sec_angle_max );
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
   for ( int ip=0; ip<NPoints; ip++ )
   {
      if ( (K2PiRatio[ip]+K2PiERatio[ip]) > ymax ) ymax = K2PiRatio[ip]+K2PiERatio[ip];
      if ( (K2PiRatio[ip]-K2PiERatio[ip]) < ymin ) ymin = K2PiRatio[ip]-K2PiERatio[ip];
      if ( ymin < 0. ) ymin = 0.;
   }
   
   TH1F* hi[NModels];
  
   TString YTitle = "K^{+}/#pi^{+} ratio";
   
   for ( int m=0; m<NModels; m++ )
   {
   
      std::string histofile = "./na61-histo/" + beam + target + "31.0GeV" + ModelName[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );
      
      std::string histoname = "kplus2piplus_" + sec_angle_min + "_" + sec_angle_max;
      std::cout << "histoname = " << histoname << std::endl;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      // FIXME !!!
      // hi[m]->GetXaxis()->SetTitle("Kinetic energy of secondary neutron (MeV)");
      // hi[m]->GetYaxis()->SetTitle("Number of neutrons per MeV");
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
         hi[m]->Draw();
	 hi[m]->GetXaxis()->SetTitle("p, GeV/c");
	 hi[m]->GetYaxis()->SetTitle(YTitle);
	 hi[m]->GetYaxis()->SetTitleOffset(1.5);
      }
      else hi[m]->Draw("same");     
   }
   
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
   
   for ( int m=0; m<NModels; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelName[m].c_str(), "L" );
   }
      
   TGraph* gr = new TGraphErrors( NPoints, SecMom, K2PiRatio, 0, K2PiERatio );
   // gr->GetYaxis()->SetRangeUser( 0., 0.02 );
   // gr->SetRangeUser( ymin, ymax*1.5 );
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
   
   gr->Draw("p");
   
   leg->AddEntry( gr, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);

   return ;

}

void drawKPlus2PiPlusRatioRegre( std::string beam, std::string target, 
		                 std::string sec_angle_min, std::string sec_angle_max,
			         std::string model )
{

/*
   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA61Data.C");
      // gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA61UtilLoaded = 1;
   }
*/   

   readKPlus2PiPlusRatio( beam, target, sec_angle_min, sec_angle_max );
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
   for ( int ip=0; ip<NPoints; ip++ )
   {
      if ( (K2PiRatio[ip]+K2PiERatio[ip]) > ymax ) ymax = K2PiRatio[ip]+K2PiERatio[ip];
      if ( (K2PiRatio[ip]-K2PiERatio[ip]) < ymin ) ymin = K2PiRatio[ip]-K2PiERatio[ip];
      if ( ymin < 0. ) ymin = 0.;
   }
   
   TH1F* hi[NModels];
  
   TString YTitle = "K^{+}/#pi^{+} ratio";
   
   for ( int m=0; m<NVersions-1; m++ )
   {
   
      std::string histofile = "";
      if ( Versions[m] == CurrentVersion )
      {
         histofile = "./na61-histo/" + beam + target + "31.0GeV" + model;
      }
      else
      {
         histofile = Versions[m] + "/na61-histo/" + beam + target + "31.0GeV" + model;
      }
      histofile += ".root";

      TFile* f = new TFile( histofile.c_str() );
      
      std::string histoname = "kplus2piplus_" + sec_angle_min + "_" + sec_angle_max;
      std::cout << "histoname = " << histoname << std::endl;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorVersion[m]);
      hi[m]->SetLineWidth(2);
      // FIXME !!!
      // hi[m]->GetXaxis()->SetTitle("Kinetic energy of secondary neutron (MeV)");
      // hi[m]->GetYaxis()->SetTitle("Number of neutrons per MeV");
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
         hi[m]->Draw();
	 hi[m]->GetXaxis()->SetTitle("p, GeV/c");
	 hi[m]->GetYaxis()->SetTitle(YTitle);
	 hi[m]->GetYaxis()->SetTitleOffset(1.5);
      }
      else hi[m]->Draw("same");     
   }
   
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
   
   for ( int m=0; m<NVersions-1; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], Versions[m].c_str(), "L" );
   }
      
   TGraph* gr = new TGraphErrors( NPoints, SecMom, K2PiRatio, 0, K2PiERatio );
   // gr->GetYaxis()->SetRangeUser( 0., 0.02 );
   // gr->SetRangeUser( ymin, ymax*1.5 );
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
   
   gr->Draw("p");
   
   leg->AddEntry( gr, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);

   return ;

}

void drawMomSpectrumRegre( std::string beam, std::string target, 
                           std::string secondary, 
		           std::string sec_angle_min, std::string sec_angle_max,
			   std::string model )
{

/*
   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA61Data.C");
      // gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA61UtilLoaded = 1;
   }
*/

   readMomSpectrum( beam, target, secondary, sec_angle_min, sec_angle_max );
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
   for ( int ip=0; ip<NPoints; ip++ )
   {
      if ( (SecSigma[ip]+SecESigma[ip]) > ymax ) ymax = SecSigma[ip]+SecESigma[ip];
      if ( (SecSigma[ip]-SecESigma[ip]) < ymin ) ymin = SecSigma[ip]-SecESigma[ip];
      if ( ymin < 0. ) ymin = 0.;
   }
   
   TH1F* hi[NVersions];
  
   TString YTitle = "d#sigma/dp, [mb/(GeV/c)]";

   for ( int iv=0; iv<NVersions; iv++ )
   {
      
      std::string histofile = "";
      
      if ( Versions[iv] == CurrentVersion )
      {
         histofile = "./na61-histo/" + beam + target + "31.0GeV" + model;
      }
      else
      {
         histofile = Versions[iv] + "/na61-histo/" + beam + target + "31.0GeV" + model;
      }
      histofile += ".root";
      
      std::cout << " Regression, histofile: " << histofile << std::endl;
      
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "Mult_" + sec_angle_min + "_" + sec_angle_max;
      
      std::cout << "Regression, histoname: " << histoname << std::endl;
      
      hi[iv] = (TH1F*)f->Get( histoname.c_str() );
      hi[iv]->SetStats(0);
      hi[iv]->SetLineColor(ColorVersion[iv]);
      hi[iv]->SetLineWidth(2);
      int nx = hi[iv]->GetNbinsX();
      for (int k=1; k <= nx; k++) 
      {
	double yy = hi[iv]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( iv == 0 ) 
      {
         hi[iv]->Draw();
	 hi[iv]->GetXaxis()->SetTitle("xF");
	 hi[iv]->GetYaxis()->SetTitle( YTitle );
	 hi[iv]->GetYaxis()->SetTitleOffset(1.5);
      }
      else hi[iv]->Draw("same");     

   }
   
   std::cout << "now do TLegend" << std::endl;

   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
   
   for ( int iv=0; iv<NVersions; iv++ )
   {
      hi[iv]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // hi[m]->SetTitle("");
      leg->AddEntry( hi[iv], Versions[iv].c_str(), "L" );
   }
      
   TGraph* gr = new TGraphErrors( NPoints, SecMom, SecSigma, 0, SecESigma );
   // gr->GetYaxis()->SetRangeUser( 0., 0.02 );
   // gr->SetRangeUser( ymin, ymax*1.5 );
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
   
   gr->Draw("p");
   
   leg->AddEntry( gr, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);

   return ;

}
