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
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"

int    NPoints    = 0;

float* SecMom     = 0;

float* SecSigma   = 0; 
float* SecESigma  = 0;

float* K2PiRatio  = 0;
float* K2PiERatio = 0;

const int NModels = 3;
std::string ModelName[3] = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };
// int         ColorModel[2]    = { 14, 7 }; // 14 = grey, 7 = light "sky"-blue
int         ColorModel[3] = { kGreen, 7, kRed }; // 14 = grey, 7 = light "sky"-blue

const int NVersions = 4;
std::string Versions[4] = { "geant4-09-06-ref07", "geant4-10-00-b01", "geant4-09-06-ref03", "geant4-09-06-p01" };
// std::string Versions[2] = { "geant4-10-00-b01", "geant4-09-06-ref03" };
std::string CurrentVersion = "geant4-10-00-b01";
int ColorVersion[5] = { kRed, kGreen, 7, kBlack, 14 };


void readMomSpectrum( std::string beam, std::string target, 
                      std::string secondary, 
		      std::string sec_angle_min, std::string sec_angle_max )
{

   std::string dirname = "./na61/";
   std::string filename = beam + "_" + target + "_" + secondary;
   filename += "_" + sec_angle_min + "_" + sec_angle_max + ".dat";
   std::string file = dirname + filename;
   
   // std::cout << " file = " << file << std::endl;
   
   ifstream infile;
   infile.open( file.c_str() );
   
   infile >> NPoints;
   
   // std::cout << "NPoints = " << NPoints << std::endl;
   
   if ( SecMom ) delete [] SecMom;
   SecMom = new float[NPoints];
   
   if ( SecSigma ) delete [] SecSigma;
   SecSigma = new float[NPoints];
   
   if ( SecESigma ) delete [] SecESigma;
   SecESigma = new float[NPoints];
   
   float pmin, pmax;
   for ( int i=0; i<NPoints; i++ )
   {  
      infile >> pmin >> pmax >> SecSigma[i] >> SecESigma[i] ;
      SecMom[i] = 0.5 * ( pmin + pmax );
      // std::cout << SecMom[i] << " " << SecSigma[i] << " " << SecESigma[i] << std::endl;
      if ( secondary != "proton" )
      {
         // scale to the total xsec (proton data come normalized)
	 SecSigma[i]  /= 229.3 ; // 229.3 is production sxec, while inelastic xsec = 257.2 
	                         // which one should be used for normalization ???        
         SecESigma[i] /= 229.3;
      }
   }
   
   infile.close();
   
   return;

}

void readKPlus2PiPlusRatio( std::string beam, std::string target, 
		            std::string sec_angle_min, std::string sec_angle_max )
{

   std::string dirname = "./na61/";
   std::string filename = beam + "_" + target + "_kplus2piplus";
   filename += "_" + sec_angle_min + "_" + sec_angle_max + ".dat";
   std::string file = dirname + filename;
   
   // std::cout << " file = " << file << std::endl;
   
   ifstream infile;
   infile.open( file.c_str() );
   
   infile >> NPoints;
   
   if ( SecMom ) delete [] SecMom;
   SecMom = new float[NPoints];

   if ( K2PiRatio ) delete [] K2PiRatio;
   K2PiRatio = new float[NPoints];
   
   if ( K2PiERatio ) delete [] K2PiERatio;
   K2PiERatio = new float[NPoints];
   
   float pmin, pmax;
   for ( int i=0; i<NPoints; i++ )
   {  
      infile >> pmin >> pmax >> K2PiRatio[i] >> K2PiERatio[i] ;
      SecMom[i] = 0.5 * ( pmin + pmax );
      // std::cout << SecMom[i] << " " << K2PiRatio[i] << " " << K2PiERatio[i] << std::endl;
   }
   
   infile.close();

   return;

}


void drawMomSpectrum( std::string beam, std::string target, 
                      std::string secondary, 
		      std::string sec_angle_min, std::string sec_angle_max)
{

   
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

void plotSecondarySum8( std::string secondary )
{

   TCanvas *myc1 = new TCanvas("myc1","",900,900);
   myc1->Divide(2,2);
   
   myc1->cd(1);
   gPad->SetLogy(); //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "0", "20" );

   myc1->cd(2);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "20", "40" );

   myc1->cd(3);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "40", "60" );

   myc1->cd(4);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "60", "100" );

   TCanvas *myc2 = new TCanvas("myc2","",900,900);
   myc2->Divide(2,2);
   
   myc2->cd(1);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "100", "140" );

   myc2->cd(2);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "140", "180" );

   myc2->cd(3);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "180", "240" );
   
   if ( secondary == "proton" ) return;

   myc2->cd(4);
   //gPad->SetLogx();
   drawMomSpectrum( "proton", "C", secondary, "240", "300" );

   return;

}

void plotSecondarySum8Regre( std::string secondary, std::string model )
{

   TCanvas *myc1 = new TCanvas("myc1","",900,900);
   myc1->Divide(2,2);
   
   myc1->cd(1);
   gPad->SetLogy(); //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "0", "20", model );

   myc1->cd(2);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "20", "40", model );

   myc1->cd(3);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "40", "60", model );

   myc1->cd(4);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "60", "100", model );

   TCanvas *myc2 = new TCanvas("myc2","",900,900);
   myc2->Divide(2,2);
   
   myc2->cd(1);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "100", "140", model );

   myc2->cd(2);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "140", "180", model );

   myc2->cd(3);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "180", "240", model );
   
   if ( secondary == "proton" ) return;

   myc2->cd(4);
   //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "240", "300", model );

   return;

}
void drawMomSpectrumRegre( std::string beam, std::string target, 
                           std::string secondary, 
		           std::string sec_angle_min, std::string sec_angle_max,
			   std::string model )
{

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
