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

// for all types of secondaries except pions
//
int    NPoints         = 0;

float* xF              = 0;

float* dXSecdxF     = 0;
float* err_dXSecdxF = 0;
float* dNdxF           = 0;
float* err_dNdxF       = 0;
float* pT              = 0;
float* err_pT          = 0;
float* pT2             = 0;
float* err_pT2         = 0;

float* XSec            = 0;
float* err_XSec        = 0;

// for pions (pi+ & pi-)
//
float* Y =0;
float* dNdY[2]     = { 0, 0 };
float* err_dNdY[2] = { 0, 0 };


const int NModels = 3;
std::string ModelName[3]  = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };
int         ColorModel[3] = { 14, 7, kRed }; // 14 = grey, 7 = light "sky"-blue
//
//int         ColorModel[4] = { kGreen, 7, kRed, kBlack }; // 14 = grey, 7 = light "sky"-blue
// int         ColorModel[4] = { kMagenta, 7, kRed, kBlack }; // 14 = grey, 7 = light "sky"-blue
int         SymbModel[3]     = { 8, 21, 23 };

const int NVersions = 4;
// ---> std::string Versions[3] = { "geant4-10-00-b01", "geant4-09-06-ref03", "geant4-09-06-p01" };
std::string Versions[4] = { "geant4-09-06-ref07", "geant4-10-00-b01", "geant4-09-06-ref03", "geant4-09-06-p01" };
// std::string Versions[2] = { "geant4-10-00-b01", "geant4-09-06-ref03" };
std::string CurrentVersion = "geant4-10-00-b01";
int ColorVersion[5] = { kRed, kGreen, 7, kBlack, 14 };

void readIntegratedSpectra( std::string beam, std::string target, std::string secondary )
{

   // std::string dirname = "../test23/na49-exp-data/";
   std::string dirname = "./na49/";
   
   std::string filename = beam + "_" + target + "_" + secondary;
   
   std::string filenam1;
   
   if ( secondary == "pion" )
   { 
      filenam1 = filename + "_dndy.dat";
   }
   else
   {
      filename1 = filename + "_integr.dat";
   }
   
   std::string file1 = dirname + filename1;
   
   ifstream infile1;
   infile1.open( file1.c_str() );

   infile1 >> NPoints;
   
   if ( secondary == "pion" )
   {
      if (Y) delete Y;
      Y = new float[NPoints];
      if ( dNdY[0] = 0 ) delete dNdY[0];
      if ( dNdY[1] = 0 ) delete dNdY[1];
      dNdY[0] = new float[NPoints]; // pi+
      dNdY[1] = new float[NPoints]; // pi-
      if ( err_dNdY[0] = 0 ) delete err_dNdY[0];
      if ( err_dNdY[1] = 0 ) delete err_dNdY[1];      
      err_dNdY[0] = new float[NPoints];
      err_dNdY[1] = new float[NPoints];
      for ( int i=0; i<NPoints; i++ )
      {
         infile1 >> Y[i] >> dNdY[0][i] >> dNdY[1][i];
	 err_dNdY[0][i] = dNdY[0][i]*0.02; // errors are said to be 2%
	 err_dNdY[1][i] = dNdY[1][i]*0.02;
      } 
      return;     
   }
   

   if (xF) delete [] xF;
   xF = new float[NPoints];
   
  
   if (dNdxF) delete dNdxF;
   dNdxF = new float[NPoints];
   
   if (err_dNdxF) delete err_dNdxF;
   err_dNdxF = new float[NPoints];
   
   if ( secondary != "neutron" )
   {
     if  (dXSecdxF) delete [] dXSecdxF;
     dXSecdxF = new float[NPoints];
     if  (err_dXSecdxF) delete [] err_dXSecdxF;
     err_dXSecdxF = new float[NPoints];
     if (pT) delete [] pT;
     pT = new float[NPoints];
     if ( err_pT ) delete [] err_pT;
     err_pT = new float[NPoints];
     if (pT2) delete [] pT2;
     pT2 = new float[NPoints];
     if ( err_pT2 ) delete [] err_pT2;
     err_pT2 = new float[NPoints];
   }
   
   std::cout << " NPoints = " << NPoints << std::endl;

   for ( int i=0; i<NPoints; i++ )
   {
      if ( secondary == "neutron" )
      {
         infile1 >> xF[i] >> dNdxF[i] >> err_dNdxF[i];
	 err_dNdxF[i] *= 0.01;
      }
      else
      {
         infile1 >> xF[i] >> dXSecdxF[i] >> err_dXSecdxF[i]  >> dNdxF[i] >> err_dNdxF[i]  
	         >> pT[i] >> err_pT[i] >> pT2[i] >> err_pT2[i];         
	 err_dXSecdxF[i]    *= 0.01;
	 err_dNdxF[i]       *= 0.01;
	 err_pT[i]          *= 0.01;
	 err_pT2[i]         *= 0.01;
      }
   }
   
/*
   if ( secondary == "neutron" ) return;
   
   infile1.close();
   
   std::string filename2 = filename + "_pt_xf.dat";
   
   std::string file2 = dirname + filename2;
     
   ifstream infile2;
   infile2.open( file2.c_str() );

   infile2 >> NPoints;
   
   if (xF) delete xF;
   xF = new float[NPoints];
   
   if (pT) delete pT;
   pT = new float[NPoints];
   
   if (XSec) delete XSect;
   XSec = new float[NPoints];
   
   if (err_XSec) delete err_XSec;
   err_XSec = new float[NPoints];

   for ( int i=0; i<NPoints; i++ )
   {
      infile2 >> xF[i] >> pT[i] >> XSec[i] >> err_XSec[i];
   }
*/
   return;

}

void plot_dNdxF( std::string beam, std::string target )
{

   TCanvas *myc = new TCanvas("myc","",1200,800);
   myc->Divide(3,2);

   myc->cd(1);
   drawIntegratedSpectrum( beam, target, "proton", "dNdxF" );
   
   myc->cd(2);
   drawIntegratedSpectrum( beam, target, "antiproton", "dNdxF" );

   myc->cd(3);
   drawIntegratedSpectrum( beam, target, "neutron", "dNdxF" );

   myc->cd(4);
   drawIntegratedSpectrum( beam, target, "piplus", "dNdxF" );

   myc->cd(5);
   drawIntegratedSpectrum( beam, target, "piminus", "dNdxF" );


   return;

}

void plot_dNdxF_Regre( std::string beam, std::string target, std::string model )
{

   TCanvas *myc = new TCanvas("myc","",1200,800);
   myc->Divide(3,2);
   
   myc->cd(1);
   drawIntergratedSpectrumRegre( beam, target, "proton", "dNdxF", model );

   myc->cd(2);
   drawIntergratedSpectrumRegre( beam, target, "antiproton", "dNdxF", model );

   myc->cd(3);
   drawIntergratedSpectrumRegre( beam, target, "neutron", "dNdxF", model );

   myc->cd(4);
   drawIntergratedSpectrumRegre( beam, target, "piplus", "dNdxF", model );

   myc->cd(5);
   drawIntergratedSpectrumRegre( beam, target, "piminus", "dNdxF", model );
   
   return;

} 

void plot_pT( std::string beam, std::string target )
{

   TCanvas *myc = new TCanvas("myc","",800,800);
   myc->Divide(2,2);

   myc->cd(1);
   drawIntegratedSpectrum( beam, target, "proton", "pT" );
   
   myc->cd(2);
   drawIntegratedSpectrum( beam, target, "antiproton", "pT" );

   myc->cd(3);
   drawIntegratedSpectrum( beam, target, "piplus", "pT" );

   myc->cd(4);
   drawIntegratedSpectrum( beam, target, "piminus", "pT" );


   return;

}

void plot_pT_Regre( std::string beam, std::string target, std::string model )
{

   TCanvas *myc = new TCanvas("myc","",800,800);
   myc->Divide(2,2);
   
   myc->cd(1);
   drawIntergratedSpectrumRegre( beam, target, "proton", "pT", model );

   myc->cd(2);
   drawIntergratedSpectrumRegre( beam, target, "antiproton", "pT", model );

   myc->cd(3);
   drawIntergratedSpectrumRegre( beam, target, "piplus", "pT", model );

   myc->cd(4);
   drawIntergratedSpectrumRegre( beam, target, "piminus", "pT", model );
   
   return;

} 


void drawIntegratedSpectrum( std::string beam, std::string target, std::string secondary, std::string histo )
{

   readIntegratedSpectra( beam, target, secondary );

   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;

   float* Value = new float[NPoints];
   float* Error = new float[NPoints];
   
   for ( int i=0; i<NPoints; i++ )
   {
      if ( histo == "dNdxF" )
      {
         Value[i] = dNdxF[i];
	 Error[i] = err_dNdxF[i];
      }
      else if ( histo == "pT" )
      {
         Value[i] = pT[i];
	 Error[i] = err_pT[i];
      }
      else if ( histo == "pT2" )
      {
         Value[i] = pT2[i];
	 Error[i] = err_pT2[i];
      }
      if ( Value[i]+Error[i] > ymax ) ymax = Value[i] + Error[i];
      if ( Value[i]-Error[i] < ymin ) ymin = Value[i] - Error[i];
      if ( ymin < 0. ) ymin = 0.; 
   }
 
   TH1F* hi[NModels];
   std::string YTitle;
   if ( histo == "dNdxF" )
   {
      YTitle = "dN/dxF";
   }
   else if ( histo == "pT" )
   {
      YTitle = "dpT/dxF, GeV/c";
   }
   else if ( histo == "pT2" )
   {
      YTitle = "dpT^2/dxF, (GeV/c)^2";
   }
   
   for ( int m=0; m<NModels; m++ )
   {

      std::string histofile = "./na49-histo/" + beam + target + "158.0GeV" + ModelName[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_" + histo ;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) 
      {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
         hi[m]->Draw();
	 hi[m]->GetXaxis()->SetTitle("xF");
	 hi[m]->GetYaxis()->SetTitle( YTitle.c_str() );
	 hi[m]->GetYaxis()->SetTitleOffset(1.5);
      }
      else hi[m]->Draw("same");     

   }
   
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
   
   for ( int m=0; m<NModels; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelName[m].c_str(), "L" );
   }
     
   TGraph* gr = new TGraphErrors( NPoints, xF, Value, 0, Error );

   // gr->GetYaxis()->SetRangeUser( 0., 2.5 );
   // gr->GetXaxis()->SetRangeUser( -0.3, 0.4 );
   // gr->SetRangeUser( 0., 2.5 );
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
    
   gr->Draw("p");
      
   leg->AddEntry( gr, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);   
   
   return;

}

void drawIntSpectrumMC2Data( std::string beam, std::string target, std::string secondary, std::string histo )
{

   readIntegratedSpectra( beam, target, secondary );

   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;

   float* Value = new float[NPoints];
   float* Error = new float[NPoints];
   
   for ( int i=0; i<NPoints; i++ )
   {
      if ( histo == "dNdxF" )
      {
	 Value[i] = dNdxF[i];
	 Error[i] = err_dNdxF[i];
      }
      else if ( histo == "pT" )
      {
	 Value[i] = pT[i];
	 Error[i] = err_pT[i];
      }
      else if ( histo == "pT2" )
      {
	 Value[i] = pT2[i];
	 Error[i] = err_pT2[i];
      }
   }

   TH1F* hi[NModels];
   
   float* MC2DataX = new float[NPoints];
   float* MC2DataY = new float[NPoints];
   float* DX = new float[NPoints];
   float* DY = new float[NPoints];
   int np = 0;
   
   TGraph* gr[NModels];

   for ( int m=0; m<NModels; m++ )
   {

      std::string histofile = "./na49-histo/" + beam + target + "158.0GeV" + ModelName[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_" + histo ;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
      int nx = hi[m]->GetNbinsX();
      
      np=0;
      for ( int k=0; k<=nx; k++ ) 
      {
         double xx1 = hi[m]->GetBinLowEdge(k);
	 double xx2 = hi[m]->GetBinWidth(k);
	 for (int kk=0; kk<NPoints; kk++ )
	 {
	    if ( xx1 < xF[kk] && xx1+xx2 > xF[kk] )
	    {
	       double yy = hi[m]->GetBinContent(k);
	       MC2DataX[np] = xF[kk];
	       DX[np] = 0.;
	       
	       MC2DataY[np] = yy / Value[kk];
	       // also need error calc here !...
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
      gr[m]->GetXaxis()->SetTitle("xF");
      gr[m]->GetYaxis()->SetTitle("MC/Data");
      gr[m]->SetMarkerColor(ColorModel[m]);  
      gr[m]->SetMarkerStyle(SymbModel[m]);
      gr[m]->SetMarkerSize(1.6);
      
      // if ( m==0 ) gr1->SetTitle(hi[m]->GetTitle());
   
   }
   
   for ( int i=0; i<NPoints; i++ )
   {
      Error[i] /= Value[i];
      Value[i] = 1.;
      if ( Value[i]+Error[i] > ymax ) ymax = Value[i] + Error[i];
      if ( Value[i]-Error[i] < ymin ) ymin = Value[i] - Error[i];
      if ( ymin < 0. ) ymin = 0.; 
   }
   TGraph* gr1 = new TGraphErrors( NPoints, xF, Value, 0, Error );
   gr1->GetYaxis()->SetRangeUser( ymin, ymax );
   // gr1->GetXaxis()->SetRangeUser( -0.3, 0.4 );
   gr1->SetMarkerColor(kBlue);
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.5);
   gr1->GetXaxis()->SetTitle("xF");
   if ( histo == "dNdxF" )
   {
      gr1->GetYaxis()->SetTitle("MC/Data (dN/dxF)");
   }
   else if ( histo == "pT" )
   {
      gr1->GetYaxis()->SetTitle("MC/Data (Average pT vs xF)" );
   }


   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->Draw("apl");    
   
   for ( int m=0; m<NModels; m++ )
   {
      // gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[m]->Draw("lpsame");
   }
  
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);   
   
   for ( int m=0; m<NModels; m++ )
   {
      leg->AddEntry( gr[m], ModelName[m].c_str(), "p" );
   }

   leg->AddEntry( gr1, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);
   
   return;

}


void drawIntergratedSpectrumRegre( std::string beam, std::string target, std::string secondary, 
                                   std::string histo, std::string model )
{

   readIntegratedSpectra( beam, target, secondary );

   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;

   float* Value = new float[NPoints];
   float* Error = new float[NPoints];
   
   for ( int i=0; i<NPoints; i++ )
   {
      if ( histo == "dNdxF" )
      {
         Value[i] = dNdxF[i];
	 Error[i] = err_dNdxF[i];
      }
      else if ( histo == "pT" )
      {
         Value[i] = pT[i];
	 Error[i] = err_pT[i];
      }
      else if ( histo == "pT2" )
      {
         Value[i] = pT2[i];
	 Error[i] = err_pT2[i];
      }
      if ( Value[i]+Error[i] > ymax ) ymax = Value[i] + Error[i];
      if ( Value[i]-Error[i] < ymin ) ymin = Value[i] - Error[i];
      if ( ymin < 0. ) ymin = 0.; 
   }
   
   TH1F* hi[NVersions];
   std::string YTitle;
   if ( histo == "dNdxF" )
   {
      YTitle = "dN/dxF";
   }
   else if ( histo == "pT" )
   {
      YTitle = "dpT/dxF, GeV/c";
   }
   else if ( histo == "pT2" )
   {
      YTitle = "dpT^2/dxF, (GeV/c)^2";
   }
   
   for ( int iv=0; iv<NVersions; iv++ )
   {
      
      std::string histofile = "";
      
      if ( Versions[iv] == CurrentVersion )
      {
         histofile = "./na49-histo/" + beam + target + "158.0GeV" + model;
      }
      else
      {
         histofile = Versions[iv] + "/na49-histo/" + beam + target + "158.0GeV" + model;
      }
      histofile += ".root";
      
      std::cout << " Regression, histofile: " << histofile << std::endl;
      
      TFile* f = new TFile( histofile.c_str() );
      std::string histoname = secondary + "_" + histo ;
      
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
	 hi[iv]->GetYaxis()->SetTitle( YTitle.c_str() );
	 hi[iv]->GetYaxis()->SetTitleOffset(1.5);
      }
      else hi[iv]->Draw("same");     

   }

   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
   leg->SetTextSize(0.017);
   
   for ( int iv=0; iv<NVersions; iv++ )
   {
      hi[iv]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); // hi[m]->SetTitle("");
      // std::string htitle = "pi- on " + target + ", " + model;
      // hi[iv]->SetTitle( htitle.c_str() );
      // leg->AddEntry( hi[iv], entry.c_str(), "L" );
      leg->AddEntry( hi[iv], Versions[iv].c_str(), "L" );
   }

   TGraph* gr = new TGraphErrors( NPoints, xF, Value, 0, Error );

   // gr->GetYaxis()->SetRangeUser( 0., 2.5 );
   // gr->GetXaxis()->SetRangeUser( -0.3, 0.4 );
   // gr->SetRangeUser( 0., 2.5 );
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
    
   gr->Draw("p");
      
   leg->AddEntry( gr, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);   

   return;

}
