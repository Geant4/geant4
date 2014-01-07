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

//const int NModels = 2;
// std::string ModelName[4]  = { "qgsp_ftfp_bert", "ftfp_bert", "ftfp", "qgsp" };  
//std::string ModelName[2]  = { "NuBeam", "ftfp_bert" };  
// std::string ModelName[5]  = { "NuBeam", "qgsp_bert", "ftfp_bert", "NuBeam-with-res-decays", "qgsp-g4lund-str-fragm"};  
// std::string ModelName[4]  = { "ftfp", "qgsp", "ftfp_bert", "qgsp_ftfp_bert" };  
// const int NModels = 3;
// std::string ModelName[3]  = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };
//int         ColorModel[5] = { kMagenta, 7, kRed, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue

// static int isReadNA49Loaded = 0;


void drawIntegratedSpectrum( std::string beam, std::string target,  
                             std::string secondary, std::string histo )
{

//   std::cout << "About to load Reader" << std::endl;
   
//   if ( isReadNA49Loaded <= 0 )
//   {
//      gROOT->LoadMacro("ReadNA49Data.C");
//      isReadNA49Loaded = 1;
//   }
   
//   std::cout << "Reader loaded" << std::endl;
 
   
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
      
      if ( secondary == "piplus" || secondary == "piminus"  || secondary == "antiproton" )
      {
         hi[m]->GetXaxis()->SetRangeUser(-0.55,0.55);
      }
      
/*
      if ( histo == "pT" && m<3 )
      {
         hi[m]->Scale(32.);
      }
*/
      
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel[m]);
      hi[m]->SetLineWidth(2);
      if ( m == 0 ) hi[m]->SetLineWidth(3.5);
      
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
   
   TLegend* leg = new TLegend(0.55, 0.65, 0.9, 0.9);
   leg->SetTextSize(0.025);
   
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
   gr->SetMarkerSize(1.0);
   // gr->SetMarkerSize(1.5);
    
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
      // gr[m]->SetTitle(hi[m]->GetTitle());
      gr[m]->SetTitle("");
      gr[m]->GetXaxis()->SetTitle("xF");
      gr[m]->GetYaxis()->SetTitle("MC/Data");
      gr[m]->SetMarkerColor(ColorModel[m]);  
      gr[m]->SetMarkerStyle(SymbModel[m]);
      gr[m]->SetMarkerSize(1.6);
      if ( secondary == "piplus" || secondary == "piminus" || secondary == "antiproton" )
      {
	 gr[m]->GetXaxis()->SetRangeUser(-0.55,0.55);
      }
      
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
   gr1->GetYaxis()->SetRangeUser( 0.5, 1.5 );
   // gr1->GetYaxis()->SetRangeUser( ymin, ymax );
   // gr1->GetXaxis()->SetRangeUser( -0.3, 0.4 );
   gr1->SetMarkerColor(kBlue);
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.0);
   // gr1->SetMarkerSize(1.5);
   gr1->SetTitle("");
   
   TAxis* xaxis = gr1->GetXaxis();
   
   xaxis->SetTitle("xF");
   
   if ( histo == "dNdxF" )
   {
      gr1->GetYaxis()->SetTitle("MC/Data (dN/dxF)");
      if ( secondary == "piplus" || secondary == "piminus" || secondary == "antiproton" )
      {
         xaxis->SetLimits( -0.55, 0.55 );
      }
   }
   else if ( histo == "pT" )
   {
      gr1->GetYaxis()->SetTitle("MC/Data (Average pT vs xF)" );
   }


   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->Draw("apl");
   // gr1->Draw("AC*");    
   
   for ( int m=0; m<NModels; m++ )
   {
      // gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[m]->Draw("lpsame");
   }
  
   TLegend* leg = new TLegend(0.1, 0.70, 0.4, 0.9);   
   
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

void drawStdDiffGraph( std::string beam, std::string target, std::string secondary, std::string histo )
{

   readIntegratedSpectra( beam, target, secondary );

   //double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   //double ymax = -1. ;
   float ymin = 10000.; // something big... don't know if I can use FLT_MAX
   float ymax = -1. ;

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
   
   float* StdDiffX = new float[NPoints];
   float* StdDiffY = new float[NPoints];
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
	       double yerror = hi[m]->GetBinError(k);
	       StdDiffX[np] = xF[kk];
	       StdDiffY[np] = ( yy - Value[kk]) / sqrt( Error[kk]*Error[kk] + yerror*yerror );
	       if ( StdDiffY[np] > ymax ) ymax = StdDiffY[np];
	       if ( StdDiffY[np] < ymin ) ymin = StdDiffY[np];
	       np++;
	       break;
	    }
	 }
      }
            
      gr[m] = new TGraph( np, StdDiffX, StdDiffY );
      gr[m]->SetTitle(hi[m]->GetTitle());
      // gr[m]->SetTitle("");
      gr[m]->GetXaxis()->SetTitle("xF");
      gr[m]->GetYaxis()->SetTitle("SSMD");
      gr[m]->GetYaxis()->SetTitleOffset(1.5);
      gr[m]->SetMarkerColor(ColorModel[m]);  
      gr[m]->SetMarkerStyle(SymbModel[m]);
      gr[m]->SetMarkerSize(1.6);
      if ( secondary == "piplus" || secondary == "piminus" || secondary == "antiproton" )
      {
	 // gr[m]->GetXaxis()->SetRangeUser(-0.55,0.55);
	 gr[m]->GetXaxis()->SetLimits(-0.55,0.55);
      }
      
      // if ( m==0 ) gr1->SetTitle(hi[m]->GetTitle());
   
   }
   
   int iymin = int(ymin) - 1;
   int iymax = int(ymax) + 1;
      
   TLine* line = new TLine();
         
   for ( int m=0; m<NModels; m++ )
   {
      if ( m == 0 )
      {
         // gr[m]->GetYaxis()->SetRangeUser( -2., 1.0 );
         gr[m]->GetYaxis()->SetRangeUser( float(iymin), float(iymax) );
	 /// gr[m]->GetYaxis()->SetLimits( ymin, ymax );
	 gr[m]->Draw("ap");
	 line->DrawLine(-0.55,0.,0.55,0.);
      }
      else
      {
         // gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
         gr[m]->Draw("psame");
      }
   }
  
   TLegend* leg = new TLegend(0.1, 0.70, 0.43, 0.9);   
   leg->SetTextSize(0.025);
   
   for ( int m=0; m<NModels; m++ )
   {
      leg->AddEntry( gr[m], ModelName[m].c_str(), "p" );
   }

   leg->Draw();
   leg->SetFillColor(kWhite);
   
   return;

}

void drawStdDiffHisto( std::string beam, std::string target, std::string secondary, std::string histo )
{

   readIntegratedSpectra( beam, target, secondary );

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
   
   float** StdDiffY = new float*[NModels];
   for ( int m=0; m<NModels; m++ )
   {
      StdDiffY[m] = new float[NPoints];
   }
   
   int np = 0;
   
   TH1F* hst[NModels];

   float ymax = -1.;
   float ymin = 1000000.;
   
   for ( int m=0; m<NModels; m++ )
   {

      std::string histofile = "./na49-histo/" + beam + target + "158.0GeV" + ModelName[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_" + histo ;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );

      int nx = hi[m]->GetNbinsX();
      
      np=0;
      for ( int k=1; k<=nx; k++ ) 
      {
         float xx1 = hi[m]->GetBinLowEdge(k);
	 float xx2 = hi[m]->GetBinWidth(k);
	 for (int kk=0; kk<NPoints; kk++ )
	 {
	    if ( xx1 < xF[kk] && xx1+xx2 > xF[kk] )
	    {	       
	       float yy = hi[m]->GetBinContent(k);
	       float yerror = hi[m]->GetBinError(k);
	       StdDiffY[m][np] = ( yy - Value[kk]) / sqrt( Error[kk]*Error[kk] + yerror*yerror );
	       if ( StdDiffY[m][np] > ymax ) ymax = StdDiffY[m][np];
	       if ( StdDiffY[m][np] < ymin ) ymin = StdDiffY[m][np];
	       np++;
	       break;
	    }
	 }
      } 

   }
   
   float hstxmax = std::max( std::fabs(ymin), std::fabs(ymax) );            
   int ihstxmax = int(hstxmax) + 1;
   ymax = float(ihstxmax);
   ymin = - ymax ;
   int nbins = (2.*float(ihstxmax) / 0.1);

   float hstymax = -1;

   for ( int m=0; m<NModels; m++ )
   {
      
      std::ostringstream osnum;
      osnum << m;
      std::string hstname = "hst" + osnum.str(); 
                   
      hst[m] = new TH1F( hstname.c_str(), hi[m]->GetTitle(), nbins, ymin, ymax );
      std::string xtitle = "SSMD (from " + histo + " spectra)";
      hst[m]->GetXaxis()->SetTitle(xtitle.c_str());
      
      for ( int k=0; k<NPoints; k++ )
      {
         hst[m]->Fill( StdDiffY[m][k] );
      }
      
      // now find out glocal ymax (to set plot's limit later)
      int nx = hst[m]->GetNbinsX();
      for ( int k=1; k<=nx; k++ )
      {
         if ( (hst[m]->GetBinContent(k)+hst[m]->GetBinError(k)) > hstymax ) hstymax = hst[m]->GetBinContent(k)+hst[m]->GetBinError(k);
      }                     
   }
   
   TLegend* leg = new TLegend(0.1, 0.70, 0.45, 0.9);  
   leg->SetTextSize(0.028); 
   
   for ( int m=0; m<NModels; ++m )
   {

      leg->AddEntry( hst[m], ModelName[m].c_str(), "L" );
      
      hst[m]->SetStats(0);
      hst[m]->SetLineColor( ColorModel[m] );
      // hst[m]->SetLineStyle(NModels-m);
      hst[m]->SetLineWidth(2.);       
      if ( m == 0 )
      {
         hst[m]->GetYaxis()->SetRangeUser( 0., hstymax );
	 hst[m]->Draw();
      }
      else
      {
         hst[m]->Draw("same");
      }
   }
   
   leg->Draw();
   leg->SetFillColor(kWhite);
   
   return;

}
