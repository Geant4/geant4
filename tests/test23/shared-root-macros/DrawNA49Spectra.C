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
                             std::string secondary, std::string histo,
			     bool NuERange=false )
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
      YTitle = "Average pT, GeV/c";
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
         if ( histo == "dNdxF" )
	 {
	    hi[m]->Draw("hist");
	 }
	 else
	 {
	    hi[m]->Draw("e");
	 }
	 // hi[m]->Draw();
	 hi[m]->GetXaxis()->SetTitle("xF");
	 hi[m]->GetXaxis()->SetTitleSize(0.05);
	 hi[m]->GetXaxis()->SetTitleOffset(0.9);
	 hi[m]->GetXaxis()->CenterTitle();
	 hi[m]->GetYaxis()->SetTitle( YTitle.c_str() );
//	 hi[m]->GetYaxis()->SetTitleOffset(1.5);
	 hi[m]->GetYaxis()->SetTitleSize(0.05);
	 hi[m]->GetYaxis()->CenterTitle();
      }
      else
      {
         if ( histo == "dNdxF" )
	 {
	    hi[m]->Draw("histsame");
	 }
	 else
	 {
	    hi[m]->Draw("esame");
	 } 
      }   

      if ( NuERange )
      {
         for ( int ir=0; ir<4; ++ir )
	 {
	    std::ostringstream osCount;
            // osCount.clear();
            osCount << ir;
            std::string hname_tmp = histoname + "_" + osCount.str();
	    TH1F* htmp = (TH1F*)f->Get( hname_tmp.c_str() );
	    htmp->SetLineColor(kBlack);
	    htmp->SetLineWidth(2);
	    htmp->SetLineStyle(ir);
	    htmp->Draw("same");
         }
      }
   
   }
   
   TLegend* leg = new TLegend(0.67, 0.74, 0.97, 0.99);
   leg->SetTextSize(0.035);
   
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
	       double yy  = hi[m]->GetBinContent(k);
	       double eyy = hi[m]->GetBinError(k);
	       MC2DataX[np] = xF[kk];
	       DX[np] = 0.;
	       
	       MC2DataY[np] = yy / Value[kk];
	       // also need error calc here !...
	       // DY[np]=0.;
	       DY[np] = Error[kk]*MC2DataY[np] / Value[kk];
	       //
	       // the right way to do it is this:
	       // enew2 = (e0**2 * c1**2 + e1**2 * c0**2) / c1**4
	       // enew = sqrt( enew2 )
	       //
	       // but in our case of 1M beam events (~330K inelastic interactions)
	       // the sim error will be about 5% vs 0.5-1.5% in the exp.data, 
	       // thus the sim error will largely dominate... 
	       // need to get at least 10 times more the statistics !...
	       //
	       // double etmp2 = ( eyy*eyy + Error[kk]*Error[kk]*MC2DataY[np]*MC2DataY[np] );
	       // DY[np] = sqrt(etmp2) / Value[kk];
	       //
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
//      gr[m]->GetXaxis()->SetTitle("xF");
//      gr[m]->GetYaxis()->SetTitle("MC/Data");
      gr[m]->SetMarkerColor(ColorModel[m]);  
      gr[m]->SetMarkerStyle(SymbModel[m]);
      gr[m]->SetMarkerSize(1.2);
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
   xaxis->SetTitleSize(0.05);
   xaxis->SetTitleOffset(0.9);
   xaxis->CenterTitle();
   
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
      gr1->GetYaxis()->SetTitle("MC/Data (Average pT, GeV/c)" );
   }
   gr1->GetYaxis()->SetTitleSize(0.05);
//   gr1->GetYaxis()->SetTitleOffset(1.5);
   gr1->GetYaxis()->CenterTitle();



   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->Draw("apl");
   // gr1->Draw("AC*");    
   
   for ( int m=0; m<NModels; m++ )
   {
      // gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[m]->Draw("lpsame");
   }
  
   TLegend* leg = new TLegend(0.1, 0.15, 0.4, 0.35);   
   
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
         hi[iv]->Draw("hist");
	 hi[iv]->GetXaxis()->SetTitle("xF");
	 hi[iv]->GetXaxis()->SetTitleSize(0.05);
	 hi[iv]->GetXaxis()->SetTitleOffset(0.9);
	 hi[iv]->GetXaxis()->CenterTitle();
	 hi[iv]->GetYaxis()->SetTitle( YTitle.c_str() );
	 hi[iv]->GetYaxis()->SetTitleOffset(1.5);
      }
      else hi[iv]->Draw("samehist");     

   }

   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
   leg->SetTextSize(0.035);
   
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

// this routine will draw a ** single ** plot for a selected xF-bin 
// (specified by icount, from 0 to 21)
//
void draw1DDiffXSec( std::string beam, std::string target, std::string secondary, int icount )
{

   TLegend* leg = new TLegend(0.55, 0.7, 0.99, 0.9);
   leg->SetTextSize(0.035);
   
   for ( int m=0; m<NModels; ++m )
   {

      std::string histofile = "./na49-histo/" + beam + target + "158.0GeV" + ModelName[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );

      std::ostringstream cnt;
      cnt << (icount+7);
      std::string hname = "pT";
      if ( secondary == "piplus" )
      {
	    hname +="pip";
      }
      else if ( secondary == "piminus" )
      {
	    hname += "pim";
      }
      hname += cnt.str();
      TH1F* h = (TH1F*)f->Get( hname.c_str() );
      h->Scale(226.);
      // h->Scale( 251. ); // is it G4's xsec p+c at 158GeV ???
      h->SetStats(0);
      h->SetLineColor(ColorModel[m]);
      h->SetLineWidth(2);
      h->SetTitle(0);
      h->GetYaxis()->SetRangeUser( 0.001, 700. );
      leg->AddEntry( h, ModelName[m].c_str(), "L" );
      if ( m == 0 ) 
      {
	    gPad->SetLogy();
	    std::string newname = "p+C ->" + secondary + " + X at 158GeV/c, xF=";
            std::ostringstream xf;
	    xf << DDiff_xF[icount];;
	    newname += xf.str();
	    h->SetTitle( newname.c_str() );
	    h->GetXaxis()->SetTitleOffset(1.2);
	    h->GetXaxis()->SetTitle( "p_{T}, GeV/c" );
	    h->GetYaxis()->SetTitleOffset(1.5);
	    h->GetYaxis()->SetTitle( "f(x_{F},p_{T}), mb/(GeV^{2}/c^{3})" );
	    h->Draw();
      }
      else
      { 
	    h->Draw("same");
      }
   } 
   
   int NPt = DDiffDataHolder[icount].GetEntries();
   double*    tmpPT    = new double[NPt];
   double*    tmpXSec  = new double[NPt];
   double*    tmpEXSec = new double[NPt];
   for ( int i=0; i<NPt; ++i )
   {
      
      TVector3* tmpVec = DDiffDataHolder[icount].At(i);
      tmpPT[i]    = tmpVec->x();
      tmpXSec[i]  = tmpVec->y();
      tmpEXSec[i] = tmpVec->z();
   }
   
   TGraph* gr = new TGraphErrors( NPt, tmpPT, tmpXSec, 0, tmpEXSec );
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(22);
   gr->Draw("psame");
   gr->SetMarkerSize(1.0);
    
   gr->Draw("p");
      
   leg->AddEntry( gr, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);  
   
   delete [] tmpPT;
   delete [] tmpXSec;
   delete [] tmpEXSec; 
               
   return;

}
