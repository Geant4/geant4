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

//const int   NModels = 4;
const int   NModels = 3;
std::string ModelName[4]  = { "NuBeam", "qgsp_bert", "ftfp_bert", "qgsp_ftfp_bert", };
int         ColorModel[4] = { kMagenta, 7, kRed, kBlack }; // 14 = grey, 7 = light "sky"-blue

// read exp.data
//
// FW data (forward = small angle with respect to projectile)
//
const int NSetsFW = 4;
//
// LA data (lateral = larger angle with respect to projectile)
//
const int NSetsLA = 9;

// general purpose counter
//
static int NSets = 0;

float** AngleBin = 0;
int*    NPoints = 0; 
float** XMin = 0;
float** XMax = 0;
float** Y = 0;
float** EY = 0;

void ReadHARPData( std::string beam, std::string target, std::string energy, 
                   std::string secondary, 
		   std::string region ) // should be either "FW" or "LA"
{

   std::string dirname = "./harp-exp-data/";
   
   std::string filename = beam + "_" + target + "_" + energy + "GeV_" + secondary + "_" + region + ".dat";
   
   std::string file = dirname + filename;
   
   ifstream infile;
   infile.open( file.c_str() );
      
   if (AngleBin) 
   {
      delete [] AngleBin[0];
      delete [] AngleBin[1];
      delete [] AngleBin;
   }    
   
   if ( NPoints ) delete [] NPoints;   
   
   int i = 0;
   
   if ( NSets > 0 )
   {
      if ( XMin )
      {
         for ( i=0; i<NSets; i++ )
         {
            delete [] XMin[i];
         }
         delete [] XMin;
      }   
      if ( XMax )
      {
         for ( i=0; i<NSets; i++ )
         {
            delete [] XMax[i];
         }
         delete [] XMax;
      }   
      if ( Y )
      {
         for ( i=0; i<NSets; i++ )
         {
            delete [] Y[i];
         }
         delete [] Y;
      }      
      if ( EY )
      {
         for ( i=0; i<NSets; i++ )
         {
            delete [] EY[i];
         }
         delete [] EY;
      }
      NSets = 0;
   }   

   if ( region == "FW" )
   {
      NSets = NSetsFW;
   }
   else if ( region == "LA" )
   {
      NSets = NSetsLA;
   }
   else
   {
      return;
   }
   
   AngleBin = new float*[2];
   AngleBin[0] = new float[NSets];
   AngleBin[1] = new float[NSets]; 

   NPoints = new int[NSets];

   XMin = new float*[NSets];
   XMax = new float*[NSets];
   Y    = new float*[NSets];
   EY   = new float*[NSets]; 

   for ( i=0; i<NSets; i++ )
   {
      infile >> AngleBin[0][i] >> AngleBin[1][i];
      infile >> NPoints[i];
      
      // std::cout << "Angle Bin: " << AngleBin[0][i] << " " << AngleBin[1][i] << std::endl;
      // std::cout << "NPoints: " << NPoints[i] << std::endl;
      
      XMin[i] = new float[NPoints[i]];
      XMax[i] = new float[NPoints[i]];
      Y[i]    = new float[NPoints[i]];
      EY[i]   = new float[NPoints[i]];
      for ( int j=0; j<NPoints[i]; j++ )
      {
         infile >> XMin[i][j] >> XMax[i][j] >> Y[i][j] >> EY[i][j];
	 Y[i][j] *= 1000.;
	 EY[i][j] *= 1000.;
      }
   }
   
   infile.close();

   return;

}

void PlotHARPAnalysis( std::string beam, std::string target, std::string energy,
                       std::string secondary,
		       std::string region )
{

   ReadHARPData( beam, target, energy, secondary, region );
   
   TCanvas* myc = 0;
   
   if ( region == "FW" )
   {   
      myc = new TCanvas( "myc", "", 800, 800 );
      myc->Divide( 2, 2 );
   }
   else if ( region == "LA" )
   {
      myc = new TCanvas( "myc", "", 1200, 1200 );
      myc->Divide( 3, 3 );
   }

   for ( int i=0; i<NSets; i++ )
   {
      myc->cd(i+1);
      PlotHARPHisto( beam, target, energy, secondary, region, i );
   }
   
   myc->cd();
   
   return;

}
void PlotHARPHisto( std::string beam, std::string target, std::string energy,
                    std::string secondary,
		    std::string region,
		    int ibin )
{
   
   // ReadHARPData( beam, target, energy, secondary, region );
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
   for ( int i=0; i<NPoints[ibin]; i++ )
   {
      if ( (Y[ibin][i]+EY[ibin][i]) > ymax ) ymax = Y[ibin][i]+EY[ibin][i];
      if ( (Y[ibin][i]-EY[ibin][i]) < ymin && (Y[ibin][i]-EY[ibin][i]) > 0. ) ymin = (Y[ibin][i]-EY[ibin][i]);
   }
   
   TH1F* hi[NModels];
   std::string YTitle;

   for ( int m=0; m<NModels; m++ )
   {

      std::string histofile = "";
      
      // histofile = "./harp-histo/";
      histofile = "./harp-histo-no-res-decays/";
      // histofile = "../t23-bld/harp-histo/";
      
      // std::string histofile = "./harp-histo/" + beam + target + energy + "GeV" + ModelName[m] + ".root"; 
      histofile += ( beam + target + energy + "GeV" + ModelName[m] + ".root" ); 
      
      // std::cout << " histofile = " << histofile << std::endl;
      
      TFile* f = new TFile( histofile.c_str() );
      
      char buf[5];
      sprintf( buf, "%i", ibin );      
      std::string histoname = secondary + "_" + region + "_";
      histoname.append( buf );
      
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
	 hi[m]->GetXaxis()->SetTitle("momentum (GeV/c)");
	 //hi[m]->GetYaxis()->SetTitle( YTitle.c_str() );
	 // hi[m]->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dpd#Theta} [mb/(GeV/c/rad)]");
	 hi[m]->GetYaxis()->SetTitle("d^{2}#sigma / dpd#Theta [mb/(GeV/c/rad)]");
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
      
   float* X = new float[NPoints[ibin]];
   for ( int i=0; i<NPoints[ibin]; i++ )
   {
      X[i] = 0.5 * (XMin[ibin][i]+XMax[ibin][i]);
      //std::cout << "X[" << i << "] = " << X[i] << std::endl;
      //std::cout << "Y[" << i << "] = " << Y[0][i] << std::endl;
   }
   
   TGraph* gr = new TGraphErrors( NPoints[ibin], X, Y[ibin], 0, EY[ibin] );
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
    
   gr->Draw("p");
      
   leg->AddEntry( gr, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);   


   return;

}

void fancyMerge( std::string beam, std::string target, std::string energy, std::string physlist, bool doScale=false )
{
      
   std::string output = beam + target + energy + "GeV" + physlist + ".root" ;
   
   targetFile = TFile::Open( output.c_str(), "RECREATE" );
   
   double scale = 1./32.;
   
   // std::string input = beam + target + model + energy + "GeV-1.root";
   // std::string input = "../t23-bld/harp-histo-no-res-decays/" + beam + target + energy + "GeV" + physlist +"-1.root";

   // std::string input = "../t23-bld/harp-histo/" + beam + target + energy + "GeV" + physlist +"-1.root";

   std::string input = "../t23-bld/na49-histo/" + beam + target + energy + "GeV" + physlist +"-1.root";
   
   TFile* iFile1 = TFile::Open( input.c_str() );
   TIter  next( iFile1->GetListOfKeys() );
   TKey*  key = (TKey*)next();
   TH1* h = 0;
   while ( key )
   {   
         if ( !(TClass::GetClass(key->GetClassName())->InheritsFrom(TH1::Class())) ) continue;
         const char* kName = key->GetName();
         h = (TH1*)key->ReadObj();
         const char* hName = h->GetName();
         std::cout << " histoname = " << hName << std::endl;
	 TH1F* h1 = h->Clone();
	 for ( int id=2; id<=32; id++ )
	 {
	    // std::string input_t = "../t23-bld/harp-histo-no-res-decays/" + beam + target + energy + "GeV" + physlist + "-" ;
	    // std::string input_t = "../t23-bld/harp-histo/" + beam + target + energy + "GeV" + physlist + "-" ;
	    std::string input_t = "../t23-bld/na49-histo/" + beam + target + energy + "GeV" + physlist + "-" ;
            char buf[5];
            sprintf( buf, "%i", id );
            input_t.append( buf ); 
            input_t += ".root"; 
	    TFile* iFile_t = TFile::Open( input_t.c_str() );
	    TH1F* h_t = (TH1F*)iFile_t->Get( h->GetName() );
	    h1->Add( h_t );  
	    iFile_t->Close();
	 }
	 if ( doScale )
	 {
	    if (!(strcmp(key->GetClassName(),"TProfile"))) h1->Scale( scale );
	 }
	 targetFile->cd();
	 h1->Write();
         key = (TKey*)next();
   }
   
   targetFile->Close();
     
   return;

}

