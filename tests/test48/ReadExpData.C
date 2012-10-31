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


// pi- beam business
const int   NPointsMadey = 30;
const int   NTargetsMadey = 7; 
float KENeut[30];  
float Value[7][30], Error[7][30];
std::string TargetsMadey[7] = { "C", "N", "O", "Al", "Cu", "Ta", "Pb" };
 
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
