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


#ifndef G4VAL_READ_PIM_C
#define G4VAL_READ_PIM_C

// pi- beam business
const int   NPointsMadey = 30;
const int   NTargetsMadey = 7; 
float KENeut[30];  
float Value[7][30], Error[7][30];
std::string TargetsMadey[7] = { "C", "N", "O", "Al", "Cu", "Ta", "Pb" };
 
void readMadey( bool useStatEr=true, bool useSysEr=true )
{

   // NOTE: Overal systematical errors is +/-6.3%
   //       It can be added to the published stat errors if requested (default)
   
   std::string fname = "./piminus/madey_spectra.dat";
   std::cout << "Reads data from file " << fname << "\n";
   ifstream infile;
   infile.open(fname.c_str());
      
   for ( int i=0; i<NPointsMadey; i++ )
   {

      KENeut[i] = 0.;
      infile >> KENeut[i];
      for ( int j=0; j<NTargetsMadey; j++ )
      {
	 Value[j][i] = 0.;
         Error[j][i] = 0.,
         float err = 0.;
	 infile >> Value[j][i] >> err;
	 float err2 = 0.;
	 if ( useStatEr) err2 = err*err;
	 if ( useSysEr ) err2 += (Value[j][i]*0.063)*(Value[j][i]*0.063); 
	 Error[j][i] = sqrt(err2);
	 // rescale as numbers are per 100 pi-'s stopped in the target
	 Value[j][i] /= 100. ;
	 Error[j][i] /= 100. ;
      }
   }
   infile.close();
      
   return;

}

TGraphErrors* getMadeyAsGraph( std::string target )
{

   readMadey();
   int TargetID = findTarget( target );
   if ( TargetID == -1  || TargetID >= NTargetsMadey ) 
   {
      stc::cout << " Invalid Target: " << target << std::endl;      
      return 0;
   }
   TGraphErrors*  gr1 = new TGraphErrors(NPointsMadey,KENeut,Value[TargetID],0,Error[TargetID]);
//   gr1->GetYaxis()->SetRangeUser(0.2,5.);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);

   return gr1;

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

#endif
