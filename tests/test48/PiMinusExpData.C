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

static int isExpReadLoaded = 0;


void PiMinusExpData( std::string target )
{

   if ( isExpReadLoaded <= 0 )
   {   
      gROOT->LoadMacro("ReadExpData.C");
      isExpReadLoaded = 1;
   }
   
   readMadey();
   std::cout << "NPMadey= " << NPointsMadey << " KENeut= " << KENeut[0] << std::endl;
   
   std::string fname = "expdata_piminus" + target + ".root";
   TFile* fout = new TFile( fname.c_str(), "RECREATE");


   int TargetID = findTarget( target );
   
   std::cout << " TargetID= " << TargetID << std::endl;
   
   TGraph*  gr = new TGraphErrors(NPointsMadey,KENeut,Value[TargetID],0,Error[TargetID]);
   gr->SetName("ExpDataGraph");
   std::string title = "piminus on " + target + ", exp.data";
   gr->SetTitle(title.c_str());
   gr->SetMarkerColor(4);  gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.6);
   gr->GetYaxis()->SetRangeUser( 0., 1. );
//    gr->Draw("apl");

   // NOTE: this histogram will be filled up with the exp numbers 
   //       but without the exp. measurement errors !!!
   //
   TH1F* hist = new TH1F("ExpDataHisto", title.c_str(), 150, 0., 150. );
   for ( int i=0; i<NPointsMadey; i++ )
   {
      hist->Fill( KENeut[i], Value[TargetID][i] );
   }
      
   gr->Write();
   hist->Write();
   
   fout->Close();
   
   return;

}
