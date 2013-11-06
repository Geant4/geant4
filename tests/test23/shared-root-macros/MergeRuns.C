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


void fancyMerge( std::string beam, std::string target, std::string energy, std::string physlist, std::string dir, bool doScale=true )
{
      
   std::string output = beam + target + energy + "GeV" + physlist + ".root" ;
   
   targetFile = TFile::Open( output.c_str(), "RECREATE" );
   
   double scale = 1./32.;
   
   // std::string input = beam + target + model + energy + "GeV-1.root";
   // std::string input = "../t23-bld/harp-histo-no-res-decays/" + beam + target + energy + "GeV" + physlist +"-1.root";

   // std::string input = "../t23-bld/harp-histo/" + beam + target + energy + "GeV" + physlist +"-1.root";

   // std::string input = "../t23-bld/na49-histo/" + beam + target + energy + "GeV" + physlist +"-1.root";
   
   std::string input = dir + "/" + beam + target + energy + "GeV" + physlist + "-1.root";
   
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
	    // std::string input_t = "../t23-bld/na49-histo/" + beam + target + energy + "GeV" + physlist + "-" ;
	    std::string input_t = dir + "/" + beam + target + energy + "GeV" + physlist + "-" ;
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

