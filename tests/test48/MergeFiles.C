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

#include "TFileMerger.h"



void crudeMerge( std::string beam, std::string target, std::string model )
{


//Note: if one merges weighted/normilized histograms, in this case 
//      the resulting histogram(s) will be 32 times the statistics;
//      at the analysis stage they'd need to be scaled by 1/32 

   TFileMerger* fm = new TFileMerger();
   
   std::string output = beam + target + model + ".root";
   
   fm->OutputFile( output.c_str() );
   
   for ( int id=1; id<=32; id++ )
   {
      std::string filename = beam + target + model + "-";
      char buf[5];
      sprintf( buf, "%i", id );
      filename.append( buf ); 
      filename += ".root";    
      // std::cout << " file name = " << file << std::endl;           
      fm->AddFile( filename.c_str() );
   }
   
   fm->Merge();
   
   return;

}

void fancyMerge( std::string beam, std::string target, std::string model, bool doScale=false )
{
      
   std::string output = beam + target + model + ".root" ;
   
   targetFile = TFile::Open( output.c_str(), "RECREATE" );
   
   double scale = 1./32.;
   
   std::string input = beam + target + model + "-1.root";
   
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
	    std::string input_t = beam + target + model + "-" ;
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
	    h1->Scale( scale );
	 }
	 targetFile->cd();
	 h1->Write();
         key = (TKey*)next();
   }
   
   targetFile->Close();
     
   return;

}


