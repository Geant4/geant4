
#include "TSystem.h"

#include "Test19Histo.hh"
#include "TstReader.hh"

#include "TestNA49Histo.hh"
#include "TestNA61Histo.hh"

#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <iomanip>


Test19Histo::Test19Histo( const TstReader* pset )
   : TstHisto(pset)
{

   if ( pset->GetExpDataSet() == "NA61" )
   {
      fHistoSet = new TestNA61Histo( fHistoTitle );
      fHistoDirName  = "na61-histo";
   }
   else if ( pset->GetExpDataSet() == "NA49" )
   {
      fHistoSet = new TestNA49Histo( fHistoTitle );
      fHistoDirName  = "na49-histo";
   }

}

TFile* Test19Histo::OpenHistoFile()
{

   // AccessPathName actually returns o (false) is directory ** exists **
   // and 1 (true) if it does NOT
   //
   if ( gSystem->AccessPathName( fHistoDirName.c_str() ) )
   {   
      gSystem->mkdir( fHistoDirName.c_str() );
   }
   
   G4String fname = fHistoDirName + "/" + fBeam + fTarget + fBeamMomentum + fModel;
   
   if ( fJobID > -1 )
   {
      char buf[5];
      sprintf( buf, "%i", fJobID );
      fname += "-";
      fname.append( buf ); 
   }  
   fname += ".root";

   return new TFile( fname.c_str(), "recreate" );

}
