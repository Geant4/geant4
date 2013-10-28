
#include "TSystem.h"

#include "Tst23Histo.hh"
#include "TstReader.hh"

#include "Tst23NA49Histo.hh"
#include "Tst23HARPHisto.hh"

#include <iostream>
#include <fstream>
#include <iomanip>


Tst23Histo::Tst23Histo( const TstReader* pset )
   : TstHisto(pset)
{
   
   if ( pset->GetExpDataSet() == "NA49" )
   {
      fHistoSet = new Tst23NA49Histo( fHistoTitle );
      fHistoDirName  = "na49-histo";
   }
//   else if ( pset->GetExpDataSet() == "NA61" )
//   {
//      fHistoSet = new TestNA49Histo( fHistoTitle );
//      fHistoDirName = "na61-histo";
//   }
   else if ( pset->GetExpDataSet() == "HARP" )
   {
      fHistoSet = new Tst23HARPHisto( fHistoTitle );
      fHistoDirName = "harp-histo";
   }

}

TFile* Tst23Histo::OpenHistoFile()
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
