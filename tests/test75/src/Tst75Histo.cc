
//#include "TSystem.h"

#include "Tst75Histo.hh"
//#include "TstReader.hh"

#include "Tst75HistoSet.hh"

//#include "G4SystemOfUnits.hh"

//#include <iostream>
//#include <fstream>
//#include <iomanip>


Tst75Histo::Tst75Histo( const TstReader* pset )
   : TstHisto(pset)
{
  
   // NOTE: This methid can be successfully called only AFTER G4ParticleTable is filled up.
   //       In principle, the way things are implemented isn't safe, because nothing guarantee 
   //       that the particle table is alset... 
   //       Will need to go over this aspect some more, and make it safer.
   //
   pset->SyncKinematics();
   
   std::ostringstream ossmom;
   ossmom << pset->GetBeamMomentum();
   
   fBeamMomentum = "";
   fBeamMomentum += ossmom.str();
   fBeamMomentum += "MeV";
     
   fHistoTitle = fBeam + "(" + ossmom.str() + "MeV) + " + fTarget + " -> X "; 
     
   fHistoSet = new Tst75HistoSet( fHistoTitle, pset->GetBeamKineticEnergy() );
   
}

TFile* Tst75Histo::OpenHistoFile()
{
   
   G4String fname = fBeam + fBeamMomentum + fTarget + fModel;
   
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

