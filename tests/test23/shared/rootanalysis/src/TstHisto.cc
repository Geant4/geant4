
#include "TstHisto.hh"

#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <iomanip>

TstHisto::TstHisto( const TstReader* pset )
    : fJobID(pset->GetJobID()), 
      fBeam(pset->GetBeamParticle()), 
      fTarget(pset->GetTargetMaterial()), 
      fModel(pset->GetPhysics()),
      fHistoSet(0), 
      fHistoFile(0)
{
   
   char ene[3];
   sprintf( ene, "%3.1f", pset->GetBeamMomentum()/GeV );
   fBeamMomentum = "";
   fBeamMomentum.append( ene );
   fBeamMomentum += "GeV";
   
   fHistoTitle = fBeam + " + " + fTarget;

}

TstHisto::~TstHisto()
{

   if ( fHistoSet )  delete fHistoSet;
   if ( fHistoFile ) delete fHistoFile;

}

void TstHisto::Write( G4int stat, G4double scale )
{

/*
   std::string fname = fBeam + fTarget + fModel;
   if ( fJobID > -1 )
   {
      char buf[5];
      sprintf( buf, "%i", fJobID );
      fname += "-";
      fname.append( buf ); 
   }  
   fname += ".root";

   // std::cout << "Writing histogram file: " << fname << std::endl;
*/

   if ( fHistoFile ) delete fHistoFile;

   fHistoFile = OpenHistoFile();

   fHistoFile->cd();
   // WriteHisto( stat );
   fHistoSet->Write( stat, scale );

   fHistoFile->Close();
   
   return;

}

TFile* TstHisto::OpenHistoFile()
{

   G4String fname = fBeam + fTarget + fBeamMomentum + fModel;
   
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
