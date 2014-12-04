//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

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
