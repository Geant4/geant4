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
#ifndef TstHisto_h
#define TstHisto_h

#include <string>
#include "TFile.h"

// #include "G4LorentzVector.hh"
#include "TstReader.hh"

#include "TstHistoSet.hh"

// fwd declaration
//class G4VParticleChange;

class TstHisto
{

public:
      
   // ctor & dtor
   //
   TstHisto( const TstReader* pset );
   virtual ~TstHisto();
      
   void FillEvt( G4VParticleChange* pc, const G4LorentzVector& labv, const G4LorentzVector& labp )
      { fHistoSet->FillEvt( pc, labv, labp ); return; }

   void Write( G4int stat=1, G4double scale=1. );
   void SetJobID( G4int id ) { fJobID=id; return; }
   
   virtual TFile* OpenHistoFile();

protected:
           
  // data members
  //
   G4int fJobID; 
   G4String        fBeam;
   G4String        fTarget;
   G4String        fBeamMomentum;
   G4String        fModel;
   G4String        fHistoTitle;
   TstHistoSet*    fHistoSet;
   TFile*          fHistoFile;

};

#endif
