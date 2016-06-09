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
// G4MuonNucleusProcess.cc
//
//     M.Takahata (Makoto.Takahata@cern.ch)

#include "G4MuonNucleusProcess.hh"


//-----------------------------------------------------------------------------
  G4MuonNucleusProcess::G4MuonNucleusProcess(const G4String& processName)
//-----------------------------------------------------------------------------
    : G4LeptonHadronProcess( processName)
  {
    theInteractionModel = new G4MuonNucleusInteractionModel;
  }


//-----------------------------------------------------------------------------
  G4MuonNucleusProcess::~G4MuonNucleusProcess()
//-----------------------------------------------------------------------------
  {
    delete theInteractionModel;
  }

//-----------------------------------------------------------------------------
  G4double
  G4MuonNucleusProcess::GetMeanFreePath( const G4Track &muonTrack,
                                         G4double,
                                         G4ForceCondition * )
//-----------------------------------------------------------------------------
  {
    G4Material *aMaterial  = muonTrack.GetMaterial();

    const G4double*  theAtomicNumDensityVector 
      = aMaterial->GetAtomicNumDensityVector();
    const G4int  theNumberOfElements
      = aMaterial->GetNumberOfElements();

    G4double macroscopicCrossSection = 0.0;
    for(G4int iel=0; iel<theNumberOfElements; iel++)
    {
       macroscopicCrossSection 
         += theAtomicNumDensityVector[iel]
           *theInteractionModel->computeMicroscopicCrossSection(muonTrack);
    }

    if( macroscopicCrossSection > 0.0 ) {
      return 1.0/macroscopicCrossSection;
    } else {
      return DBL_MAX;
    }

  }
