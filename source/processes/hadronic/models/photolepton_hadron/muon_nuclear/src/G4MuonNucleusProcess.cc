//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
    theInteractionModel = chooseInteractionModel();
  }


//-----------------------------------------------------------------------------
  G4MuonNucleusProcess::~G4MuonNucleusProcess()
//-----------------------------------------------------------------------------
  {
    delete theInteractionModel;
  }


//-----------------------------------------------------------------------------
  G4LeptonHadronInteractionModel*
  G4MuonNucleusProcess::chooseInteractionModel()
//-----------------------------------------------------------------------------
  {
    G4MuonNucleusInteractionModel* aModel = new G4MuonNucleusInteractionModel;
    return aModel;
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
