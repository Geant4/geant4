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
                                         G4double previousStepSize,
                                         G4ForceCondition *condition )
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
