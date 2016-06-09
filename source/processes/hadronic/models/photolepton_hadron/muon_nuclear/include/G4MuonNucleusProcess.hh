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
// G4MuonNucleusProcess.hh
//
//     M.Takahata (Makoto.Takahata@cern.ch)

#ifndef G4MuonNucleusProcess_h
#define G4MuonNucleusProcess_h 1

#include "globals.hh"
#include "G4LeptonHadronProcess.hh"
#include "G4MuonNucleusInteractionModel.hh"


  class G4MuonNucleusProcess : public G4LeptonHadronProcess
  {
    public:

      G4MuonNucleusProcess(const G4String& processName ="MuonNucleus");
      virtual ~G4MuonNucleusProcess();

      G4double GetMeanFreePath(const G4Track &muonTrack,
                               G4double previousStepSize,
                               G4ForceCondition *condition);

    private:
      G4MuonNucleusProcess(const G4MuonNucleusProcess &right);

  };

#endif
