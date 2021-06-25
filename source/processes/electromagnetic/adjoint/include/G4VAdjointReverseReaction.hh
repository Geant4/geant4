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
////////////////////////////////////////////////////////////////////////////////
//  Class:    G4VAdjointReverseReaction
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
//  Abstract class for adjoint/reverse discrete scattering
////////////////////////////////////////////////////////////////////////////////

#ifndef G4VAdjointReverseReaction_h
#define G4VAdjointReverseReaction_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"

class G4AdjointCSManager;
class G4ParticleChange;
class G4ParticleDefinition;
class G4Track;
class G4VEmAdjointModel;
class G4VParticleChange;

class G4VAdjointReverseReaction : public G4VDiscreteProcess
{
 public:
  explicit G4VAdjointReverseReaction(G4String process_name,
                                     G4bool whichScatCase);

  ~G4VAdjointReverseReaction() override;

  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  G4VParticleChange* PostStepDoIt(const G4Track&,
                                  const G4Step&) override;

  G4VAdjointReverseReaction(G4VAdjointReverseReaction&) = delete;
  G4VAdjointReverseReaction& operator=(
    const G4VAdjointReverseReaction& right) = delete;

 protected:
  G4double GetMeanFreePath(const G4Track& track,
                           G4double previousStepSize,
                           G4ForceCondition* condition) override;

  G4VEmAdjointModel* fAdjointModel = nullptr;
  G4bool fIsScatProjToProj;

 private:

  G4ParticleChange* fParticleChange;
  G4AdjointCSManager* fCSManager;

  G4int fTrackId = 0;

  G4bool fIsFwdCSUsed = false;
};

#endif
