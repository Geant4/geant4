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
/// \file include/SpecializedTrackingManager.hh
/// \brief Definition of the SpecializedTrackingManager class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SpecializedTrackingManager_h
#define SpecializedTrackingManager_h 1

#include "G4VTrackingManager.hh"

class G4Region;

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SpecializedTrackingManager : public G4VTrackingManager
{
 public:
  SpecializedTrackingManager();
  ~SpecializedTrackingManager();

  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  void PreparePhysicsTable(const G4ParticleDefinition&) override;

  void HandOverOneTrack(G4Track* aTrack) override;
  void FlushEvent() override;

 private:
  void StepInBackRegion(G4Track* aTrack);
  void StepOutside(G4Track* aTrack);

  std::vector<G4Track*> fBufferedTracks;
  G4Region* fBackRegion = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
