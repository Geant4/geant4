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
//
/// file:SteppingAction.hh
/// brief:
#ifndef MOLECULAR_STEPPING_ACTION_HH
#define MOLECULAR_STEPPING_ACTION_HH

#include "G4UserSteppingAction.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"

class EventAction;

class DNAGeometry;

class G4VPhysicalVolume;

class PlacementVolumeInfo;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction : public G4UserSteppingAction
{
 public:
  explicit SteppingAction(EventAction*);

  ~SteppingAction() override;

  void UserSteppingAction(const G4Step*) override;

 protected:
  void ProcessDNARegionHit(const G4Step*);

  void DoChromosomeRegionHit(const G4ThreeVector&, const G4double&,
                             const G4bool&);

  void DoChromosomeDNAHit(const G4ThreeVector&, const G4ThreeVector&,
                          const G4double&, const G4double&, const G4String&,
                          const G4String&);

 private:
  EventAction* fEventAction;
  DNAGeometry* fDNAGeometry;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_STEPPING_ACTION_HH
