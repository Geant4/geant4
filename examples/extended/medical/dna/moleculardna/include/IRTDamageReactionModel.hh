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
/// \file DetectorMessenger.hh
/// \brief
#ifndef MolecularIRTDamageReactionModel_hh
#define MolecularIRTDamageReactionModel_hh

#include "globals.hh"
#include "G4VDNAHitModel.hh"
#include "G4String.hh"
#include "G4Integrator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"

class DNAGeometry;

class G4DNAComponentNode;

class G4MolecularConfiguration;

class G4DNAMolecularReactionTable;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class IRTDamageReactionModel : public G4VDNAHitModel
{
 public:
  explicit IRTDamageReactionModel(const G4String& name);

  ~IRTDamageReactionModel() override = default;

  IRTDamageReactionModel& operator=(const IRTDamageReactionModel& right) =
    delete;

  IRTDamageReactionModel(const IRTDamageReactionModel&) = delete;

  using DNANode =
    std::variant<const G4DNAComponentNode*, /*for dnadamage chain*/
                 const G4VPhysicalVolume* /*for molecularDNA chain*/>;

  G4double CalculateReactionTime(const G4Track& trackA, DNANode&) override;

  G4bool DoReaction(const G4Track& track, const G4double&,
                    const DNANode&) override;

 private:
  G4double GetTimeToEncounter(const G4MolecularConfiguration* molA,
                              const G4MolecularConfiguration* molB,
                              const G4double& distance);

  void MakeReaction(const G4Track& track);

  void RecordDNADamage() const;

  const G4VPhysicalVolume* fpDNAPhyVolume = nullptr;
  const G4Track* fpTrack                  = nullptr;
  G4double fminTimeStep                   = DBL_MAX;
  G4double fReactionTime                  = DBL_MAX;
  DNAGeometry* fpDNAGeometry;
  const G4DNAMolecularReactionTable* fMolecularReactionTable;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif