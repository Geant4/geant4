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
// -------------------------------------------------------------------
//
// Geant4 Class header file
//
// File name:     G4PolarizedAnnihilation
//
// Author:   A. Schaelicke on base of Vladimir Ivanchenko / Michel Maire code
//
// Class Description:
//  Polarized process of e+ annihilation into 2 gammas
//
// -------------------------------------------------------------------

#ifndef G4PolarizedAnnihilation_h
#define G4PolarizedAnnihilation_h 1

#include "globals.hh"
#include "G4eplusAnnihilation.hh"

class G4PolarizedAnnihilationModel;

class G4PolarizedAnnihilation : public G4eplusAnnihilation
{
 public:
  explicit G4PolarizedAnnihilation(const G4String& name = "pol-annihil");

  virtual ~G4PolarizedAnnihilation() override;

  virtual void ProcessDescription(std::ostream&) const override;
  virtual void DumpInfo() const override { ProcessDescription(G4cout); };

  virtual G4double GetMeanFreePath(const G4Track& track,
                                   G4double previousStepSize,
                                   G4ForceCondition* condition) override;

  virtual G4double PostStepGetPhysicalInteractionLength(
    const G4Track& track, G4double previousStepSize,
    G4ForceCondition* condition) override;

  virtual void BuildPhysicsTable(const G4ParticleDefinition&) override;

  G4PolarizedAnnihilation& operator=(const G4PolarizedAnnihilation& right) =
    delete;
  G4PolarizedAnnihilation(const G4PolarizedAnnihilation&) = delete;

 private:
  void CleanTables();

  void BuildAsymmetryTables(const G4ParticleDefinition& part);

  G4double ComputeAsymmetry(G4double energy, const G4MaterialCutsCouple* couple,
                            const G4ParticleDefinition& particle, G4double cut,
                            G4double& tasm);

  G4double ComputeSaturationFactor(const G4Track& aTrack);

  G4PolarizedAnnihilationModel* fEmModel;

  // table for cross section asymmetry
  G4PhysicsTable* fAsymmetryTable;
  // table for transverse cross section asymmetry
  G4PhysicsTable* fTransverseAsymmetryTable;
};

#endif
