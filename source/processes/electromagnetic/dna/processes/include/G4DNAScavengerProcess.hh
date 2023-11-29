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

#ifndef G4DNASCAVENGERPROCESS_HH
#define G4DNASCAVENGERPROCESS_HH

#include "G4VITProcess.hh"
class G4DNAMolecularReactionData;
class G4MolecularConfiguration;
class G4DNABoundingBox;
class G4DNAScavengerMaterial;
class G4DNAScavengerProcess : public G4VITProcess
{
 public:
  using MolType = const G4MolecularConfiguration*;
  using Data    = const G4DNAMolecularReactionData;
  explicit G4DNAScavengerProcess(const G4String& aName,
                                 const G4DNABoundingBox& box,
                                 G4ProcessType type = fUserDefined);
  ~G4DNAScavengerProcess() override;
  G4DNAScavengerProcess(const G4DNAScavengerProcess&) = delete;
  G4DNAScavengerProcess& operator=(const G4DNAScavengerProcess&) = delete;
  void StartTracking(G4Track*) override;
  void SetReaction(MolType, Data* pData);

 public:
  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  G4double PostStepGetPhysicalInteractionLength(
    const G4Track& track, G4double previousStepSize,
    G4ForceCondition* condition) override;

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

  G4double AtRestGetPhysicalInteractionLength(const G4Track&,
                                              G4ForceCondition*) override
  {
    return -1.0;
  }

  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) override
  {
    return nullptr;
  }

  //  no operation in  AlongStepDoIt
  G4double AlongStepGetPhysicalInteractionLength(const G4Track&, G4double,
                                                 G4double, G4double&,
                                                 G4GPILSelection*) override
  {
    return -1.0;
  }

  //  no operation in  AlongStepDoIt
  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) override
  {
    return nullptr;
  }

 protected:
  struct G4DNAScavengerProcessState : public G4ProcessState
  {
    G4DNAScavengerProcessState();
    ~G4DNAScavengerProcessState() override { ; }
    G4double fPreviousTimeAtPreStepPoint;
    G4bool fIsInGoodMaterial;
  };

 protected:
  G4bool fIsInitialized;
  G4double fReturnedValue;
  G4ParticleChange fParticleChange;
  const G4MolecularConfiguration* fpMolecularConfiguration;
  std::map<MolType /*MolConf*/, std::map<MolType /*molConfMat*/, Data*>>
    fConfMap;
  std::vector<MolType> fpMaterialVector;
  MolType fpMaterialConf;
  const G4DNABoundingBox* fpBoundingBox;
  G4DNAScavengerMaterial* fpScavengerMaterial;
};
#endif  // FLASH1_G4DNASCAVENGERPROCESS_HH
