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
// Geant4 header G4EmDNAPhysicsActivator
//
// Author V.Ivanchenko 
//
// Build DNA physics on top of standard physics for given region
// 

#ifndef G4EmDNAPhysicsActivator_h
#define G4EmDNAPhysicsActivator_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EmParameters;
class G4ParticleDefinition;
class G4Region;

class G4EmDNAPhysicsActivator : public G4VPhysicsConstructor
{
public:

  explicit G4EmDNAPhysicsActivator(G4int ver = 1);

  ~G4EmDNAPhysicsActivator() override = default;

  void ConstructParticle() override;
  void ConstructProcess() override;

private:

  void DeactivateElectronProcesses(const G4double emaxDNA,
                                   const G4double emax,
                                   const G4Region*);

  void DeactivateHadronProcesses(G4ParticleDefinition*,
                                 const G4double emaxDNA,
                                 const G4double emax,
                                 const G4Region*);

  void DeactivateIonProcesses(G4ParticleDefinition*,
                              const G4double emaxDNA,
                              const G4double emax,
                              const G4Region*);

  void DeactivateNuclearStopping(const G4ParticleDefinition*,
                                 const G4double emaxDNA,
                                 const G4Region*);
 
  G4bool IsVerbose() const;

  G4int verbose;
  G4EmParameters* theParameters;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


