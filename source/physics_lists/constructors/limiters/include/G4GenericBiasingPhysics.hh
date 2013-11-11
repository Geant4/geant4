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
// $Id: $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4GenericBiasingPhysics_h
#define G4GenericBiasingPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4GenericBiasingPhysics : public G4VPhysicsConstructor
{
public:
  
  G4GenericBiasingPhysics(const G4String& name = "BiasingP");
  virtual ~G4GenericBiasingPhysics();

public:
  // -- Used to select particles and processes to be under biasing:
  // ---- Put under biasing all physics processes of given particleName:
  void PhysicsBias(const G4String& particleName);
  // ---- Put under biasing processes in processToBiasNames of given particleName:
  void PhysicsBias(const G4String& particleName, const std::vector< G4String >& processToBiasNames);
  // ---- Allow for non physics biasing for particle:
  void NonPhysicsBias(const G4String& particleName);
  // ---- Put under biasing all physics processes and allow for non physics biasing:
  void Bias(const G4String& particleName);
  // ---- Put under biasing processes in processToBiasNames of given particleName:
  void Bias(const G4String& particleName, const std::vector< G4String >& processToBiasNames);
  
  
public:
  
  // This method is dummy for physics
  virtual void ConstructParticle();
  
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();
  
private:
  
  // hide assignment operator
  G4GenericBiasingPhysics & operator=(const G4GenericBiasingPhysics &right);
  G4GenericBiasingPhysics(const G4GenericBiasingPhysics&);

  // -- Particles under biasing:
  std::vector< G4String >  fBiasedParticles;
  std::vector< G4bool >   fBiasAllProcesses;
  // -- Related biased processes:
  std::vector< std::vector< G4String > > fBiasedProcesses;
  // -- non physics biased particles:
  std::vector< G4String > fNonPhysBiasedParticles;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
