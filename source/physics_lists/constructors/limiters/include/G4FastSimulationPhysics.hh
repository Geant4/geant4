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

#ifndef G4FastSimulationPhysics_h
#define G4FastSimulationPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4FastSimulationPhysics : public G4VPhysicsConstructor
{
public:
  
  G4FastSimulationPhysics(const G4String& name = "FastSimP");
  virtual ~G4FastSimulationPhysics();

public:
  // --------------------------------------
  // -- Fast simulation activation methods:
  // --------------------------------------
  // -- Used to select particles for which fast simulation has to be activated:
  // ---- Activate fast simulation for one particle with fast simulation attached to mass geometry:
  void ActivateFastSimulation(const G4String& particleName);
  // ---- Activate fast simulation for one particle with fast simulation attached to a parallel geometry:
  void ActivateFastSimulation(const G4String& particleName, const G4String& parallelGeometryName);

  // -- Information about particles under fast simulation:
  void BeVerbose() { fVerbose = true; }
  
public:
  
  // This method is dummy for physics
  virtual void ConstructParticle();
  
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();
  
private:
  
  // hide assignment operator
  G4FastSimulationPhysics & operator=(const G4FastSimulationPhysics &right);
  G4FastSimulationPhysics(const G4FastSimulationPhysics&);

  // -- Particles under fast simulation:
  std::vector< G4String >  fParticlesUnderFastSimulation;
  // -- And their related possible parallel geometries:
  std::vector< G4String >  fGeometries;

  // -- Report:
  G4bool fVerbose;


  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
