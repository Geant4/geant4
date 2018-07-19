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
// $Id: G4StepLimiterPhysics.hh 107562 2017-11-22 15:39:24Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4StepLimiterPhysics_h
#define G4StepLimiterPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4StepLimiter;
class G4UserSpecialCuts;

class G4StepLimiterPhysics : public G4VPhysicsConstructor
{
public:

  G4StepLimiterPhysics(const G4String& name = "stepLimiter");
  virtual ~G4StepLimiterPhysics();

public:

  // This method is dummy for physics
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();

  // Option to apply the step limit to all particles
  // (by default the step limit is applied to charged particles only)
  void   SetApplyToAll(G4bool option) { fApplyToAll= option; }
  G4bool GetApplyToAll() const { return fApplyToAll; }

private:

   // hide assignment operator
  G4StepLimiterPhysics & operator=(const G4StepLimiterPhysics &right);
  G4StepLimiterPhysics(const G4StepLimiterPhysics&);

  G4bool  fApplyToAll;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
