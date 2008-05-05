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
// -------------------------------------------------------------
//      GEANT 4 class
//
//      ---------- Test30Material -------
//                by Vladimir Ivanchenko, 12 March 2002
//
//    Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Test30VSecondaryGenerator_Test30VSecondaryGenerator_h
#define Test30VSecondaryGenerator_Test30VSecondaryGenerator_h 1

#include <string>
#include <vector>
#include "globals.hh"
#include "G4LorentzVector.hh"
#include "G4Nucleus.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Track;
class G4ParticleDefinition;
class G4HadronicInteraction;
class G4Material;
class G4VParticleChange;
class G4ParticleChange;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test30VSecondaryGenerator
{
public:

  Test30VSecondaryGenerator(G4HadronicInteraction* hadi, G4Material* mat);

  virtual ~Test30VSecondaryGenerator();

  virtual G4HadFinalState* Secondaries(const G4Track& track);

  const G4String GeneratorName() const {return generatorName;};

  G4double GetMass() {return mass;};

  void SetA(G4int A) {targetN = A;};

protected:

  G4String generatorName;

private:

  // hide assignment operator as private
  Test30VSecondaryGenerator(const Test30VSecondaryGenerator&);
  Test30VSecondaryGenerator& operator = (const Test30VSecondaryGenerator &right);

  G4HadronicInteraction* hInteraction;
  G4Material* material;
  G4Nucleus targetNucleus;
  G4double mass;
  G4HadFinalState* result;
  G4int targetN;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif  // Test30VSecondaryGenerator_Test30VSecondaryGenerator_h

