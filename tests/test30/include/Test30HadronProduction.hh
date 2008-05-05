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
//      ---------- Test30HadronProduction -------
//                by Vladimir Ivanchenko, 12 March 2002 
// 
//    Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#ifndef Test30HadronProduction_Test30HadronProduction_h
#define Test30HadronProduction_Test30HadronProduction_h 1


#include "Test30VSecondaryGenerator.hh"

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4HadFinalState.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Step;
class G4Track;
class G4VParticleChange;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test30HadronProduction : public G4VDiscreteProcess
{
public:

  Test30HadronProduction(const G4String& processName = "HadronProduction" );

  ~Test30HadronProduction();

  void SetSecondaryGenerator(Test30VSecondaryGenerator*);

  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
						G4double   previousStepSize,
						G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  G4bool IsApplicable(const G4ParticleDefinition&) {return true;};

  G4double GetMass() {return theGenerator->GetMass();};

  void SetA(G4int A) {if(theGenerator) theGenerator->SetA(A);};

protected:

  G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*)
  {return DBL_MAX;};

private:

  // hide assignment operator as private
  Test30HadronProduction(const Test30HadronProduction&);
  Test30HadronProduction& operator = (const Test30HadronProduction &right);

  void InitializeMe();

  Test30VSecondaryGenerator* theGenerator;
  G4VParticleChange          theChange;

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif  // Test30HadronProduction_Test30HadronProduction_h

