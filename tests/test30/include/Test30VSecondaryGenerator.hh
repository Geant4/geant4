//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
#include "G4VParticleChange.hh"
#include "G4Nucleus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Track;
class G4ParticleDefinition;
class G4HadronicInteraction;
class G4Material;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test30VSecondaryGenerator 
{
public:     

  Test30VSecondaryGenerator(G4HadronicInteraction* hadi, G4Material* mat);

  virtual ~Test30VSecondaryGenerator();

  virtual G4VParticleChange* Secondaries(const G4Track& track);

  const G4String GeneratorName() const {return generatorName;};
			   
  G4double GetMass() {return mass;};
	 
protected:


  G4VParticleChange theParticleChange;
  // the G4VParticleChange object which is modified and returned

  G4String generatorName;

private:
  
  // hide assignment operator as private 
  Test30VSecondaryGenerator(const Test30VSecondaryGenerator&);
  Test30VSecondaryGenerator& operator = (const Test30VSecondaryGenerator &right);

  G4HadronicInteraction* hInteraction;
  G4Material* material;
  G4Nucleus targetNucleus;
  G4double mass;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif  // Test30VSecondaryGenerator_Test30VSecondaryGenerator_h

