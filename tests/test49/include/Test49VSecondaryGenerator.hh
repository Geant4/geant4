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
// ---------------------------------------------------------------------
// Short description of the G4 Class made by M.K. (asFarAsHeNnderstands)
// =====================================================================
// This is a tool of the HadronicModelCollection, which does not take the
// Material from the G4Track, but defines it independentli and extracts Z
// and (in a wrong way by rounding the mean value) N of the target nucleus.
// When the SecondaryGenerator and the Nucleus are set in the conctructer,
// one can get secondaries in a form of the HadronicInternal G4HadFinalState,
// which is converted to the standard G4VParticleChange in the Test49HadronProduction 
// class. The projectile parameters are still taken from the G4Track information.
// *** M.K. Why the material is not taken from the G4Track information? ***
// *** the "mass" returns a mass of the selected nucleus (?), which is fixed ***
// *** Hadronic application does not randomize the Isotope of the Element! ***
//
// ======================================================================================
//
//      ---------- Test49VSecondaryGenerator -------
// 
//    Converted from Test29 to Test49 by Mikhail Kossov, 7 Jan 2005 
//
//---------------------------------------------------------------------------------------

#ifndef Test49VSecondaryGenerator_Test49VSecondaryGenerator_h
#define Test49VSecondaryGenerator_Test49VSecondaryGenerator_h 1

#include <string>
#include <vector>
#include "globals.hh"
#include "G4LorentzVector.hh"
#include "G4Nucleus.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"

class G4Track;
class G4ParticleDefinition;
class G4HadronicInteraction;
class G4Material;
class G4VParticleChange;
class G4ParticleChange;

class Test49VSecondaryGenerator
{
public:

  Test49VSecondaryGenerator(G4HadronicInteraction* hadi, G4Material* mat);
  virtual ~Test49VSecondaryGenerator();

  // Modifier
  virtual G4HadFinalState* Secondaries(const G4Track& track); // Starts the generator

  //Selectors
  const G4String GeneratorName() const {return generatorName;};
  G4double GetMass() {return mass;};                          // Mass of the target nucleus

protected:

  G4String generatorName; // If GeneratorName() exists, why it is not private?

private:

  // hide assignment operator as private
  Test49VSecondaryGenerator(const Test49VSecondaryGenerator&);
  Test49VSecondaryGenerator& operator = (const Test49VSecondaryGenerator &right);

  G4HadronicInteraction* hInteraction;
  G4Material* material;
  G4Nucleus targetNucleus;
  G4double mass;
  G4HadFinalState* result;

};

#endif  // Test49VSecondaryGenerator_Test49VSecondaryGenerator_h

