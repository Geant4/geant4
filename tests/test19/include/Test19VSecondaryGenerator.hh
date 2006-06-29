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
// ---------------------------------------------------------------------
// Short description of the G4 Class made by M.K. (asFarAsHeNnderstands)
// =====================================================================
// This is a tool of the HadronicModelCollection, which does not take the
// Material from the G4Track, but defines it independentli and extracts Z
// and (in a wrong way by rounding the mean value) N of the target nucleus.
// When the SecondaryGenerator and the Nucleus are set in the conctructer,
// one can get secondaries in a form of the HadronicInternal G4HadFinalState,
// which is converted to the standard G4VParticleChange in the Test19HadronProduction 
// class. The projectile parameters are still taken from the G4Track information.
// *** M.K. Why the material is not taken from the G4Track information? ***
// *** the "mass" returns a mass of the selected nucleus (?), which is fixed ***
// *** Hadronic application does not randomize the Isotope of the Element! ***
//
// ======================================================================================
//
//      ---------- Test19VSecondaryGenerator -------
// 
//    Converted from Test29 to Test19 by Mikhail Kossov, 7 Jan 2005 
//
//---------------------------------------------------------------------------------------

#ifndef Test19VSecondaryGenerator_Test19VSecondaryGenerator_h
#define Test19VSecondaryGenerator_Test19VSecondaryGenerator_h 1

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

class Test19VSecondaryGenerator
{
public:

  Test19VSecondaryGenerator(G4HadronicInteraction* hadi, G4Material* mat);
  virtual ~Test19VSecondaryGenerator();

  // Modifier
  virtual G4HadFinalState* Secondaries(const G4Track& track); // Starts the generator

  //Selectors
  const G4String GeneratorName() const {return generatorName;};
  G4double GetMass() {return mass;};                          // Mass of the target nucleus

protected:

  G4String generatorName; // If GeneratorName() exists, why it is not private?

private:

  // hide assignment operator as private
  Test19VSecondaryGenerator(const Test19VSecondaryGenerator&);
  Test19VSecondaryGenerator& operator = (const Test19VSecondaryGenerator &right);

  G4HadronicInteraction* hInteraction;
  G4Material* material;
  G4Nucleus targetNucleus;
  G4double mass;
  G4HadFinalState* result;

};

#endif  // Test19VSecondaryGenerator_Test19VSecondaryGenerator_h

