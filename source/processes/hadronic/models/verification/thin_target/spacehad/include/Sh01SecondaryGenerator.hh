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
//    Sh01SecondaryGenerator
//                
// 
//    Modified:
//
// 06.03.03 V. Grichine (based on test30 og V. Ivanchenko)
//
////////////////////////////////////////////////////////////////////////////

#ifndef Sh01SecondaryGenerator_Sh01SecondaryGenerator_h
#define Sh01SecondaryGenerator_Sh01SecondaryGenerator_h 1

#include <string>
#include <vector>
#include "globals.hh"
#include "G4LorentzVector.hh"
#include "G4VParticleChange.hh"
#include "G4Nucleus.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"


///////////////////////////////////////////////////////////////////////

class G4Track;
class G4ParticleDefinition;
class G4HadronicInteraction;
class G4Material;
class G4VParticleChange;
class G4ParticleChange;

/////////////////////////////////////////////////////////////////////

class Sh01SecondaryGenerator 
{
public:     

  Sh01SecondaryGenerator(G4HadronicInteraction* hadi, G4Material* mat);

  virtual ~Sh01SecondaryGenerator();

  // virtual G4VParticleChange* Secondaries(const G4Track& track);
  virtual G4HadFinalState* Secondaries(const G4Track& track);

  const G4String GeneratorName() const {return generatorName;};
			   
  G4double GetMass() {return mass;};
	 
protected:


  //  G4VParticleChange theParticleChange;
  // the G4VParticleChange object which is modified and returned

  G4String generatorName;

private:
  
  // hide assignment operator as private 
  Sh01SecondaryGenerator(const Sh01SecondaryGenerator&);
  Sh01SecondaryGenerator& operator = (const Sh01SecondaryGenerator &right);

  G4HadronicInteraction* hInteraction;
  G4Material* material;
  G4Nucleus targetNucleus;
  G4double mass;
  G4HadFinalState* result;

};

//////////////////////////////////////////////////////////////////////

#endif  // Sh01SecondaryGenerator_Sh01SecondaryGenerator_h

