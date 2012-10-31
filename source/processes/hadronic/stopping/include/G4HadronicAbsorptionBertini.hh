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
// 20121017  M. Kelsey -- Add local cache of Bertini pointer

#ifndef G4HadronicAbsorptionBertini_h
#define G4HadronicAbsorptionBertini_h 1

// Class Description:
//
// Intermediate base (or concrete) class for hadronic absorption at rest. 
// Physics lists should reference the concrete subclasses for pi-, K-, Sigma-

#include "globals.hh"
#include "G4HadronStoppingProcess.hh"
#include <iosfwd>

class G4VParticleChange;
class G4ParticleDefinition;
class G4CascadeInterface;
class G4Track;


class G4HadronicAbsorptionBertini : public G4HadronStoppingProcess { 
public:
  // May instantiate this class directly for all three pi-, K-, Sigma-
  G4HadronicAbsorptionBertini(G4ParticleDefinition* pdef=0);
  virtual ~G4HadronicAbsorptionBertini() {;}
  
  G4bool IsApplicable(const G4ParticleDefinition&);

  void ProcessDescription(std::ostream& outFile) const;

private:
  // hide assignment operator as private 
  G4HadronicAbsorptionBertini& operator=(const G4HadronicAbsorptionBertini&);
  G4HadronicAbsorptionBertini(const G4HadronicAbsorptionBertini&);
  
private:
  G4ParticleDefinition* pdefApplicable;
  G4CascadeInterface* theCascade;
};

#endif

