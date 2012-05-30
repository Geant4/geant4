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
//---------------------------------------------------------------------
// Class Description:
//
// Intermediate class for hadronic absorption at rest using Bertini
// Physics lists should reference the concrete subclasses for pi-, K-, Sigma-

#include "G4HadronicAbsorptionBertini.hh"
#include "G4CascadeInterface.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include <iostream>


// Constructor

G4HadronicAbsorptionBertini::
G4HadronicAbsorptionBertini(G4ParticleDefinition* pdef)
  : G4HadronStoppingProcess("hBertiniCaptureAtRest"), pdefApplicable(pdef) {
  G4CascadeInterface* cascade = new G4CascadeInterface;
  cascade->SetMinEnergy(0.);			// Ensure it gets used at rest
  cascade->usePreCompoundDeexcitation();
  RegisterMe(cascade);
}


// Applies to constructor-specified particle, or to all known cases

G4bool G4HadronicAbsorptionBertini::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( (0==pdefApplicable && (&particle == G4PionMinus::Definition() ||
				  &particle == G4KaonMinus::Definition() ||
				  &particle == G4SigmaMinus::Definition()))
	   || (&particle == pdefApplicable) 
	   );
}


// Documentation of purpose

void 
G4HadronicAbsorptionBertini::ProcessDescription(std::ostream& os) const {
  os << "Stopping and absorption of charged hadrons (pi-, K-, or Sigma-)\n"
     << "using Bertini-like intranuclear cascade.\n"
     << "Native PreCompound model is used for nuclear de-excitation"
     << std::endl;
}
