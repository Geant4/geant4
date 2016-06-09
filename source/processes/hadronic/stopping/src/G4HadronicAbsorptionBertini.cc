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
//
// 20120905  M. Kelsey -- Drop explicit list of "allowed" particles; Bertini
//		can handle anything, or return no-interaction if not.
// 20121017  M. Kelsey -- Use Bertini's IsApplicable to check particle allowed

#include "G4HadronicAbsorptionBertini.hh"
#include "G4CascadeInterface.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include <iostream>


// Constructor

G4HadronicAbsorptionBertini::
G4HadronicAbsorptionBertini(G4ParticleDefinition* pdef)
  : G4HadronStoppingProcess("hBertiniCaptureAtRest"), pdefApplicable(pdef) {
  theCascade = new G4CascadeInterface;
  theCascade->SetMinEnergy(0.);			// Ensure it gets used at rest
  theCascade->usePreCompoundDeexcitation();
  RegisterMe(theCascade);			// Transfers ownership
}


// Applies to constructor-specified particle, or to all known cases

G4bool G4HadronicAbsorptionBertini::IsApplicable(const G4ParticleDefinition& particle)
{
  // Exclusive match (if registered for specific projectile
  if (pdefApplicable) return (&particle == pdefApplicable);

  // Any negative particles known to Bertini, excluding nuclei
  return (G4HadronStoppingProcess::IsApplicable(particle) &&
	  particle.GetAtomicMass() <= 1 &&
	  theCascade->IsApplicable(&particle));
}


// Documentation of purpose

void 
G4HadronicAbsorptionBertini::ProcessDescription(std::ostream& os) const {
  os << "Stopping and absorption of charged hadrons (pi-, K-, Sigma-)\n"
     << "using Bertini-like intranuclear cascade.\n"
     << "Native PreCompound model is used for nuclear de-excitation"
     << std::endl;
}
