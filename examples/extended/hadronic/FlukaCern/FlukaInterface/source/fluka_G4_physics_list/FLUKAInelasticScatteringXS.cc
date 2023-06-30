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
#ifdef G4_USE_FLUKA


#include "FLUKAInelasticScatteringXS.hh"

// G4
#include "G4DynamicParticle.hh"
#include "G4ElementTable.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4IsotopeList.hh"
#include "G4SystemOfUnits.hh"

#include "fluka_interface.hh"


#define DEBUG false


// ***************************************************************************
// FLUKA hadron nucleus inelastic XS.
// ***************************************************************************
FLUKAInelasticScatteringXS::FLUKAInelasticScatteringXS() 
  : G4VCrossSectionDataSet("FLUKAInelasticScatteringXS")
{
  verboseLevel = 0;
  if (verboseLevel > 0) {
    G4cout << "FLUKAInelasticScatteringXS::FLUKAInelasticScatteringXS constructor." << G4endl;
  }
  
  SetForAllAtomsAndEnergies(true);
}


// ***************************************************************************
// Basic description.
// ***************************************************************************
void FLUKAInelasticScatteringXS::CrossSectionDescription(std::ostream& outFile) const {
  outFile << "FLUKAInelasticScatteringXS: Fluka9x/200x"
          << " to compute inelastic scattering lengths and cross sections.\n";
}


// ***************************************************************************
// Make an element non-applicable,
// so that the element's considered isotopic composition is always the G4 one.
// ***************************************************************************
G4bool FLUKAInelasticScatteringXS::IsElementApplicable(const G4DynamicParticle*, 
						       G4int, 
						       const G4Material*) {

  // For ions, get kinetic energy per nucleon
  /*if (std::abs(particle->GetDefinition()->GetBaryonNumber()) > 1) {
    const G4double kineticEnergyPerNucleon = particle->GetKineticEnergy() / 
    std::abs(particle->GetDefinition()->GetBaryonNumber());
    if (kineticEnergyPerNucleon > GetMaxKinEnergy()) {
		
    G4cout << "Particle kinetic energy = " << particle->GetKineticEnergy() << G4endl;
    G4cout << "GetMaxKinEnergy() = " << GetMaxKinEnergy() << G4endl;
    G4cout << "particle->GetDefinition()->GetBaryonNumber() = " << particle->GetDefinition()->GetBaryonNumber() << G4endl;
    G4cout << "No inelastic interaction model will be available for this ion" << G4endl;
    G4cout << "########################################"
    << "########################################" << G4endl;
    return false;
    }
    }*/

  return false;
}


// ***************************************************************************
// The XS is eventually called at the isotope level.
// ***************************************************************************
G4bool FLUKAInelasticScatteringXS::IsIsoApplicable(const G4DynamicParticle*,
						   G4int, G4int,
						   const G4Element*, 
						   const G4Material*) {

  return true;
}


// ***************************************************************************
// Access point to FLUKA XS at the ELEMENT level 
// (to make it available, depending on the choice made in sElementApplicable).
// ***************************************************************************
G4double FLUKAInelasticScatteringXS::GetElementCrossSection(const G4DynamicParticle* aParticle,
							    G4int targetZ, 
							    const G4Material*) {

  // GET INELASTIC XS FROM FLUKA INTERFACE
  const G4double xs = fluka_interface::computeInelasticScatteringXS(aParticle, 
                                                                    targetZ);
#if DEBUG
  G4cout << "GetElementCrossSection"
         << ", Ekin[MeV] = " << aParticle->GetKineticEnergy() / CLHEP::MeV 
         << ", targetZ = " << targetZ 
         << ", xs[barn] = " << xs / CLHEP::barn 
         << G4endl;
  G4cout << G4endl;
#endif
	
  return xs;
}


// ***************************************************************************
// Access point to FLUKA XS at the ISOTOPE level
// ***************************************************************************
G4double FLUKAInelasticScatteringXS::GetIsoCrossSection(const G4DynamicParticle* aParticle, 
							G4int targetZ, G4int targetA,
							const G4Isotope*, 
							const G4Element*,
							const G4Material*) {
 
  // GET INELASTIC XS FROM FLUKA INTERFACE
  const G4double xs = fluka_interface::computeInelasticScatteringXS(aParticle, 
                                                                    targetZ, 
                                                                    targetA);

#if DEBUG
  G4cout << "GetIsoCrossSection"
         << ", Ekin[MeV] = " << aParticle->GetKineticEnergy() / CLHEP::MeV 
         << ", targetZ = " << targetZ << ", targetA = " << targetA
         << ", xs[barn] = " << xs / CLHEP::barn 
         << G4endl;
  G4cout << G4endl;
#endif

  return xs;
}


#endif // G4_USE_FLUKA
