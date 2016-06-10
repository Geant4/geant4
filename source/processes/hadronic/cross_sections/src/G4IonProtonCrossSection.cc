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
// $Id: G4IonProtonCrossSection.cc 90447 2015-05-29 07:41:53Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:  G4IonProtonCrossSection
//
// Author  Ivantchenko, Geant4, 30 July 2010
//
// Modifications:
//

#include "G4IonProtonCrossSection.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4HadTmpUtil.hh"

using namespace std;

G4IonProtonCrossSection::G4IonProtonCrossSection() 
  : G4VCrossSectionDataSet("AxenWellischIonH") 
{
  theForward = new G4ProtonInelasticCrossSection();
}

G4IonProtonCrossSection::~G4IonProtonCrossSection()
{}

G4bool 
G4IonProtonCrossSection::IsElementApplicable(const G4DynamicParticle* dp, 
					     G4int Z, const G4Material*)
{
  return ((1 == Z) && (dp->GetDefinition()->GetPDGCharge()/eplus > 1.5)); 
}


G4double 
G4IonProtonCrossSection::GetElementCrossSection(
                                   const G4DynamicParticle* dp, 
				   G4int, const G4Material*)
{
  const G4ParticleDefinition* p = dp->GetDefinition();
  G4double e = dp->GetKineticEnergy()*CLHEP::proton_mass_c2/p->GetPDGMass();
  return 
    theForward->GetProtonCrossSection(e, G4lrint(p->GetPDGCharge()/eplus));
}

void 
G4IonProtonCrossSection::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4IonProtonCrossSection calculates the inelastic cross section\n"
          << "for any ion projectile with Z >=2 only on hydrogen target.\n"
          << "It uses the inverse kinematics and the Axen-Wellisch\n"
          << "inelastic cross section (G4ProtonInelasticCrossSection).\n"; 
}

