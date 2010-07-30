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
// $Id: G4IonProtonCrossSection.cc,v 1.2 2010-07-30 14:33:33 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4ProtonInelasticCrossSection.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"

using namespace std;

G4IonProtonCrossSection::G4IonProtonCrossSection() 
  : G4VCrossSectionDataSet("AxenWellischIonH") 
{
  theForward = new G4ProtonInelasticCrossSection();
}

G4IonProtonCrossSection::~G4IonProtonCrossSection()
{
  delete theForward;
}

G4bool 
G4IonProtonCrossSection::IsApplicable(const G4DynamicParticle* dp, 
				      const G4Element* elm)
{
  return IsZAApplicable(dp, elm->GetZ(), elm->GetN());
}

G4bool 
G4IonProtonCrossSection::IsZAApplicable(const G4DynamicParticle* dp,
					G4double Z, G4double A)
{
  G4bool result = false;
  if(Z < 1.1 && A < 1.1 && dp->GetDefinition()->GetPDGCharge()/eplus > 2.5) 
    { result = true;}
  return result;
}

G4double 
G4IonProtonCrossSection::GetCrossSection(const G4DynamicParticle* dp, 
					 const G4Element*, G4double)
{
  return GetIsoZACrossSection(dp);
}

G4double 
G4IonProtonCrossSection::GetIsoZACrossSection(const G4DynamicParticle* dp, 
					      G4double /*ZZ*/, G4double /*AA*/, 
					      G4double /*temperature*/)
{
  const G4ParticleDefinition* p = dp->GetDefinition();
  G4double e = dp->GetKineticEnergy()*proton_mass_c2/p->GetPDGMass();
  return theForward->GetCrossSection(e, p->GetBaryonNumber(),
				     p->GetPDGCharge()/eplus);
}

void G4IonProtonCrossSection::BuildPhysicsTable(const G4ParticleDefinition&)
{}

void G4IonProtonCrossSection::DumpPhysicsTable(const G4ParticleDefinition&)
{
  G4cout << "G4IonProtonCrossSection: " << GetName() << " uses formula"
	 <<G4endl;
}

