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

#include "G4IonProtonCrossSection.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleInelasticXS.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Log.hh"

G4IonProtonCrossSection::G4IonProtonCrossSection() 
  : G4VCrossSectionDataSet("InvProtonXS") 
{
  xsProton = new G4ParticleInelasticXS(G4Proton::Proton());
  xsDeuteron = new G4ParticleInelasticXS(G4Deuteron::Deuteron());
  xsTriton = new G4ParticleInelasticXS(G4Triton::Triton());
  xsHe3 = new G4ParticleInelasticXS(G4He3::He3());
  xsAlpha = new G4ParticleInelasticXS(G4Alpha::Alpha());
}

G4IonProtonCrossSection::~G4IonProtonCrossSection()
{}

G4double 
G4IonProtonCrossSection::GetProtonCrossSection(const G4DynamicParticle* dp, G4int Z)
{
  return xsProton->GetElementCrossSection(dp, Z);
}

G4double 
G4IonProtonCrossSection::GetDeuteronCrossSection(const G4DynamicParticle* dp, G4int Z)
{
  return xsDeuteron->GetElementCrossSection(dp, Z);
}

G4double 
G4IonProtonCrossSection::GetTritonCrossSection(const G4DynamicParticle* dp, G4int Z)
{
  return xsTriton->GetElementCrossSection(dp, Z);
}

G4double 
G4IonProtonCrossSection::GetHe3CrossSection(const G4DynamicParticle* dp, G4int Z)
{
  return xsHe3->GetElementCrossSection(dp, Z);
}

G4double 
G4IonProtonCrossSection::GetAlphaCrossSection(const G4DynamicParticle* dp, G4int Z)
{
  return xsAlpha->GetElementCrossSection(dp, Z);
}

G4bool G4IonProtonCrossSection::IsElementApplicable(const G4DynamicParticle*,
                                                    G4int Z, const G4Material*)
{
  return (Z <= 2);
}

G4bool G4IonProtonCrossSection::IsIsoApplicable(const G4DynamicParticle*,
                                                G4int Z, G4int A,    
			                        const G4Element*,
			                        const G4Material*)
{
  return (Z <= 2 && A <= 4);
}

G4double 
G4IonProtonCrossSection::GetElementCrossSection(const G4DynamicParticle* dp, 
				                G4int Z, const G4Material*)
{
  G4int A = (1 == Z) ? 1 : 4;
  return GetIsoCrossSection(dp, Z, A);
}

G4double 
G4IonProtonCrossSection::GetIsoCrossSection(const G4DynamicParticle* dp,
                                            G4int Z, G4int A,  
			                    const G4Isotope*,
			                    const G4Element*,
			                    const G4Material*)
{
  const G4ParticleDefinition* p = dp->GetDefinition();
  G4double e = dp->GetKineticEnergy()*CLHEP::proton_mass_c2/p->GetPDGMass();
  G4double res = 0.0;
  if(1 == Z && 1 == A) {
    res = xsProton->IsoCrossSection(e, G4Log(e), Z, A);
  } else if(1 == Z && 2 == A) {
    res = xsDeuteron->IsoCrossSection(e, G4Log(e), Z, A);
  } else if(1 == Z && 3 == A) {
    res = xsTriton->IsoCrossSection(e, G4Log(e), Z, A);
  } else if(2 == Z && 3 == A) {
    res = xsHe3->IsoCrossSection(e, G4Log(e), Z, A);
  } else if(2 == Z && 4 == A) {
    res = xsAlpha->IsoCrossSection(e, G4Log(e), Z, A);
  }
  return res;
}

void G4IonProtonCrossSection::BuildPhysicsTable(const G4ParticleDefinition&)
{
  xsProton->BuildPhysicsTable(*G4Proton::Proton());
  xsDeuteron->BuildPhysicsTable(*G4Deuteron::Deuteron());
  xsTriton->BuildPhysicsTable(*G4Triton::Triton());
  xsHe3->BuildPhysicsTable(*G4He3::He3());
  xsAlpha->BuildPhysicsTable(*G4Alpha::Alpha());
} 

void 
G4IonProtonCrossSection::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4IonProtonCrossSection calculates the inelastic cross section\n"
          << "for any ion projectile with Z >=2 only on hydrogen target.\n"
          << "It uses the inverse kinematics and the G4ParticleInelasticXS\n"
          << "cross section.\n"; 
}
