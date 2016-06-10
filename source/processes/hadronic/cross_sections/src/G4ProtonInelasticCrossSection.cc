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
// By JPW, working, but to be cleaned up. @@@
// G.Folger, 29-sept-2006: extend to 1TeV, using a constant above 20GeV
// 22 Dec 2006 - DHW added isotope dependence
// G.Folger, 25-Nov-2009: extend to 100TeV, using a constant above 20GeV
// V.Ivanchenko, 18-Aug-2011: migration to new design and cleanup;
//                            make it applicable for Z>1
//

#include "G4ProtonInelasticCrossSection.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4HadTmpUtil.hh"
#include "G4NistManager.hh"

G4ProtonInelasticCrossSection::G4ProtonInelasticCrossSection()
  : G4VCrossSectionDataSet("Axen-Wellisch"),thEnergy(19.8*CLHEP::GeV)
{
  nist = G4NistManager::Instance();
  theProton = G4Proton::Proton();
}

G4ProtonInelasticCrossSection::~G4ProtonInelasticCrossSection()
{}

G4bool 
G4ProtonInelasticCrossSection::IsElementApplicable(
                             const G4DynamicParticle* aPart, 
			     G4int Z, const G4Material*)
{
  return ((1 < Z) && (aPart->GetDefinition() == theProton));
}

G4double G4ProtonInelasticCrossSection::GetElementCrossSection(
			     const G4DynamicParticle* aPart, 
			     G4int Z, const G4Material*)
{
  return GetProtonCrossSection(aPart->GetKineticEnergy(), Z);
}

G4double G4ProtonInelasticCrossSection::GetProtonCrossSection(
			     G4double kineticEnergy, G4int Z)
{
  if(kineticEnergy <= 0.0) { return 0.0; }
 
  // constant cross section above ~20GeV
  if (kineticEnergy > thEnergy) { kineticEnergy = thEnergy; }

  G4double a = nist->GetAtomicMassAmu(Z);
  G4double a13 = G4Pow::GetInstance()->powA(a,-0.3333333333);
  G4int nOfNeutrons = G4lrint(a) - Z;
  kineticEnergy /=GeV;
  G4double alog10E = std::log10(kineticEnergy);

  static const G4double nuleonRadius=1.36e-15;
  static const G4double fac=CLHEP::pi*nuleonRadius*nuleonRadius;

  G4double b0   = 2.247-0.915*(1 - a13);
  G4double fac1 = b0*(1 - a13);
  G4double fac2 = 1.;
  if(nOfNeutrons > 1) { fac2=G4Log((G4double(nOfNeutrons))); }
  G4double crossSection = 1.0E31*fac*fac2*(1. + 1./a13 - fac1);

  // high energy correction
  crossSection *= (1 - 0.15*G4Exp(-kineticEnergy))/(1.0 - 0.0007*a);

  // first try on low energies: rise
  G4double ff1= 0.70-0.002*a;  // slope of the drop at medium energies.
  G4double ff2= 1.00+1/a;  // start of the slope.
  G4double ff3= 0.8+18/a-0.002*a; // stephight

  G4double ff4= 1.0 - (1.0/(1+G4Exp(-8*ff1*(alog10E + 1.37*ff2))));

  crossSection *= (1 + ff3*ff4);

  // low energy return to zero

  ff1=1. - 1./a - 0.001*a;       // slope of the rise
  ff2=1.17 - 2.7/a - 0.0014*a;   // start of the rise
 
  ff4=-8.*ff1*(alog10E + 2.0*ff2);
 
  crossSection *= millibarn/(1. + G4Exp(ff4));
  return crossSection;
}
