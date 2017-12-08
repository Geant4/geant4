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
// ABLAXX statistical de-excitation model
// Jose Luis Rodriguez, CEA (translation from ABLA07 and contact person)
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Davide Mancusi, CEA (contact person INCL)
// Aatos Heikkinen, HIP (project coordination)
//
#define ABLAXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifdef ABLAXX_IN_GEANT4_MODE

#include "G4AblaInterface.hh"
#include "G4ParticleDefinition.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4DynamicParticle.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
//#include "G4INCLParticleTable.hh"
//#include "G4INCLGlobals.hh"
#include <iostream>
#include <cmath>

G4AblaInterface::G4AblaInterface() :
  G4VPreCompoundModel(NULL, "ABLA"),
  ablaResult(new G4VarNtp),
  volant(new G4Volant),
  theABLAModel(new G4Abla(volant, ablaResult)),
  eventNumber(0)
{
  theABLAModel->initEvapora();
  theABLAModel->SetParameters();
}

G4AblaInterface::~G4AblaInterface() {
  delete volant;
  delete ablaResult;
  delete theABLAModel;
}

G4ReactionProductVector *G4AblaInterface::DeExcite(G4Fragment &aFragment) {
  volant->clear();
  ablaResult->clear();

  const G4int ARem = aFragment.GetA_asInt();
  const G4int ZRem = aFragment.GetZ_asInt();
  const G4double eStarRem = aFragment.GetExcitationEnergy() / MeV;
  const G4double jRem = aFragment.GetAngularMomentum().mag() / hbar_Planck;
  const G4LorentzVector &pRem = aFragment.GetMomentum();
  const G4double pxRem = pRem.x() / MeV;
  const G4double pyRem = pRem.y() / MeV;
  const G4double pzRem = pRem.z() / MeV;

  eventNumber++;

  theABLAModel->DeexcitationAblaxx(ARem, ZRem,
                          eStarRem,
                          jRem,
                          pxRem,
                          pyRem,
                          pzRem,
                          eventNumber);

  G4ReactionProductVector *result = new G4ReactionProductVector;

  for(int j = 0; j < ablaResult->ntrack; ++j) { // Copy ABLA result to the EventInfo
    G4ReactionProduct *product = toG4Particle(ablaResult->avv[j],
                                              ablaResult->zvv[j],
                                              ablaResult->enerj[j],
                                              ablaResult->plab[j]*std::sin(ablaResult->tetlab[j]*pi/180.0)*std::cos(ablaResult->philab[j]*pi/180.0),
                                              ablaResult->plab[j]*std::sin(ablaResult->tetlab[j]*pi/180.0)*std::sin(ablaResult->philab[j]*pi/180.0),
                                              ablaResult->plab[j]*std::cos(ablaResult->tetlab[j]*pi/180.0));
    if(product)
      result->push_back(product);
  }
  return result;
}

G4ParticleDefinition *G4AblaInterface::toG4ParticleDefinition(G4int A, G4int Z) const {
  if     (A == 1 && Z == 1)  return G4Proton::Proton();
  else if(A == 1 && Z == 0)  return G4Neutron::Neutron();
  else if(A == -1 && Z == 1)  return G4PionPlus::PionPlus();
  else if(A == -1 && Z == -1) return G4PionMinus::PionMinus();
  else if(A == -1 && Z == 0)  return G4PionZero::PionZero();
  else if(A == 0 && Z == 0)  return G4Gamma::Gamma();
  else if(A == 2 && Z == 1)  return G4Deuteron::Deuteron();
  else if(A == 3 && Z == 1)  return G4Triton::Triton();
  else if(A == 3 && Z == 2)  return G4He3::He3();
  else if(A == 4 && Z == 2)  return G4Alpha::Alpha();
  else if(A > 0 && Z > 0 && A >= Z) { // Returns ground state ion definition
    return G4IonTable::GetIonTable()->GetIon(Z, A, 0);
  } else { // Error, unrecognized particle
    G4cout << "Can't convert particle with A=" << A << ", Z=" << Z << " to G4ParticleDefinition, trouble ahead" << G4endl;
    return 0;
  }
}

G4ReactionProduct *G4AblaInterface::toG4Particle(G4int A, G4int Z,
						 G4double kinE,
						 G4double px,
                                                 G4double py, G4double pz) const {
  G4ParticleDefinition *def = toG4ParticleDefinition(A, Z);
  if(def == 0) { // Check if we have a valid particle definition
    return 0;
  }
  const double energy = kinE * MeV;
  const G4ThreeVector momentum(px, py, pz);
  const G4ThreeVector momentumDirection = momentum.unit();
  G4DynamicParticle p(def, momentumDirection, energy);
  G4ReactionProduct *r = new G4ReactionProduct(def);
  (*r) = p;
  return r;
}

void G4AblaInterface::ModelDescription(std::ostream& outFile) const {
   outFile << "ABLA++ does not provide an implementation of the ApplyYourself method!\n\n";
}

void G4AblaInterface::DeExciteModelDescription(std::ostream& outFile) const {
   outFile << "ABLA++ is a statistical model for nuclear de-excitation. It simulates\n"
     << "evaporation of neutrons, protons and alpha particles, as well as fission\n"
     << "where applicable. The code included in Geant4 is a C++ translation of the\n"
     << "original Fortran code. More details about the physics are available in the\n"
     << "the Geant4 Physics Reference Manual and in the reference articles.\n\n"
     << "Reference:\n"
     << "A. Kelic, M. V. Ricciardi, and K. H. Schmidt, in Proceedings of Joint ICTP-IAEA Advanced Workshop on Model Codes for Spallation Reactions, ICTP Trieste, Italy, 4–8 February 2008, edited by D. Filges, S. Leray, Y. Yariv, A. Mengoni, A. Stanculescu, and G. Mank (IAEA INDC(NDS)-530, Vienna, 2008), pp. 181–221.\n\n"; }

#endif // ABLAXX_IN_GEANT4_MODE
