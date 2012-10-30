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
// INCL++ intra-nuclear cascade model
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
// INCL++ revision: v5.1.5
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLNuclearDensity.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  NuclearDensity::NuclearDensity(G4int A, G4int Z, InverseInterpolationTable *rpCorrelationTable) :
    theA(A),
    theZ(Z),
    theMaximumRadius((*rpCorrelationTable)(1.)),
    rFromP(rpCorrelationTable),
    // The interpolation table for local-energy look-ups is simply obtained by
    // inverting the r-p correlation table.
    tFromR(new InverseInterpolationTable(rFromP->getNodeValues(), rFromP->getNodeAbscissae()))
  {
    DEBUG("Interpolation table for local energy (A=" << theA << ", Z=" << theZ << ") initialised:"
        << std::endl
        << tFromR->print()
        << std::endl);
    initCentralRadius();
    initializeTransmissionRadii();
  }

  NuclearDensity::~NuclearDensity() {
    // We don't delete the rFromP table, which is cached in the
    // NuclearDensityFactory
    delete tFromR;
  }

  NuclearDensity::NuclearDensity(const NuclearDensity &rhs) :
    theA(rhs.theA),
    theZ(rhs.theZ),
    theMaximumRadius(rhs.theMaximumRadius),
    theCentralRadius(rhs.theCentralRadius),
    transmissionRadius(rhs.transmissionRadius),
    // rFromP is owned by NuclearDensityFactory, so shallow copy is sufficient
    rFromP(rhs.rFromP)
  {
    // deep copy for tFromR
    delete tFromR;
    tFromR = new InverseInterpolationTable(*(rhs.tFromR));
  }

  NuclearDensity &NuclearDensity::operator=(const NuclearDensity &rhs) {
    NuclearDensity temporaryDensity(rhs);
    swap(temporaryDensity);
    return *this;
  }

  void NuclearDensity::swap(NuclearDensity &rhs) {
    std::swap(theA, rhs.theA);
    std::swap(theZ, rhs.theZ);
    std::swap(theMaximumRadius, rhs.theMaximumRadius);
    std::swap(theCentralRadius, rhs.theCentralRadius);
    std::swap(transmissionRadius, rhs.transmissionRadius);
    std::swap(rFromP, rhs.rFromP);
    std::swap(tFromR, rhs.tFromR);
 }

  void NuclearDensity::initializeTransmissionRadii() {
    const G4double r0 = 1.12;
    const G4double theNucleonTransmissionRadius = r0*Math::pow13((G4double)theA) + 0.88;

    transmissionRadius[Proton] = theNucleonTransmissionRadius;
    transmissionRadius[Neutron] = theNucleonTransmissionRadius;
    transmissionRadius[PiPlus] = theCentralRadius;
    transmissionRadius[PiZero] = theCentralRadius;
    transmissionRadius[PiMinus] = theCentralRadius;
    transmissionRadius[DeltaPlusPlus] = theNucleonTransmissionRadius;
    transmissionRadius[DeltaPlus] = theNucleonTransmissionRadius;
    transmissionRadius[DeltaZero] = theNucleonTransmissionRadius;
    transmissionRadius[DeltaMinus] = theNucleonTransmissionRadius;
    transmissionRadius[Composite] = theCentralRadius;
  }

  G4double NuclearDensity::getMaxRFromP(G4double p) const {
    return (*rFromP)(p);
  }

  G4double NuclearDensity::getMaxTFromR(G4double r) const {
    return (*tFromR)(r);
  }

}
