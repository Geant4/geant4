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
// File name:     G4VPolarizedXS
//
// Author:        Andreas Schaelicke
//
// Class Description:
//   (virtual) interface class
//   provides readable but efficient routines to determine
//   polarization for the final state of a given process
//   empoying the differential cross section

#include "G4VPolarizedXS.hh"

#include "Randomize.hh"

G4VPolarizedXS::G4VPolarizedXS()
  : fXmin(0)
  , fXmax(1.)
  , fYmin(1.)
  , fA(1.)
  , fZ(1.)
  , fCoul(0.)
{}

G4VPolarizedXS::~G4VPolarizedXS() {}

void G4VPolarizedXS::Initialize(G4double, G4double, G4double,
                                const G4StokesVector&, const G4StokesVector&,
                                G4int)
{}

G4StokesVector G4VPolarizedXS::GetPol2()
{
  // neglects correlation effects!
  G4double invXsecTotal =
    1. / XSection(G4StokesVector::ZERO, G4StokesVector::ZERO);
  G4double xsPol1 = XSection(G4StokesVector::P1, G4StokesVector::ZERO);
  G4double xsPol2 = XSection(G4StokesVector::P2, G4StokesVector::ZERO);
  G4double xsPol3 = XSection(G4StokesVector::P3, G4StokesVector::ZERO);
  return G4StokesVector(G4ThreeVector(
    invXsecTotal * xsPol1, invXsecTotal * xsPol2, invXsecTotal * xsPol3));
}

G4StokesVector G4VPolarizedXS::GetPol3()
{
  // neglects correlation effects!
  G4double invXsecTotal =
    1. / XSection(G4StokesVector::ZERO, G4StokesVector::ZERO);
  G4double xsPol1 = XSection(G4StokesVector::ZERO, G4StokesVector::P1);
  G4double xsPol2 = XSection(G4StokesVector::ZERO, G4StokesVector::P2);
  G4double xsPol3 = XSection(G4StokesVector::ZERO, G4StokesVector::P3);
  return G4StokesVector(G4ThreeVector(
    invXsecTotal * xsPol1, invXsecTotal * xsPol2, invXsecTotal * xsPol3));
}

// minimal energy fraction in TotalXSection
G4double G4VPolarizedXS::GetXmin(G4double /*y*/) { return fXmin; }

// maximal energy fraction in TotalXSection
G4double G4VPolarizedXS::GetXmax(G4double /*y*/) { return fXmax; }

G4double G4VPolarizedXS::TotalXSection(G4double, G4double, G4double,
                                       const G4StokesVector&,
                                       const G4StokesVector&)
{
  G4ExceptionDescription ed;
  ed << "WARNING virtual function G4VPolarizedXS::TotalXSection() "
        "called.\n";
  G4Exception("G4VPolarizedXS::TotalXSection", "pol032", FatalException, ed);
  return 0.;
}
