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
//
// G4PSStepChecker3D
#include "G4PSStepChecker3D.hh"

//////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for Debug.
//
// Created: 2011-03-24  Tsukasa ASO
//
///////////////////////////////////////////////////////////////////////////////

G4PSStepChecker3D::G4PSStepChecker3D(G4String name, G4int ni, G4int nj,
                                     G4int nk, G4int depi, G4int depj,
                                     G4int depk)
  : G4PSStepChecker(name)
  , fDepthi(depi)
  , fDepthj(depj)
  , fDepthk(depk)
{
  SetNijk(ni, nj, nk);
}

G4int G4PSStepChecker3D::GetIndex(G4Step* aStep)
{
  const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int i                       = touchable->GetReplicaNumber(fDepthi);
  G4int j                       = touchable->GetReplicaNumber(fDepthj);
  G4int k                       = touchable->GetReplicaNumber(fDepthk);

  G4int N = i * fNj * fNk + j * fNk + k;

  G4cout << " depi= " << fDepthi << " depj= " << fDepthj << " depk= " << fDepthk
         << G4endl;
  G4cout << "    i= " << i << "   j= " << j << "    k= " << k << G4endl;
  G4cout << "    N= " << N << "  Nx= " << fNi << " Nj= " << fNj
         << " Nk= " << fNk << G4endl;

  return i * fNj * fNk + j * fNk + k;
}
