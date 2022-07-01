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
/// \file materials/src/G4LatticePhysical.cc
/// \brief Implementation of the G4LatticePhysical class
//
//
// 20131115  Save rotation results in local variable, report verbosely
// 20131116  Replace G4Transform3D with G4RotationMatrix

#include "G4LatticePhysical.hh"
#include "G4LatticeLogical.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"


// Unit vectors defined for convenience (avoid memory churn)

namespace {
  G4ThreeVector xhat(1,0,0), yhat(0,1,0), zhat(0,0,1), nullVec(0,0,0);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LatticePhysical::G4LatticePhysical(const G4LatticeLogical* Lat,
				     const G4RotationMatrix* Rot)
  : verboseLevel(0), fTheta(0), fPhi(0), fLattice(Lat) {
  SetPhysicalOrientation(Rot);
}

G4LatticePhysical::~G4LatticePhysical() {;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LatticePhysical::SetPhysicalOrientation(const G4RotationMatrix* Rot) {
  if(Rot == nullptr)
  {  // No orientation specified
    fLocalToGlobal = fGlobalToLocal = G4RotationMatrix::IDENTITY;
  }
  else
  {
    fLocalToGlobal = fGlobalToLocal = *Rot;		// Frame rotation
    fGlobalToLocal.invert();
  }

  if(verboseLevel != 0)
  {
    G4cout << "G4LatticePhysical::SetPhysicalOrientation " << *Rot
	   << "\nfLocalToGlobal: " << fLocalToGlobal
	   << "\nfGlobalToLocal: " << fGlobalToLocal
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LatticePhysical::SetLatticeOrientation(G4double t_rot, G4double p_rot) {
  fTheta = t_rot;
  fPhi = p_rot;

  if(verboseLevel != 0)
  {
    G4cout << "G4LatticePhysical::SetLatticeOrientation " << fTheta << " "
           << fPhi << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LatticePhysical::SetMillerOrientation(G4int l, G4int k, G4int n) {
  fTheta = halfpi - std::atan2(n+0.000001,l+0.000001);
  fPhi = halfpi - std::atan2(l+0.000001,k+0.000001);

  if(verboseLevel != 0)
  {
    G4cout << "G4LatticePhysical::SetMillerOrientation(" << l << k << n
           << ") : " << fTheta << " " << fPhi << G4endl;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

///////////////////////////////
//Loads the group velocity in m/s
/////////////////////////////
G4double G4LatticePhysical::MapKtoV(G4int polarizationState,
				    G4ThreeVector k) const {
  if(verboseLevel > 1)
  {
    G4cout << "G4LatticePhysical::MapKtoV " << k << G4endl;
  }

  k.rotate(yhat,fTheta).rotate(zhat, fPhi);
  return fLattice->MapKtoV(polarizationState, k);
}

///////////////////////////////
//Loads the normalized direction vector along VG
///////////////////////////////
G4ThreeVector G4LatticePhysical::MapKtoVDir(G4int polarizationState,
					    G4ThreeVector k) const {
  if(verboseLevel > 1)
  {
    G4cout << "G4LatticePhysical::MapKtoVDir " << k << G4endl;
  }

  k.rotate(yhat,fTheta).rotate(zhat,fPhi);

  G4ThreeVector VG = fLattice->MapKtoVDir(polarizationState, k);  

  return VG.rotate(zhat,-fPhi).rotate(yhat,-fTheta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Apply orientation transforms to specified vector

G4ThreeVector 
G4LatticePhysical::RotateToGlobal(const G4ThreeVector& dir) const {
  if (verboseLevel>1) {
    G4cout << "G4LatticePhysical::RotateToGlobal " << dir
	   << "\nusing fLocalToGlobal " << fLocalToGlobal
	   << G4endl;
  }

  G4ThreeVector result = fLocalToGlobal*dir;
  if(verboseLevel > 1)
  {
    G4cout << " result " << result << G4endl;
  }

  return result;
}

G4ThreeVector 
G4LatticePhysical::RotateToLocal(const G4ThreeVector& dir) const {
  if (verboseLevel>1) {
    G4cout << "G4LatticePhysical::RotateToLocal " << dir
	   << "\nusing fGlobalToLocal " << fGlobalToLocal
	   << G4endl;
  }

  G4ThreeVector result = fGlobalToLocal*dir;
  if(verboseLevel > 1)
  {
    G4cout << " result " << result << G4endl;
  }

  return result;
}
