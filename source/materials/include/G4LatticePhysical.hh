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
/// \file materials/include/G4LatticePhysical.hh
/// \brief Definition of the G4LatticePhysical class
//
//
// 20131114  Add verbosity for diagnostic output
// 20131116  Replace G4Transform3D with G4RotationMatrix

#ifndef G4LatticePhysical_h
#define G4LatticePhysical_h 1

#include "G4LatticeLogical.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"


class G4LatticePhysical {
public:
 G4LatticePhysical(
   const G4LatticeLogical* Lat = nullptr,
   const G4RotationMatrix* Rot = nullptr);  // Use FRAME rotation
 virtual ~G4LatticePhysical();

 void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

 // NOTE:  Pass by value to allow in-situ rotations
 G4double MapKtoV(G4int, G4ThreeVector) const;
 G4ThreeVector MapKtoVDir(G4int, G4ThreeVector) const;

 void SetLatticeLogical(const G4LatticeLogical* Lat) { fLattice = Lat; }
 void SetPhysicalOrientation(const G4RotationMatrix* Rot);  // FRAME rotation
 void SetLatticeOrientation(G4double, G4double);
 void SetMillerOrientation(G4int, G4int, G4int);

public:  
  const G4LatticeLogical* GetLattice() const { return fLattice; }

  G4double GetScatteringConstant() const { return fLattice->GetScatteringConstant(); }
  G4double GetAnhDecConstant() const { return fLattice->GetAnhDecConstant(); }
  G4double GetLDOS() const           { return fLattice->GetLDOS(); }
  G4double GetSTDOS() const          { return fLattice->GetSTDOS(); }
  G4double GetFTDOS() const          { return fLattice->GetFTDOS(); }
  G4double GetBeta() const           { return fLattice->GetBeta(); }
  G4double GetGamma() const          { return fLattice->GetGamma(); }
  G4double GetLambda() const         { return fLattice->GetLambda(); }
  G4double GetMu() const             { return fLattice->GetMu(); }

  // Apply orientation transforms to specified vector
  G4ThreeVector RotateToGlobal(const G4ThreeVector& dir) const;
  G4ThreeVector RotateToLocal(const G4ThreeVector& dir) const;

private:
  G4int verboseLevel;			// Enable diagnostic output

  G4double fTheta, fPhi;		// Lattice orientation within object
  const G4LatticeLogical* fLattice;	// Underlying lattice parameters

  G4RotationMatrix fLocalToGlobal;
  G4RotationMatrix fGlobalToLocal;
};

#endif
