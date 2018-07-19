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
// $Id:$
//
//
// Implementation of G4LogicalCrystalVolume
//
// 21-04-16, created by E.Bagli
// 
// --------------------------------------------------------------------

#include "G4LogicalCrystalVolume.hh"
#include "G4ExtendedMaterial.hh"
#include "G4CrystalExtension.hh"
#include "G4VMaterialExtension.hh"

std::vector<G4LogicalVolume*> G4LogicalCrystalVolume::fLCVvec;

// --------------------------------------------------------------------

G4LogicalCrystalVolume::
G4LogicalCrystalVolume(G4VSolid* pSolid, G4ExtendedMaterial* pMaterial,
                       const G4String& name, G4FieldManager* pFieldMgr,
                       G4VSensitiveDetector* pSDetector,
                       G4UserLimits* pULimits, G4bool optimise,
                       G4int h, G4int k, G4int l, G4double rot)
: G4LogicalVolume(pSolid,pMaterial,name,pFieldMgr,pSDetector,pULimits,optimise),
  hMiller(1), kMiller(0), lMiller(0), fRot(0), verboseLevel(0)
{
   SetMillerOrientation(h, k, l, rot);
   fLCVvec.push_back(this);
}

// --------------------------------------------------------------------

G4LogicalCrystalVolume::~G4LogicalCrystalVolume()
{
   fLCVvec.erase( std::remove(fLCVvec.begin(),fLCVvec.end(), this ),
                  fLCVvec.end() );
}

// --------------------------------------------------------------------

G4bool G4LogicalCrystalVolume::IsLattice(G4LogicalVolume* aLV)
{
  return std::find(fLCVvec.begin(), fLCVvec.end(), aLV) != fLCVvec.end();
}

// --------------------------------------------------------------------

const G4CrystalExtension* G4LogicalCrystalVolume::GetCrystal() const
{
  return dynamic_cast<G4CrystalExtension*>(dynamic_cast<G4ExtendedMaterial*>(GetMaterial())
                                ->RetrieveExtension("crystal"));
}				

// --------------------------------------------------------------------

const G4ThreeVector& G4LogicalCrystalVolume::GetBasis(G4int i) const
{
  return GetCrystal()->GetUnitCell()->GetBasis(i);
}

// --------------------------------------------------------------------

void G4LogicalCrystalVolume::SetMillerOrientation(G4int h,
                                                  G4int k,
                                                  G4int l,
                                                  G4double rot)
{
  // Align Miller normal vector (hkl) with +Z axis, and rotation about axis

  if (verboseLevel)
  {
    G4cout << "G4LatticePhysical::SetMillerOrientation(" << h << " "
           << k << " " << l << ", " << rot/CLHEP::deg << " deg)" << G4endl;
  }
    
   hMiller = h;
   kMiller = k;
   lMiller = l;
   fRot = rot;
    
   G4ThreeVector norm = (h*GetBasis(0)+k*GetBasis(1)+l*GetBasis(2)).unit();
    
   if (verboseLevel>1) G4cout << " norm = " << norm << G4endl;
    
   // Aligns geometry +Z axis with lattice (hkl) normal
   fOrient = G4RotationMatrix::IDENTITY;
   fOrient.rotateZ(rot).rotateY(norm.theta()).rotateZ(norm.phi());
   fInverse = fOrient.inverse();
    
   if (verboseLevel>1) G4cout << " fOrient = " << fOrient << G4endl;
    
   // FIXME:  Is this equivalent to (phi,theta,rot) Euler angles???
}

// --------------------------------------------------------------------

// Rotate input vector between lattice and solid orientations

const G4ThreeVector&
G4LogicalCrystalVolume::RotateToLattice(G4ThreeVector& dir) const
{
  return dir.transform(fOrient);
}

const G4ThreeVector&
G4LogicalCrystalVolume::RotateToSolid(G4ThreeVector& dir) const
{
  return dir.transform(fInverse);
}

// --------------------------------------------------------------------
