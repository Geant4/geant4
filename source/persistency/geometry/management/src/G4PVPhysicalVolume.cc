//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PVPhysicalVolume.cc,v 1.5 2001/07/11 10:02:20 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// 
//                         Takashi.Sasaki@kek.jp

#include "G4PVPhysicalVolume.hh"

#include "globals.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PLogicalVolume.hh"

G4PVPhysicalVolume::G4PVPhysicalVolume()
: flogical(0), fname(0), fmother(0)
{
  frot.resize(9);
  ftrans.resize(3);
}

G4PVPhysicalVolume::G4PVPhysicalVolume(
               const G4VPhysicalVolume *PhysVol,
               const HepRef(G4PLogicalVolume) persLogVol)
{
    const G4RotationMatrix* aRot = PhysVol->GetRotation();
    if ( aRot )
      SetRotation( PhysVol->GetRotation() );

    const G4ThreeVector atrans = PhysVol->GetTranslation();
    SetTranslation( atrans );

    flogical = persLogVol;
    fname = PhysVol->GetName();	
//    fmother = LookUp(PhysVol->GetMother());
}

G4PVPhysicalVolume::~G4PVPhysicalVolume()
{
}

void G4PVPhysicalVolume::SetMother(const HepRef(G4PVPhysicalVolume) persMother)
{
  fmother = persMother;
}

G4VPhysicalVolume* G4PVPhysicalVolume::MakeTransientObject(
                       G4LogicalVolume* aLogical,
                       G4VPhysicalVolume* aMother )
{
  G4RotationMatrix* pRot = GetRotation();
  G4ThreeVector tlate = GetTranslation();
  G4String pName;

  // Need to call concrete Phys Vol class according to the physics
  // volume type.  For now, a simple physical volume is assumed for
  // quick solution.

  G4VPhysicalVolume* aPhysVol = new G4PVPlacement(
                     pRot, tlate, pName = fname, aLogical, 0, false, 0);
  aPhysVol->SetMother(aMother);

  return aPhysVol;
}

G4PString G4PVPhysicalVolume::GetName()
{
  return fname;
}

void G4PVPhysicalVolume::SetName(const G4PString aName)
{
  fname = aName;
}

HepRef(G4PLogicalVolume) G4PVPhysicalVolume::GetLogicalVolume()
{
  return flogical;
}

G4RotationMatrix* G4PVPhysicalVolume::GetRotation()
{
  // create unit matrix with default constructor
  G4RotationMatrix* aRot = new G4RotationMatrix();

  // setup the transient rotation matrix
  aRot->rotateAxes( Hep3Vector(frot[0],frot[3],frot[6]),
                    Hep3Vector(frot[1],frot[4],frot[7]),
                    Hep3Vector(frot[2],frot[5],frot[8]) );
  return aRot;
}

G4ThreeVector& G4PVPhysicalVolume::GetTranslation()
{
  G4ThreeVector atrans(ftrans[0],ftrans[1],ftrans[2]);  // local variable
  ftransvector = G4ThreeVector(atrans);       // cast by copy constructor
  return ftransvector;
}

void G4PVPhysicalVolume::SetRotation(const G4RotationMatrix* aRot)
{
  frot.resize(9);

  frot[0] = aRot->xx();
  frot[1] = aRot->xy();
  frot[2] = aRot->xz();
  frot[3] = aRot->yx();
  frot[4] = aRot->yy();
  frot[5] = aRot->yz();
  frot[6] = aRot->zx();
  frot[7] = aRot->zy();
  frot[8] = aRot->zz();

  frot.update();
}

void G4PVPhysicalVolume::SetTranslation(G4ThreeVector atrans)
{
  ftrans.resize(3);

  ftrans[0] = atrans.x();
  ftrans[1] = atrans.y();
  ftrans[2] = atrans.z();

  ftrans.update();
}

