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
// $Id: G4PPVReplica.cc,v 1.3 2001/07/11 10:02:20 gunter Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
//
//

#include "G4PPVReplica.hh"

#include "HepODBMS/clustering/HepClustering.h"

#include "globals.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4PVReplica.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PLogicalVolume.hh"

G4PPVReplica::G4PPVReplica( G4VPhysicalVolume* PhysVol,
		                    HepRef(G4PLogicalVolume) persLogVol)
 :  G4PVPhysicalVolume(PhysVol, persLogVol)
{
  EAxis axis;
  G4int nReplicas;
  G4double width;
  G4double offset;
  G4bool consuming;

  PhysVol->GetReplicationData(axis, nReplicas, width, offset, consuming);
  faxis      = axis;
  fnReplicas = nReplicas;
  fwidth     = width;
  foffset    = offset;
  fcopyNo    = PhysVol->GetCopyNo();
}

G4PPVReplica::~G4PPVReplica()
{
}

G4VPhysicalVolume* G4PPVReplica::MakeTransientObject(
                             G4LogicalVolume* aLogical,
                             G4VPhysicalVolume* aMother)
{
  const G4String pName = (const char *) GetName();
  EAxis pAxis = faxis;
  G4int nReplicas = fnReplicas;
  G4double width = fwidth;
  G4double offset = foffset;

  G4VPhysicalVolume* aPhysVol = new G4PVReplica(
        pName, aLogical, aMother, pAxis, nReplicas, width, offset);

  return aPhysVol;
}

G4bool G4PPVReplica::IsMany() const
{
  return false; 
}

G4int G4PPVReplica::GetCopyNo() const
{
  return fcopyNo;
}

void  G4PPVReplica::SetCopyNo(G4int newCopyNo)
{
  fcopyNo= newCopyNo;
}
