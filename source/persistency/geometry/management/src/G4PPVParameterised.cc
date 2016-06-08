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
// $Id: G4PPVParameterised.cc,v 1.3.2.1 2001/06/28 19:11:28 gunter Exp $
// GEANT4 tag $Name:  $
//
//
//

#include "G4PPVParameterised.hh"

#include "HepODBMS/clustering/HepClustering.h"

#include "globals.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4PVParameterised.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PLogicalVolume.hh"

G4PPVParameterised::G4PPVParameterised( G4VPhysicalVolume* PhysVol,
		                    HepRef(G4PLogicalVolume) persLogVol)
 :  G4PPVReplica(PhysVol, persLogVol)
{
  // G4VPVParameterisation class has no data member
  //  fparam = 0;
}

G4PPVParameterised::~G4PPVParameterised()
{
}

G4VPhysicalVolume* G4PPVParameterised::MakeTransientObject(
                             G4LogicalVolume* aLogical,
                             G4VPhysicalVolume* aMother)
{
  const G4String pName = (const char *) GetName();
  EAxis pAxis = faxis;
  G4int nReplicas = fnReplicas;

  G4VPhysicalVolume* aPhysVol = new G4PVParameterised(
       pName, aLogical, aMother, pAxis, nReplicas, 0);

  return aPhysVol;
}
