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
// $Id: G4VPVParameterisation.cc,v 1.4 2003/11/02 14:01:23 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// Default implementations for Parameterisations that do not
// parameterise solid and/or material.
// --------------------------------------------------------------------

#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

G4VSolid*
G4VPVParameterisation::ComputeSolid(const G4int,
                                    G4VPhysicalVolume *pPhysicalVol) 
{
  return pPhysicalVol->GetLogicalVolume()->GetSolid();
}
       
G4Material*
G4VPVParameterisation::ComputeMaterial(const G4int,
                                       G4VPhysicalVolume *pPhysicalVol) 
{
  return pPhysicalVol->GetLogicalVolume()->GetMaterial();
}
