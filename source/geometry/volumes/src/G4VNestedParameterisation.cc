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
// $Id: G4VNestedParameterisation.cc,v 1.6 2005/11/24 18:19:53 japost Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// class G4VNestedParameterisation implementation
//
// --------------------------------------------------------------------

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VTouchable.hh"

#include "G4VNestedParameterisation.hh" 

G4VNestedParameterisation::G4VNestedParameterisation() :
  G4VPVParameterisation(), G4VVolumeMaterialScanner()
{
}

G4VNestedParameterisation::~G4VNestedParameterisation() {}

G4VSolid* G4VNestedParameterisation::ComputeSolid(const G4int, 
                                     G4VPhysicalVolume  *pvol)
{ 
  return pvol->GetLogicalVolume()->GetSolid(); 
}

G4bool G4VNestedParameterisation::IsNested() const 
{ 
  return true;
}

G4VVolumeMaterialScanner* G4VNestedParameterisation::GetMaterialScanner() 
{ 
  return this; 
} 

G4Material* 
G4VNestedParameterisation::ComputeMaterial(const G4int repNo, 
					   G4VPhysicalVolume *currentVol,
					   const G4VTouchable *parentTouch)
{
    return ComputeMaterial( currentVol, repNo, parentTouch );
}
