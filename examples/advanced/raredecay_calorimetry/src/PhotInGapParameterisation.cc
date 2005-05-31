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
// $Id: PhotInGapParameterisation.cc,v 1.2 2005-05-31 15:23:01 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#define debug

#include "PhotInGapParameterisation.hh"

PhotInGapParameterisation::PhotInGapParameterisation()
:numberOfSlabs(PhotInNOfSlabs),gapMaterial(0),hyTot(1.),hzLay(1.),sampFract(.5) {}

PhotInGapParameterisation::~PhotInGapParameterisation() {}

void PhotInGapParameterisation::ComputeTransformation(const G4int copyNo,
                                                      G4VPhysicalVolume* physVol) const
{
  G4double halfZ = hzLay*sampFract;     // half thichness of a slab of active part of layer
  G4double halfY = hyTot/numberOfSlabs; // half width of each slab after the subdivision
  G4double Zposition= halfZ - hzLay;    // Z-Position of all Active Gaps (rest is absorber)
  G4double Yposition= (copyNo+copyNo)*halfY - (hyTot-halfY); // Trans (y) position of slabs
#ifdef debug
  G4cout<<"PhotInGapParameterisation::ComputeTransformation: hZ="<<halfZ<<", hY="<<halfY
        <<", Z="<<Zposition<<", Y="<<Yposition<<G4endl;
#endif
  G4ThreeVector origin(0,Yposition,Zposition); // 3-vector for positioning of the slabs
  physVol->SetTranslation(origin);      // Positioning of the slab
  physVol->SetRotation(0);              // without rotation (which is possible!)
  physVol->SetName(PhotInSlabName);             // Name of the positioned slab              
}

G4Material* PhotInGapParameterisation::ComputeMaterial (const G4int, G4VPhysicalVolume*)
{
  return gapMaterial;                   // Material can be different for different slabs
}
  
