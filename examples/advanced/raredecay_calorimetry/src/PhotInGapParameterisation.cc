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
// $Id: PhotInGapParameterisation.cc,v 1.5 2006/06/29 16:25:15 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//

//#define debug

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
  physVol->SetName(PhotInSlabName);     // Name of the positioned slab
}

G4Material* PhotInGapParameterisation::
                      ComputeMaterial(const G4int, G4VPhysicalVolume*, const G4VTouchable*)
{
  return gapMaterial;                   // Material can be different for different slabs
}
  
