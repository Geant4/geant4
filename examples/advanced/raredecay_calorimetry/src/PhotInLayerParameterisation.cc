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
// $Id: PhotInLayerParameterisation.cc,v 1.2 2005-12-09 16:44:21 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//#define debug

#include "PhotInLayerParameterisation.hh"

PhotInLayerParameterisation::PhotInLayerParameterisation()
:numberOfLayers(PhotInNOfLayers),absMaterial(0),hzTot(1.) {}

PhotInLayerParameterisation::~PhotInLayerParameterisation() {}

void PhotInLayerParameterisation::ComputeTransformation(const G4int copyNo,
                                                        G4VPhysicalVolume* physVol) const
{
  G4double halfZ = hzTot/numberOfLayers; // half width of one layer
  G4double Zposition= (copyNo+copyNo)*halfZ - (hzTot-halfZ); // Z-positions of hte layers
#ifdef debug
  G4cout<<"PhotInLayerParameterisation::ComputeTransformation: hZ="<<halfZ
        <<", Z="<<Zposition<<G4endl;
#endif
  G4ThreeVector origin(0,0,Zposition);   // 3-vector for positioning of the layers
  physVol->SetTranslation(origin);       // positioning of the layer
  physVol->SetRotation(0);               // without rotation (can be turned!)
  physVol->SetName(PhotInAbsorberName);  // Name of the layer which is primarily absorber
}

G4Material* PhotInLayerParameterisation::
                     ComputeMaterial(const G4int, G4VPhysicalVolume*, const G4VTouchable*)
{
  return absMaterial;                    // MaterialOfAbsorber in diff. layers can be diff.
}
  
