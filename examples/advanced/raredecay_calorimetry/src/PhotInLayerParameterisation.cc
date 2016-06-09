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
// $Id: PhotInLayerParameterisation.cc,v 1.3 2006/06/29 16:25:17 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
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
  
