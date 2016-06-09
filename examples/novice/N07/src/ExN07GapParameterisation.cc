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
// $Id: ExN07GapParameterisation.cc,v 1.2 2003/05/28 11:24:35 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//

#include "ExN07GapParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"

ExN07GapParameterisation::ExN07GapParameterisation()
:numberOfLayers(40),absMaterial(0),gapMaterial(0)
{;}

ExN07GapParameterisation::~ExN07GapParameterisation()
{;}

void ExN07GapParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4double halfZ = 1.*m/G4double(2*numberOfLayers);
  G4double Zposition= G4double(copyNo)*2.*halfZ - (1.*m-halfZ);
  G4ThreeVector origin(0,0,Zposition);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
  if(copyNo%2==0)
  { physVol->SetName("Abs"); }
  else
  { physVol->SetName("Gap"); }
}

G4Material*
ExN07GapParameterisation::
ComputeMaterial (const G4int copyNo, G4VPhysicalVolume*)
{
  if(copyNo%2==0) { return absMaterial; }
  return gapMaterial;
}
  
