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
// $Id: BrachyDetectorConstructionLeipzig.hh,v 1.4 2006-05-12 14:57:54 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//    ************************************************
//    *                                              *
//    *      BrachyDetectorConstructionLeipzig.hh    *
//    *                                              *
//    ************************************************
//
// Management of Leipzig applicator 
//
#ifndef BrachyDetectorConstructionLeipzig_H
#define BrachyDetectorConstructionLeipzig_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class BrachyMaterial;

class BrachyDetectorConstructionLeipzig
{
public:
  BrachyDetectorConstructionLeipzig();
  ~BrachyDetectorConstructionLeipzig(); 
  void ConstructLeipzig(G4VPhysicalVolume*);

private:
  G4VPhysicalVolume* capsulePhys;
  G4VPhysicalVolume* capsuleTipPhys;
  G4VPhysicalVolume* iridiumCorePhys;
  G4VPhysicalVolume* applicator1Phys;
  G4VPhysicalVolume* applicator2Phys;

  BrachyMaterial* pMaterial;   
};

#endif
