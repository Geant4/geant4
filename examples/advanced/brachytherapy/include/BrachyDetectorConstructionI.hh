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
// $Id: BrachyDetectorConstructionI.hh,v 1.4 2006-05-12 13:49:29 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstructionI.hh     *
//    *                                      *
//    ****************************************
// Model of the Iodium source
//
#ifndef BrachyDetectorConstructionI_H
#define BrachyDetectorConstructionI_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class BrachyMaterial;

class BrachyDetectorConstructionI
{
public:
  BrachyDetectorConstructionI();
  ~BrachyDetectorConstructionI();
  void  ConstructIodium(G4VPhysicalVolume*);// Construct iodium source

private:
  G4VPhysicalVolume* capsulePhys;
  G4VPhysicalVolume* capsuleTipPhys1;
  G4VPhysicalVolume* capsuleTipPhys2;
  G4VPhysicalVolume* iodiumCorePhys;
  G4VPhysicalVolume* markerPhys;

  BrachyMaterial* pMaterial;   
};
#endif








