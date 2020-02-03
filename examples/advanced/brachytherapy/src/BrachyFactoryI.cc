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
// Code developed by:
//  S.Guatelli
//
//    *******************************
//    *                             *
//    *    BrachyFactoryI.cc
//    *                             *
//    *******************************
//
//
#include "BrachyFactoryI.hh"
#include "BrachyDetectorConstructionI.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"  
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include "G4RunManager.hh" 

BrachyFactoryI:: BrachyFactoryI()
{
  iodiumSource = new BrachyDetectorConstructionI();
}

BrachyFactoryI::~BrachyFactoryI()
{
  delete iodiumSource;
}

void BrachyFactoryI::CreateSource(G4VPhysicalVolume* mother)
{
  iodiumSource -> ConstructIodium(mother);
}

void BrachyFactoryI::CleanSource()
{
 iodiumSource -> CleanIodium(); 
 iodiumSource = 0;
}
