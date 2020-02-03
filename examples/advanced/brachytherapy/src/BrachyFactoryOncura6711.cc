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
//  S.Guatelli, D. Cutajar, A. Le
//
//    ********************************
//    *                              *
//    *  BrachyFactoryOncura6711.cc  *
//    *                              *
//    ********************************
//
// 
//
#include "globals.hh"
#include "BrachyFactoryOncura6711.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"  
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh" 
#include "BrachyDetectorMessenger.hh"
#include "BrachyDetectorConstructionOncura6711.hh"

BrachyFactoryOncura6711:: BrachyFactoryOncura6711()
{
  Oncura6711IodineSource = new  BrachyDetectorConstructionOncura6711(); 
}

BrachyFactoryOncura6711:: ~BrachyFactoryOncura6711()
{
  delete Oncura6711IodineSource;
}
 
void BrachyFactoryOncura6711::CreateSource(G4VPhysicalVolume* mother)
{
  Oncura6711IodineSource -> ConstructOncura6711(mother);
}

void BrachyFactoryOncura6711::CleanSource()
{
  Oncura6711IodineSource -> CleanOncura6711();
  Oncura6711IodineSource = 0;
}
