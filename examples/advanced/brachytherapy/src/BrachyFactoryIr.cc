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
//    *    BrachyFactoryIr          *
//    *                             *
//    *******************************
//
// $Id: BrachyFactoryIr.cc,v 1.6 2006-06-29 15:48:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "BrachyFactoryIr.hh"
#include "BrachyPrimaryGeneratorActionIr.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"  
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh" 
#include "BrachyDetectorMessenger.hh"
#include "BrachyDetectorConstructionIr.hh"

BrachyFactoryIr:: BrachyFactoryIr()
{
  iridiumSource = new  BrachyDetectorConstructionIr(); 
  iridiumPrimaryParticle = new BrachyPrimaryGeneratorActionIr();
}

BrachyFactoryIr:: ~BrachyFactoryIr()
{
  delete iridiumSource;
  delete iridiumPrimaryParticle;
}
 
void BrachyFactoryIr::CreatePrimaryGeneratorAction(G4Event* anEvent)
{
  iridiumPrimaryParticle -> GeneratePrimaries(anEvent);
}

void BrachyFactoryIr::CreateSource(G4VPhysicalVolume* mother)
{
  iridiumSource -> ConstructIridium(mother);
}

void BrachyFactoryIr::CleanSource()
{
  iridiumSource -> CleanIridium();
}
