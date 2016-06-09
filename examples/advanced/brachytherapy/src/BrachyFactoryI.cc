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
// Code developed by:
//  S.Guatelli
//
//    *******************************
//    *                             *
//    *    BrachyFactoryI.cc
//    *                             *
//    *******************************
//
// $Id: BrachyFactoryI.cc,v 1.4 2003/05/22 17:20:43 guatelli Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
//
#include "BrachyFactoryI.hh"
#include "BrachyPrimaryGeneratorActionI.hh"
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

G4VUserPrimaryGeneratorAction*  BrachyFactoryI::CreatePrimaryGeneratorAction()
{ 
  G4VUserPrimaryGeneratorAction* iodiumPrimaryParticle = 
                                          new BrachyPrimaryGeneratorActionI();
  return iodiumPrimaryParticle;
}

void BrachyFactoryI::CreateSource(G4VPhysicalVolume* mother)
{
  iodiumSource->ConstructIodium(mother);
}

void BrachyFactoryI::CleanSource()
{;}
