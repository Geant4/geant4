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
//    *    BrachyFactoryIr          *
//    *                             *
//    *******************************
//
// $Id: BrachyFactoryIr.cc,v 1.2 2002-11-18 15:18:38 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "BrachyFactoryIr.hh"
#include"BrachyPrimaryGeneratorActionIr.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"  
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4RadioactiveDecay.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include <math.h>
#include "G4RunManager.hh" 
#include"BrachyDetectorMessenger.hh"
#include"BrachyDetectorConstructionIr.hh"
BrachyFactoryIr:: BrachyFactoryIr()
{
   pIridio=new  BrachyDetectorConstructionIr();

}

BrachyFactoryIr:: ~BrachyFactoryIr()
{
  delete pIridio;

}
 

G4VUserPrimaryGeneratorAction*  BrachyFactoryIr::CreatePrimaryGeneratorAction()

{ 
 
   
  G4VUserPrimaryGeneratorAction*    pIridium =new BrachyPrimaryGeneratorActionIr();
 if(pIridium) return pIridium ;


}

void BrachyFactoryIr::CreateSource(G4VPhysicalVolume* mother)
{

  pIridio -> ConstructIridium(mother);

}
void BrachyFactoryIr::CleanSource()
{

  pIridio -> CleanIridium();

}
