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
// $Id: BrachyFactoryI.cc,v 1.2 2002-11-18 15:18:38 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "BrachyFactoryI.hh"
#include"BrachyPrimaryGeneratorActionI.hh"
#include"BrachyDetectorConstructionI.hh"
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

BrachyFactoryI:: BrachyFactoryI()
{
   pIodio=new  BrachyDetectorConstructionI();
}

BrachyFactoryI:: ~BrachyFactoryI()
{
 
delete pIodio;
}
 

G4VUserPrimaryGeneratorAction*  BrachyFactoryI::CreatePrimaryGeneratorAction()

{ 
  
   
  G4VUserPrimaryGeneratorAction*    pIodium =new BrachyPrimaryGeneratorActionI();
 if(pIodium) return pIodium ;


 

}

void BrachyFactoryI::CreateSource(G4VPhysicalVolume* mother)
{

  pIodio -> ConstructIodium(mother);

}

void BrachyFactoryI::CleanSource()
{;}

