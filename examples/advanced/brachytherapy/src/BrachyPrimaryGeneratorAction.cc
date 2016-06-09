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
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed by:
// S.Guatelli
//
//    ********************************************
//    *                                          *
//    *    BrachyPrimaryGeneratorAction.cc       *
//    *                                          *
//    ********************************************
//
// $Id$
//
#include "globals.hh"
#include "BrachyPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"  
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include "BrachyFactory.hh"
#include "BrachyFactoryLeipzig.hh"
#include "BrachyFactoryIr.hh"
#include "BrachyFactoryI.hh"
#include "BrachyPrimaryGeneratorMessenger.hh"

BrachyPrimaryGeneratorAction::BrachyPrimaryGeneratorAction()
{

 primaryMessenger = new BrachyPrimaryGeneratorMessenger(this);
 // Default source: iridium source 
 factory = new BrachyFactoryIr();
}

BrachyPrimaryGeneratorAction::~BrachyPrimaryGeneratorAction()
{
 delete factory;
 delete primaryMessenger;
}

void BrachyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  factory -> CreatePrimaryGeneratorAction(anEvent);
}

void BrachyPrimaryGeneratorAction::SwitchEnergy(G4String sourceChoice)
{
  G4int flag = 0;

  // Switch the energy spectrum of the photons delivered by the radiative source	
  if (sourceChoice == "Iodium")
    {
      flag=1;
      if (factory) delete factory;
    }
  switch(flag)
    {
    case 1:
      factory = new BrachyFactoryI;
      break;
    default:   
      factory = new BrachyFactoryIr; 
    }      
}
