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
// $Id: RemSimPrimaryGeneratorAction.cc,v 1.11 2005/05/27 14:21:42 guatelli Exp $// Author: Susanna Guatelli, guatelli@ge.infn.it

#include "RemSimPrimaryGeneratorAction.hh"
#include "RemSimPrimaryGeneratorMessenger.hh"
#include "RemSimBasicGenerator.hh"
#include "RemSimInterplanetarySpaceConfiguration.hh"
#include "RemSimVPrimaryGeneratorFactory.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

RemSimPrimaryGeneratorAction::RemSimPrimaryGeneratorAction()
{
  value = "Basic";
  primaryFactory = new RemSimBasicGenerator();
  messenger = new RemSimPrimaryGeneratorMessenger(this);
}

RemSimPrimaryGeneratorAction::~RemSimPrimaryGeneratorAction()
{
  delete messenger;
  delete primaryFactory;
}

G4double RemSimPrimaryGeneratorAction::GetInitialEnergy()
{
  G4double initialEnergy = 0.;

  initialEnergy = primaryFactory -> GetInitialEnergy();   
    
  return initialEnergy;
}

void RemSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
 if(value == "Interplanetary") 
    { 
      primaryFactory ->  SetMoon(false);
     }
  else if(value == "Moon")        
    {
      primaryFactory ->  SetMoon(true);
    }

  primaryFactory -> GeneratePrimaries(anEvent); 
}

void RemSimPrimaryGeneratorAction::SelectPrimaries(G4String val)
{ 
  value = val;

  if(value == "Basic") 
      G4cout<< "The configuration is the basic generator" <<G4endl;

  else if(value == "Interplanetary") 
    {
      G4cout<< "The configuration is the interplanetary space configuration"
	    <<G4endl;
      G4cout<< 
	"Remember to type /run/data file.txt with the energy spectrum of primary particles!!!" <<G4endl;
      delete primaryFactory;
      primaryFactory = 0;
      primaryFactory = new RemSimInterplanetarySpaceConfiguration();
    }

  else if (value == "Moon") 
    {
      G4cout<< "The configuration is the Moon Configuration" <<G4endl;
      G4cout<< 
	"Remember to type /run/data file.txt with the energy spectrum of primary particles!!!" <<G4endl;
      delete primaryFactory;
      primaryFactory =0;
      primaryFactory = new RemSimInterplanetarySpaceConfiguration();
    }

  else G4cout << "This Generator is not defined!" <<G4endl;  
}

void RemSimPrimaryGeneratorAction::Read(G4String fileName)
{
    primaryFactory -> Read(fileName); 
}
