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
// $Id: RemSimPrimaryGeneratorAction.cc,v 1.8 2004/05/27 10:33:12 guatelli Exp $// Author: Susanna Guatelli, guatelli@ge.infn.it

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
  primaryFactory1 = new RemSimBasicGenerator();
  primaryFactory2 = new RemSimInterplanetarySpaceConfiguration();
  messenger = new RemSimPrimaryGeneratorMessenger(this);
}

RemSimPrimaryGeneratorAction::~RemSimPrimaryGeneratorAction()
{
  delete messenger;
  delete primaryFactory2;
  delete primaryFactory1;
}

G4double RemSimPrimaryGeneratorAction::GetInitialEnergy()
{
  G4double initialEnergy = 0.;

  if(value == "Basic") 
      initialEnergy = primaryFactory1 -> GetInitialEnergy();   
    
  else if(value == "Interplanetary")  
      initialEnergy = primaryFactory2 -> GetInitialEnergy();    

  return initialEnergy;
}

void RemSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(value == "Basic") 
       primaryFactory1 -> GeneratePrimaries(anEvent);
   
  else if(value == "Interplanetary") 
    { 
      primaryFactory2 ->  SetMoon(false);
      primaryFactory2 -> GeneratePrimaries(anEvent);    
     }
  else if(value == "Moon")        
    {
      primaryFactory2 ->  SetMoon(true);
      primaryFactory2 -> GeneratePrimaries(anEvent);
    }
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
    }

  else if (value == "Moon") 
    {
      G4cout<< "The configuration is the Moon Configuration" <<G4endl;
      G4cout<< 
	"Remember to type /run/data file.txt with the energy spectrum of primary particles!!!" <<G4endl;
    }

  else G4cout << "This Generator is not defined!" <<G4endl;  
}
