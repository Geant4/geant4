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
// Author: S. Guatelli,University of Wollongong, Australia susanna@uow.edu.au
//
#include "G4HumanPhantomActionInitialization.hh"
#include "G4HumanPhantomPrimaryGeneratorAction.hh"
#include "G4HumanPhantomSteppingAction.hh"
#include "G4RunManager.hh"
#include "G4HumanPhantomRunAction.hh"
#include "G4HumanPhantomEventAction.hh"
#include "G4GeneralParticleSource.hh"

G4HumanPhantomActionInitialization::G4HumanPhantomActionInitialization():
G4VUserActionInitialization()
{
}


G4HumanPhantomActionInitialization::~G4HumanPhantomActionInitialization()
{
}

void G4HumanPhantomActionInitialization::BuildForMaster() const
{
	
}


void G4HumanPhantomActionInitialization::Build() const
{   
// Instantiate the analysis manager
G4HumanPhantomAnalysisManager* analysis = new G4HumanPhantomAnalysisManager();
 
SetUserAction(new G4HumanPhantomPrimaryGeneratorAction);


SetUserAction(new G4HumanPhantomRunAction(analysis));

  
G4HumanPhantomEventAction* eventAction = new G4HumanPhantomEventAction();
SetUserAction(eventAction);

SetUserAction(new G4HumanPhantomSteppingAction()); 
	
}  

