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
// $Id: RE01PhysicsList.cc,v 1.2 2006-06-29 17:44:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "RE01PhysicsList.hh"

#include "RE01DecayPhysics.hh"
#include "RE01BosonPhysics.hh"
#include "RE01LeptonPhysics.hh"
#include "RE01HadronPhysics.hh"
#include "RE01IonPhysics.hh"

RE01PhysicsList::RE01PhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  // SetVerboseLevel(1);

  // Particle decays
  RegisterPhysics( new RE01DecayPhysics("decay"));

  // Bosons (gamma + geantinos)
  RegisterPhysics( new RE01BosonPhysics("boson"));

  // Leptons
  RegisterPhysics( new RE01LeptonPhysics("lepton"));

  // Hadron Physics
  RegisterPhysics( new RE01HadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics( new RE01IonPhysics("ion"));
}

RE01PhysicsList::~RE01PhysicsList()
{;}

void RE01PhysicsList::SetCuts()
{
  // Use default cut values gamma and e processes
  SetCutsWithDefault();   
}



