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
// $Id: RE01PhysicsList.cc,v 1.1 2004/11/26 07:37:42 asaim Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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



