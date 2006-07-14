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
// $Id: A01Transport.cc,v 1.1 2006-07-14 14:43:36 asaim Exp $
// --------------------------------------------------------------
//
// 22-Nov-2004 Construt ALL Particles by T. Koi


#include "A01Transport.hh"

#include "globals.hh"


A01Transport::A01Transport(const G4String& name)
                     :  G4VPhysicsConstructor(name)
{
}

A01Transport::~A01Transport()
{
}

void A01Transport::ConstructParticle()
{
}

#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

void A01Transport::ConstructProcess()
{
  ////G4Transportation* theTransportationProcess= new G4Transportation();
  G4CoupledTransportation* theTransportationProcess= new G4CoupledTransportation();

  // loop over all particles in G4ParticleTable
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (!particle->IsShortLived()) {
      // Add transportation process for all particles other than  "shortlived"
      if ( pmanager == 0) {
        // Error !! no process manager
        G4String particleName = particle->GetParticleName();
        G4Exception("G4VUserPhysicsList::AddTransportation","No process manager",
                    RunMustBeAborted, particleName );
      } else {
        // add transportation with ordering = ( -1, "first", "first" )
        pmanager ->AddProcess(theTransportationProcess);
        pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
        pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
      }
    } else {
      // shortlived particle case
    }
  }

}


