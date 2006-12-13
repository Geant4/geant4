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
// $Id: A01Transport.cc,v 1.2 2006-12-13 15:49:28 gunter Exp $
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


