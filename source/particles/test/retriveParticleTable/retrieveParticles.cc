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
// $Id: retrieveParticles.cc,v 1.3 2009-09-15 14:38:47 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "G4ios.hh"
#include "globals.hh"
#include "tst2ParticleConstructor.hh"
#include "G4ParticlePropertyTable.hh"
#include "G4TextPPRetriever.hh"
#include "G4SimplePPReporter.hh"
#include "G4StateManager.hh"
#include <fstream>
#include <iomanip>

int main(int ,char** ) 
{
  // set the initial application state
  G4StateManager::GetStateManager()->SetNewState(G4State_PreInit);

  // create all particles
  tst2ParticleConstructor pConstructor;
  pConstructor.ConstructParticle();

  // Retrive Particle Properties from File 
  G4TextPPRetriever fPPRetriever;
  fPPRetriever.Retrieve("leptons");

  // particleContainer 
  G4SimplePPReporter* pLepton = new G4SimplePPReporter();

  // pointer to the particle table
  G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  G4ParticleTable::G4PTblDicIterator* theParticleIterator;
  theParticleIterator = theParticleTable->GetIterator();

  // loop over all particles in G4ParticleTable 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String type = particle->GetParticleType();
    if  (type=="lepton") {
    // add property data to PPDcontainer
      pLepton->FillList(particle->GetParticleName());;
    } 
  }
  pLepton->Print();

  return EXIT_SUCCESS;
}




