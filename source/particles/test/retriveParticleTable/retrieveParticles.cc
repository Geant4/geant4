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
// $Id: retrieveParticles.cc,v 1.1 2004-03-11 09:52:18 kurasige Exp $
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
#include <fstream>
#include <iomanip>

int main(int ,char** ) {
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




