// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: reportParticles.cc,v 1.1 1999-06-17 04:50:01 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "G4ios.hh"
#include "globals.hh"
#include "tst2ParticleConstructor.hh"
#include "tst2ParticleContainer.hh"
#include "tst2SimpleReporter.hh"
#include "tst2HtmlReporter.hh"
#include <fstream.h>
#include <iomanip.h>


int main(int argc,char** argv) {
  // create all particles
  tst2ParticleConstructor pConstructor;
  pConstructor.ConstructParticle();

  // particleContainer 
  tst2ParticleContainer* leptonContainer = new tst2ParticleContainer();
  tst2ParticleContainer* mesonContainer = new tst2ParticleContainer();
  tst2ParticleContainer* baryonContainer = new tst2ParticleContainer();
  tst2ParticleContainer* ionContainer = new tst2ParticleContainer();
  tst2ParticleContainer* quarkContainer = new tst2ParticleContainer();
  tst2ParticleContainer* otherContainer = new tst2ParticleContainer();

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
      leptonContainer->Insert(particle); 
    } else  if  (type=="meson") {
      mesonContainer->Insert(particle); 
    } else  if  (type=="baryon") {
      baryonContainer->Insert(particle); 
    } else  if  (type=="nucleus") {
      ionContainer->Insert(particle); 
    } else  if  ((type=="quarks")||(type=="diquarks")||(type=="gluons")) {
      quarkContainer->Insert(particle); 
    } else    {
      otherContainer->Insert(particle); 
    } 
  }

  // tst2SimpleReporter reporter;
  //reporter.Print(*quarkContainer);
  //reporter.Print(*leptonContainer);
  //reporter.Print(*mesonContainer);
  //reporter.Print(*baryonContainer);
  //reporter.Print(*ionContainer);
  //reporter.Print(*otherContainer);
  
  tst2HtmlReporter html;
  html.Print( *quarkContainer,  "quarks  17/June/99");
  html.Print( *leptonContainer, "leptons 17/June/99");
  html.Print( *mesonContainer,  "mesons 17/June/99");
  html.Print( *baryonContainer, "baryons 17/June/99");
  html.Print( *ionContainer,    "ions 17/June/99");
  html.Print( *otherContainer,  "others 17/June/99");
  return EXIT_SUCCESS;
}


