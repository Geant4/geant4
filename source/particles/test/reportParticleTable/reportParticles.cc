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
// $Id: reportParticles.cc,v 1.5 2002-03-26 07:34:42 kurasige Exp $
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
#include "g4std/fstream"
#include "g4std/iomanip"


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
  html.Print( *quarkContainer,  "quarks  27/Mar/02");
  html.Print( *leptonContainer, "leptons 27/Mar/02");
  html.Print( *mesonContainer,  "mesons 27/Mar/02");
  html.Print( *baryonContainer, "baryons 27/Mar/02");
  html.Print( *ionContainer,    "ions 27/Mar/02");
  html.Print( *otherContainer,  "others 27/Mar/02");
  return EXIT_SUCCESS;
}




