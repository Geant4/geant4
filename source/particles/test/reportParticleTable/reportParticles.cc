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
// $Id: reportParticles.cc,v 1.8 2003-09-21 19:38:03 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "G4ios.hh"
#include "globals.hh"
#include "tst2ParticleConstructor.hh"
#include "G4ParticlePropertyTable.hh"
#include "G4SimplePPReporter.hh"
#include "G4HtmlPPReporter.hh"
#include <fstream>
#include <iomanip>

#define USE_HTML_REPORTER 1

int main(int argc,char** argv) {
  // create all particles
  tst2ParticleConstructor pConstructor;
  pConstructor.ConstructParticle();

  // particleContainer 
  G4VParticlePropertyReporter* aPPR;
#ifdef USE_HTML_REPORTER
  G4VParticlePropertyReporter* pLepton = new G4HtmlPPReporter();
  G4VParticlePropertyReporter* pMeson  = new G4HtmlPPReporter();
  G4VParticlePropertyReporter* pBaryon = new G4HtmlPPReporter();
  G4VParticlePropertyReporter* pIon    = new G4HtmlPPReporter();
  G4VParticlePropertyReporter* pQuark  = new G4HtmlPPReporter();
  G4VParticlePropertyReporter* pOther  = new G4HtmlPPReporter();
#else
  G4VParticlePropertyReporter* pLepton = new G4SimplePPReporter();
  G4VParticlePropertyReporter* pMeson  = new G4SimplePPReporter();
  G4VParticlePropertyReporter* pBaryon = new G4SimplePPReporter();
  G4VParticlePropertyReporter* pIon    = new G4SimplePPReporter();
  G4VParticlePropertyReporter* pQuark  = new G4SimplePPReporter();
  G4VParticlePropertyReporter* pOther  = new G4SimplePPReporter();
#endif
  // pointer to the particle table
  G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  G4ParticleTable::G4PTblDicIterator* theParticleIterator;
  theParticleIterator = theParticleTable->GetIterator();

  G4ParticlePropertyTable* thePPTable = G4ParticlePropertyTable::GetParticlePropertyTable();

  // loop over all particles in G4ParticleTable 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String type = particle->GetParticleType();
    if  (type=="lepton") {
      aPPR = pLepton;
    } else  if  (type=="meson") {
      aPPR = pMeson;

    } else  if  (type=="baryon") {
      aPPR = pBaryon;

    } else  if  (type=="nucleus") {
      aPPR = pIon;

    } else  if  ((type=="quarks")||(type=="diquarks")||(type=="gluons")) {
      aPPR = pQuark;

    } else    {
      aPPR = pOther;
    } 
    // add property data to PPDcontainer
    aPPR->FillList(particle->GetParticleName()); 
}
  pLepton->Print("leptons");
  pMeson->Print("mesons");
  pBaryon->Print("baryons");
  pIon->Print("ions");
  pQuark->Print("quarks");
  pOther->Print("others");

  return EXIT_SUCCESS;
}




