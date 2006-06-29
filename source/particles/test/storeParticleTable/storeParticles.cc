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
// $Id: storeParticles.cc,v 1.2 2006-06-29 19:27:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "G4ios.hh"
#include "globals.hh"
#include "tst2ParticleConstructor.hh"
#include "G4ParticlePropertyTable.hh"
#include "G4TextPPReporter.hh"
#include <fstream>
#include <iomanip>

int main(int argc,char** argv) {
  // create all particles
  tst2ParticleConstructor pConstructor;
  pConstructor.ConstructParticle();

  // particleContainer 
  G4VParticlePropertyReporter* aPPR;

  G4VParticlePropertyReporter* pLepton = new G4TextPPReporter();
  G4VParticlePropertyReporter* pMeson  = new G4TextPPReporter();
  G4VParticlePropertyReporter* pBaryon = new G4TextPPReporter();
  G4VParticlePropertyReporter* pIon    = new G4TextPPReporter();
  G4VParticlePropertyReporter* pQuark  = new G4TextPPReporter();
  G4VParticlePropertyReporter* pOther  = new G4TextPPReporter();

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




