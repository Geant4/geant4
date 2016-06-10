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
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VParticlePropertyReporter.cc 91885 2015-08-10 07:05:56Z gcosmo $
//
// 
// ---------------------------------------------------------------
#include "G4VParticlePropertyReporter.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

////////////////////////////
G4VParticlePropertyReporter::G4VParticlePropertyReporter()
{
  pPropertyTable =   G4ParticlePropertyTable::GetParticlePropertyTable();
}

/////////////////////////////
G4VParticlePropertyReporter::~G4VParticlePropertyReporter()
{
  pList.clear();
  pPropertyTable->Clear();
}    

///////////////////////////
G4bool G4VParticlePropertyReporter::FillList(G4String name)
{
  G4ParticlePropertyData* pData = pPropertyTable->GetParticleProperty(name);
  G4bool result = false;
  if (pData != 0) {
    //the particle exists
    pList.push_back(pData);
    result = true;
  } else {
    // pointer to the particle table
    G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator* theParticleIterator;
    theParticleIterator = theParticleTable->GetIterator();
    
    // loop over all particles in G4ParticleTable 
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){ // Loop checking, 09.08.2015, K.Kurashige
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4String type = particle->GetParticleType();
      pData =pPropertyTable->GetParticleProperty(particle);
      if ( name == "all" ) {
	pList.push_back(pData);
        result = true;
      } else if ( name == type ) {
	pList.push_back(pData);
        result = true;
      } 
    }
  }
  return result;
}

void G4VParticlePropertyReporter::Clear()
{
  pList.clear();
}








