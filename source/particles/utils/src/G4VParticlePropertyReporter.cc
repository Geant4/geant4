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
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VParticlePropertyReporter.cc,v 1.1 2003/09/21 19:38:50 kurasige Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// 
// ---------------------------------------------------------------
#include "G4VParticlePropertyReporter.hh"

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
    while( (*theParticleIterator)() ){
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








