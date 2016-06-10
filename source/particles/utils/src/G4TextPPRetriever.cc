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
// $Id: G4TextPPRetriever.cc 91885 2015-08-10 07:05:56Z gcosmo $
//
// 
// ---------------------------------------------------------------
#include "G4TextPPRetriever.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayTable.hh"  
#include "G4VDecayChannel.hh"  
#include "G4Tokenizer.hh"
#include <iomanip>
#include <fstream>       

//////////////////////////////
G4TextPPRetriever::G4TextPPRetriever():G4VParticlePropertyRetriever()
{ 
}

////////////////////////////
G4TextPPRetriever::~G4TextPPRetriever()
{
}    

/////////////////////
void G4TextPPRetriever::Retrieve(const G4String& option)
{
  SparseOption( option );
 
 // pointer to the particle table
  G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  G4ParticleTable::G4PTblDicIterator* theParticleIterator;
  theParticleIterator = theParticleTable->GetIterator();
    
  // loop over all particles in G4ParticleTable 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){ // Loop checking, 09.08.2015, K.Kurashige
    G4ParticleDefinition* particle = theParticleIterator->value();
    ModifyPropertyTable(particle);
  }
}    


void G4TextPPRetriever::SparseOption(const G4String& option)
{
  G4Tokenizer savedToken( option );
  
  // 1st option : base directory
  baseDir = savedToken();
  if (!baseDir.isNull()) {
    if(baseDir(baseDir.length()-1)!='/') {
      baseDir += "/";
    }
  }
}



G4bool  G4TextPPRetriever::ModifyPropertyTable(const G4ParticleDefinition* particle)
{
  G4String name = particle->GetParticleName();
  
  //--- open file -----
  G4String fileName = baseDir + name + ".txt";
  // exception
  if (name == "J/psi") fileName = baseDir +"jpsi.txt";

  std::ifstream inFile(fileName, std::ios::in );
  if (!inFile) return false;
  
  // GetParticleProperty
  G4ParticlePropertyData* pData = pPropertyTable->GetParticleProperty(name);

  // particle name  encoding
  G4String name_t;
  G4int    encoding;
  inFile >> name_t >> encoding;
  if ( (name != name_t) || (encoding !=  pData->GetPDGEncoding()) ){
    G4cout << "G4TextPPRetriever::ModifyPropertyTable:   ";
    G4cout << "particle name or encoding mismatch for " << name ;
    G4cout << G4endl;
    return false;
  }

  // IJPC
  G4int  iIsoSpin, iSpin, iParity, iConj;
  inFile >>  iIsoSpin >> iSpin >> iParity >> iConj;   
  if  (  ( iIsoSpin != pData->GetPDGiIsospin()) ||
         ( iSpin    != pData->GetPDGiSpin())    ||
	 ( iParity  != pData->GetPDGiParity())  ||
	 ( iConj    != pData->GetPDGiConjugation()) ){
    G4cout << "G4TextPPRetriever::ModifyPropertyTable:   ";
    G4cout << "IJPC mismatch for " << name ;
    G4cout << G4endl;
    return false;
  }

  // mass, width, charge 
  G4double mass, width, charge;
  inFile >> mass >> width >>  charge;
  mass *= GeV;
  width *= GeV;
  charge *= eplus;
  if (mass  != pData->GetPDGMass()){ pData->SetPDGMass(mass);}
  if (width != pData->GetPDGWidth()){ pData->SetPDGWidth(width);}
  if (charge != pData->GetPDGCharge()){ pData->SetPDGCharge(charge);}

  // life time
  G4double tlife;
  inFile >> tlife;
  tlife *= second;
  if (tlife  != pData->GetPDGLifeTime()){ pData->SetPDGLifeTime(tlife);}

  pPropertyTable->SetParticleProperty(*pData);

  // Decay Table  
  G4DecayTable* dcyTable = particle->GetDecayTable(); 
  if (dcyTable == 0) return true;
  
  G4int idx =0;
  while (!inFile.eof() ) { // Loop checking, 09.08.2015, K.Kurashige
    G4double br;
    G4int    n_daughters;
    inFile >> br >> n_daughters;

    G4VDecayChannel * channel = dcyTable->GetDecayChannel(idx);

    if ( n_daughters == channel->GetNumberOfDaughters()) {
      channel->SetBR(br);
    }

    idx += 1;
    if (idx>= dcyTable->entries()) break;
  }
  return true;
}





