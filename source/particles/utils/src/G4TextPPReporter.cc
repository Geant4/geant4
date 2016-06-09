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
// $Id: G4TextPPReporter.cc,v 1.1 2004/03/11 09:47:45 kurasige Exp $
// GEANT4 tag $Name: geant4-06-01 $
//
// 
// ---------------------------------------------------------------
#include "G4TextPPReporter.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4DecayTable.hh"  
#include "G4Tokenizer.hh"
#include <iomanip>
#include <fstream>       

//////////////////////////////
G4TextPPReporter::G4TextPPReporter():G4VParticlePropertyReporter()
{ 
}

////////////////////////////
G4TextPPReporter::~G4TextPPReporter()
{
}    

/////////////////////
void G4TextPPReporter::Print(const G4String& option)
{
  SparseOption( option );

  for (size_t i=0; i< pList.size(); i++){
    G4ParticleDefinition* particle  = G4ParticleTable::GetParticleTable()->FindParticle( pList[i]->GetParticleName() ); 

    GeneratePropertyTable(particle);
  }
}    


void G4TextPPReporter::SparseOption(const G4String& option)
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



void  G4TextPPReporter::GeneratePropertyTable(G4ParticleDefinition* particle)
{
  G4String name = particle->GetParticleName();
  
  //--- open file -----
  G4String fileName = baseDir + name + ".txt";
  // exception
  if (name == "J/psi") fileName = baseDir +"jpsi.txt";

  std::ofstream outFile(fileName, std::ios::out );
  outFile.setf( std::ios:: scientific, std::ios::floatfield );
  outFile << std::setprecision(7) << G4endl;

  // particle name  encoding
  outFile << name << " "
	  <<  particle->GetPDGEncoding()  << G4endl;

  // IJPC
  outFile << particle->GetPDGiIsospin()       << " "
	  << particle->GetPDGiSpin()          << " "
	  << particle->GetPDGiParity()        << " "
	  << particle->GetPDGiConjugation()   << G4endl;

  // mass, width, charge 
  outFile <<  particle->GetPDGMass()/GeV      << " "
	  <<  particle->GetPDGWidth()/GeV     << " " 
	  <<  particle->GetPDGCharge()/eplus  << G4endl;

  // life time
  outFile << particle->GetPDGLifeTime()/second << G4endl;

// Decay Table  
  G4DecayTable* dcyTable = particle->GetDecayTable(); 
  if (dcyTable != 0) { 
    for (G4int i=0; i< dcyTable->entries(); i++){
      G4VDecayChannel * channel = dcyTable->GetDecayChannel(i);
      // column 1  : BR
      outFile << channel->GetBR() << " ";
      // column 2.  : daughters
      outFile << channel->GetNumberOfDaughters() << " ";
      // column 3 : Kinematics
      outFile << channel->GetKinematicsName() << " ";
      // daughters
      for (G4int j=0; j< channel->GetNumberOfDaughters(); j++){
        outFile << channel->GetDaughter(j)->GetParticleName() << " ";
      }
      outFile << G4endl;
    }
  }
}





