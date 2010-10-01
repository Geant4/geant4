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
// $Id: checkParticles.cc,v 1.13 2010-10-01 02:41:44 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "G4ios.hh"
#include "globals.hh"
#include "tstParticleConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4StateManager.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

void CheckMass(const char*);
void CheckWidth(const char*);

const double ALLOWANCE = 1.0e-5;

int main(int argc,char** argv) {
  // PDG computer-readble file for particle properties
  //  see http://www-pdg.lbl.gov/computer_read.html
  std::ifstream pdgFile;
  G4String pdgFileName = "mass_width_2010.mcd";
  if (argc > 1) pdgFileName =  argv[1];
 
  // open the pdg file
  pdgFile.open((const char*)pdgFileName);
  if( pdgFile.fail() ) {
    G4cerr << "pdg data file <" << pdgFileName << "> could not open." << G4endl;
	exit(1);
  }
  // set the initial application state
  G4StateManager::GetStateManager()->SetNewState(G4State_PreInit);


  // set Readiness
  G4ParticleTable::GetParticleTable()->SetReadiness();

  // create all particles
  tstParticleConstructor pConstructor;
  pConstructor.ConstructParticle();

  // read mass and width from pdg file and check with Geant4
  char aLine[256];
  const int LineLength = 255;

  while(1){
	// read one line
    pdgFile.getline( aLine, LineLength );
    if (pdgFile.eof() ) {
		break;
    } else if( pdgFile.bad() ) {
        G4cout << "Cannot read " << pdgFileName << "." << G4endl;
        break;
    } 

    aLine[LineLength] = '\0';
    if ( aLine[0] == 'M' ){
      CheckMass( &aLine[1] );
	} else if ( aLine[0] == 'W' ){
      CheckWidth( &aLine[1] );
	}
  }
  return EXIT_SUCCESS;
}

#include "G4ParticleTable.hh"

void CheckMass(const char* inputString)
{
  const char* t = inputString;
  std::istringstream is((char*)t);
 
  G4int     encoding;
  G4double  mass;
  char      massRange[20];
  char      name[20];

  is >> encoding;

  t = &inputString[33];
  std::istringstream is2((char*)t);
  is2 >> mass >> massRange >> name;

  // get a pointer to G4ParticleDefinition with encoding
  G4ParticleDefinition* particle = 
    G4ParticleTable::GetParticleTable()->FindParticle(encoding);
  if (particle == 0) {
    G4cout << "Can not find the particle with " << encoding;
    G4cout << " for " << name << G4endl;
    return;
  } 
 
  if ( std::abs(particle->GetPDGMass()/GeV-mass)/mass > ALLOWANCE ) {
	G4cout << "Mass is inconsistent for " << particle->GetParticleName();
    G4cout << " mass = " << particle->GetPDGMass()/GeV << G4endl;
    G4cout <<  encoding << " " << mass <<  " ";
	G4cout << massRange <<  " " << name << G4endl;
  }

  // anti_particle
  G4int anti_encoding = particle->GetAntiPDGEncoding();
  if ( (anti_encoding !=0) && (encoding != anti_encoding) ) {
	particle =
      G4ParticleTable::GetParticleTable()->FindParticle(anti_encoding);
    if (particle == 0) {
      G4cout << "Can not find the anti_particle with " << anti_encoding;
      G4cout << " for " << name << G4endl;
      return;
    } 
    if ( std::abs(particle->GetPDGMass()/GeV-mass)/mass > ALLOWANCE ) {
	  G4cout << "Mass is inconsistent for " << particle->GetParticleName();
      G4cout << " mass = " << particle->GetPDGMass()/GeV << G4endl;
      G4cout <<  encoding << " " << mass <<  " ";
	  G4cout << massRange <<  " " << name << G4endl;
    }
  }
} 

void CheckWidth(const char* inputString)
{
  const char* t = inputString;
  std::istringstream is((char*)t);
 
  G4int     encoding;
  G4double  width;
  char      widthRange[20];
  char      name[20];

  is >> encoding;

   t = &inputString[33];
   std::istringstream is2((char*)t);
   is2 >> width >> widthRange >> name;

  // get a pointer to G4ParticleDefinition with encoding
  G4ParticleDefinition* particle = 
    G4ParticleTable::GetParticleTable()->FindParticle(encoding);
  if (particle == 0) {
    G4cout << "Can not find the particle with " << encoding;
    G4cout << " for " << name << G4endl;
    return;
  } 
 
  if ( std::abs(particle->GetPDGWidth()/GeV-width)/width > ALLOWANCE ) {
	G4cout << "Width is inconsistent for " << particle->GetParticleName();
    G4cout << " width = " << particle->GetPDGWidth()/GeV << G4endl;
    G4cout <<  encoding << " " << width <<  " ";
	G4cout << widthRange <<  " " << name << G4endl;
  }

  // anti_particle
  G4int anti_encoding = particle->GetAntiPDGEncoding();
  if ((anti_encoding !=0) && (encoding != anti_encoding)) {
	particle =
      G4ParticleTable::GetParticleTable()->FindParticle(anti_encoding);
    if (particle == 0) {
      G4cout << "Can not find the anti_particle with " << anti_encoding;
      G4cout << " for " << name << G4endl;
      return;
    } 
   if ( std::abs(particle->GetPDGWidth()/GeV-width)/width > ALLOWANCE ) {
	G4cout << "Width is inconsistent for " << particle->GetParticleName();
	G4cout << " width = " << particle->GetPDGWidth()/GeV << G4endl;
	G4cout <<  encoding << " " << width <<  " ";
	G4cout << widthRange <<  " " << name << G4endl;
  }
 }

} 

