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
// $Id: checkParticles.cc,v 1.4 2003-06-16 16:58:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "G4ios.hh"
#include "globals.hh"
#include "tstParticleConstructor.hh"
#include <fstream>
#include <iomanip>

void CheckMass(const char*);
void CheckWidth(const char*);

int main(int argc,char** argv) {
  // PDG computer-readble file for particle properties
  //  see http://www-pdg.lbl.gov/computer_read.html
  std::ifstream pdgFile;
  G4String pdgFileName = "garren_98.mc";
  if (argc > 1) pdgFileName =  argv[1];
 
  // open the pdg file
  pdgFile.open((const char*)pdgFileName);
  if( pdgFile.fail() ) {
    G4cerr << "pdg data file <" << pdgFileName << "> could not open." << G4endl;
	exit(1);
  }

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

#include <strstream>

#include "G4ParticleTable.hh"

void CheckMass(const char* inputString)
{
  const char* t = inputString;
  std::istrstream is((char*)t);
 
  G4int     encoding;
  G4double  mass;
  char      massRange[20];
  char      name[20];

  is >> encoding;

  t = &inputString[33];
  std::istrstream is2((char*)t);
  is2 >> mass >> massRange >> name;

  // get a pointer to G4ParticleDefinition with encoding
  G4ParticleDefinition* particle = 
    G4ParticleTable::GetParticleTable()->FindParticle(encoding);
  if (particle == 0) {
    G4cout << "Can not find the particle with " << encoding;
    G4cout << " for " << name << G4endl;
    return;
  } 
 
  if ( abs(particle->GetPDGMass()/GeV-mass)/mass > 1.0e-4 ) {
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
    if ( abs(particle->GetPDGMass()/GeV-mass)/mass > 1.0e-3 ) {
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
  std::istrstream is((char*)t);
 
  G4int     encoding;
  G4double  width;
  char      widthRange[20];
  char      name[20];

  is >> encoding;

   t = &inputString[33];
   std::istrstream is2((char*)t);
   is2 >> width >> widthRange >> name;

  // get a pointer to G4ParticleDefinition with encoding
  G4ParticleDefinition* particle = 
    G4ParticleTable::GetParticleTable()->FindParticle(encoding);
  if (particle == 0) {
    G4cout << "Can not find the particle with " << encoding;
    G4cout << " for " << name << G4endl;
    return;
  } 
 
  if ( abs(particle->GetPDGWidth()/GeV-width)/width > 1.0e-3 ) {
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
   if ( abs(particle->GetPDGWidth()/GeV-width)/width > 1.0e-4 ) {
	G4cout << "Width is inconsistent for " << particle->GetParticleName();
    G4cout << " width = " << particle->GetPDGWidth()/GeV << G4endl;
    G4cout <<  encoding << " " << width <<  " ";
	G4cout << widthRange <<  " " << name << G4endl;
  }
 }

} 

