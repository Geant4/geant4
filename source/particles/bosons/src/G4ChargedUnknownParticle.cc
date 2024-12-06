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
// --------------------------------------------------------------------------
//
//      GEANT4 source file 
//
//      File name:     G4ChargedUnknownParticle.cc
//
//      Author:        A.Ribon
// 
//      Creation date: August 2024
//
//      Description:   This class is similar to G4UnknownParticle,
//                     i.e. representing particles with valid PDG code but
//                     unknown to Geant4 - e.g. produced by MC event
//                     generators - but with the extra requirement of having
//                     a non-zero electric charge.
//                     While for G4UnknownParticle is possible to assign
//                     transportation and decay processes, for
//                     G4ChargedUnknownParticle is also possible to assign
//                     ionisation and multiple scattering.
//
//      Modifications:
//      
// --------------------------------------------------------------------------
//

#include "G4ChargedUnknownParticle.hh"
#include "G4ParticleTable.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"


G4ChargedUnknownParticle* G4ChargedUnknownParticle::theInstance = nullptr;


G4ChargedUnknownParticle* G4ChargedUnknownParticle::Definition() {
  if ( theInstance != nullptr ) return theInstance;
  const G4String name = "chargedunknown";
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle( name );
  if ( anInstance == nullptr ) {
    //    Arguments for constructor are as follows
    //               name             mass          width         charge
    //             2*spin           parity  C-conjugation
    //          2*Isospin       2*Isospin3       G-parity
    //               type    lepton number  baryon number   PDG encoding
    //             stable         lifetime    decay table
    //             shortlived      subType    anti_encoding
   anInstance = new G4ParticleDefinition(
	         name,         0.0*MeV,       0.0*MeV,  +1.0*eplus, 
		    0,               0,             0,          
		    0,               0,             0,             
	   "geantino",               0,             0,           0,
		 true,            -1.0,          nullptr,
		false,      "geantino",            0
		);
  }
  theInstance = static_cast< G4ChargedUnknownParticle* >( anInstance );
  return theInstance;
}


G4ChargedUnknownParticle* G4ChargedUnknownParticle::ChargedUnknownParticleDefinition() {
  return Definition();
}


G4ChargedUnknownParticle* G4ChargedUnknownParticle::ChargedUnknownParticle() {
  return Definition();
}
