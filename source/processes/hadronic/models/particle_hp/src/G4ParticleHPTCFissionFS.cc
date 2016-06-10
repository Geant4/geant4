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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPTCFissionFS.hh"

  void G4ParticleHPTCFissionFS::Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String &, G4ParticleDefinition* projectile )
  {
    G4String aString = "/TC/";
    G4ParticleHPFissionBaseFS::Init(A, Z, M, dirName, aString, projectile);
  }
  
  G4DynamicParticleVector * G4ParticleHPTCFissionFS::ApplyYourself(G4int NNeutrons)
  {  
    G4DynamicParticleVector * aResult;
//    G4cout <<"G4ParticleHPTCFissionFS::ApplyYourself +"<<G4endl;
    aResult = G4ParticleHPFissionBaseFS::ApplyYourself(NNeutrons);    
    return aResult;
  }
