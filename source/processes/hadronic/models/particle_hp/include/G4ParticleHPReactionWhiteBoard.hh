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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPReactionWhiteBoard_h
#define G4ParticleHPReactionWhiteBoard_h 1

// Class Description
// Keep Information of Current Interaction
// Class Description - End

// 121031 First implementation done by T. Koi (SLAC/PPA)
//
#include "globals.hh"
#include <map>

class G4ParticleHPReactionWhiteBoard 
{
   public:
      G4ParticleHPReactionWhiteBoard();
      ~G4ParticleHPReactionWhiteBoard();

   public:
      void Dump();

      void SetTargZ( G4int Z ){ targZ=Z; };
      void SetTargA( G4int A ){ targA=A; };
      void SetTargM( G4int M ){ targM=M; };
      G4int GetTargZ(){ return targZ; };
      G4int GetTargA(){ return targA; };
      G4int GetTargM(){ return targM; };
       
      bool AddRecord( std::pair<G4String,G4String> );
      G4String GetValue( G4String );
      // "0" or "0.0" will retrun for invalid key by following two methods
      G4int GetValueInInt( G4String );
      G4double GetValueInDouble( G4String );

   private:
      G4int targZ;
      G4int targA;
      G4int targM;
      std::map< G4String,G4String > mapStringPair;
      //,,,
};

#endif
