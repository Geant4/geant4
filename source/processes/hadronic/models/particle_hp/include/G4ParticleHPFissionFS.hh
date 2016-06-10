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
#ifndef G4ParticleHPFissionFS_h
#define G4ParticleHPFissionFS_h 1

#include "globals.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"
#include "G4ParticleHPFinalState.hh"
#include "G4ParticleHPNames.hh"

#include "G4ParticleHPFCFissionFS.hh"
#include "G4ParticleHPSCFissionFS.hh"
#include "G4ParticleHPTCFissionFS.hh"
#include "G4ParticleHPLCFissionFS.hh"
#include "G4ParticleHPFSFissionFS.hh"

#include "G4ParticleHPFFFissionFS.hh"

class G4ParticleHPFissionFS : public G4ParticleHPFinalState
{
  public:
  
  G4ParticleHPFissionFS(){ hasXsec = false; produceFissionFragments = false; }
  ~G4ParticleHPFissionFS(){}
  void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & aFSType, G4ParticleDefinition* );
  G4HadFinalState * ApplyYourself(const G4HadProjectile & theTrack);
  G4ParticleHPFinalState * New() 
  {
   G4ParticleHPFissionFS * theNew = new G4ParticleHPFissionFS;
   return theNew;
  }
        
  private:
  
  G4ParticleHPFSFissionFS theFS;
  G4ParticleHPFCFissionFS theFC;
  G4ParticleHPSCFissionFS theSC;
  G4ParticleHPTCFissionFS theTC;
  G4ParticleHPLCFissionFS theLC;
    
  G4ParticleHPFFFissionFS theFF;
  G4bool produceFissionFragments;
};
#endif
