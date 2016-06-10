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
// $Id: GFlashParticleBounds.hh 68057 2013-03-13 14:46:00Z gcosmo $
//
//
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  GFlashParticleBounds
//
//  Class description:
//
//  GFlash parameterisation particle bounds.

//
// Author: Joanna Weng - 9.11.04
//---------------------------------------------------------------
#ifndef GFlashParticleBounds_h
#define GFlashParticleBounds_h 

#include  "G4ParticleDefinition.hh"

class  GFlashParticleBounds
{
  public:  // with description

    GFlashParticleBounds();
    ~GFlashParticleBounds();
  
    // methods to get/set ELE/Gamma max & min energy bounds

    G4double GetMinEneToParametrise(G4ParticleDefinition &particleType);
    G4double GetMaxEneToParametrise(G4ParticleDefinition &particleType); 
    G4double GetEneToKill(G4ParticleDefinition &particleType) ;
  
    void SetMinEneToParametrise(G4ParticleDefinition &particleType,
                                G4double enemin);
    void SetMaxEneToParametrise(G4ParticleDefinition &particleType,
                                G4double enemax);
    void SetEneToKill(G4ParticleDefinition &particleType,
                                G4double enekill);

  private:
  
    // electron and positron
    G4double EMinEneToParametrise;
    G4double EMaxEneToParametrise;
    G4double EEneToKill;
};
#endif

