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
// $Id: GFlashParticleBounds.hh,v 1.3 2005/10/04 09:08:33 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

