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
// This software was developed by Lawrence Livermore National Laboratory.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// 3. The name of the author may not be used to endorse or promote products
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Copyright (c) 2006 The Regents of the University of California.
// All rights reserved.
// UCRL-CODE-224807
//
// $Id: G4FissionLibrary.hh 67966 2013-03-13 09:38:38Z gcosmo $
//

#ifndef G4FissionLibrary_h
#define G4FissionLibrary_h 1

#include "globals.hh"
#include "G4HadProjectile.hh"
#include "G4DynamicParticleVector.hh"
#include "G4HadFinalState.hh"
#include "G4ParticleHPFinalState.hh"
#include "G4ParticleHPNames.hh"
#include "G4ParticleHPParticleYield.hh"
#include "G4ParticleHPFissionERelease.hh"
#include "G4ParticleHPEnergyDistribution.hh"
#include "G4ParticleHPPhotonDist.hh"
#include "G4ParticleHPAngular.hh"
#include "G4Nucleus.hh"
#include "Randomize.hh"
#include "G4LLNLFission.hh"
#include "G4fissionEvent.hh"

class G4FissionLibrary : public G4ParticleHPFinalState
{
  public:
  
  G4FissionLibrary();
  ~G4FissionLibrary();

  //void Init (G4double A, G4double Z, G4String & dirName, G4String &);
  //void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String &);
  void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String &, G4ParticleDefinition* );
  G4HadFinalState * ApplyYourself(const G4HadProjectile & theTrack);
  G4ParticleHPFinalState * New() ;

  private:
  G4fissionEvent* fe;
  G4int theIsotope; // used to call G4fissionEvent
  G4double targetMass;
  void SampleMult(const G4HadProjectile & theTrack, G4int* nPrompt,
                                   G4int* gPrompt, G4double eKinetic);
  inline G4ParticleHPFissionERelease * GetEnergyRelease() {
    return &theEnergyRelease;
  }
  G4ParticleHPParticleYield theFinalStateNeutrons;
  G4ParticleHPEnergyDistribution thePromptNeutronEnDis;
  G4ParticleHPEnergyDistribution theDelayedNeutronEnDis;
  G4ParticleHPAngular theNeutronAngularDis;

  G4ParticleHPPhotonDist theFinalStatePhotons;
  G4ParticleHPFissionERelease theEnergyRelease;
};
#endif
