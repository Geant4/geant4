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
// $Id: G4EvtBiasMechanism.hh,v 1.4 2001-07-11 10:08:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// 
// ------------------------------------------------------------
//  Class Description
//  This is a simple example for Event Biasing
//    This mechanism is applied to only particleToBeBiased.
//    If a track currently processed is the specifed particle type,  
//    ApplyMath methods multiply all secondaries in particle change 
//    by factor of MultiplicationForSecondaries with their weights of 
//    (parent weight)/(MultiplicationForSecondaries)

#ifndef G4EvtBiasMechanism_h
#define G4EvtBiasMechanism_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEvtBiasMechanism.hh"

class G4VParticleChange;

class G4EvtBiasMechanism :public G4VEvtBiasMechanism
{
public: //with description
  // constructors
  G4EvtBiasMechanism(const G4String& name = "EBSample", G4int mulFactor=10);
  G4EvtBiasMechanism(const G4EvtBiasMechanism&);

public: 
  virtual ~G4EvtBiasMechanism();

public: //with description
  // virtual methods derived from G4VEvtBiasMechanism
  virtual G4VParticleChange* ApplyMath( G4VParticleChange*, const G4Step& );
  virtual G4bool             IsApplicable(G4ParticleDefinition*) const;

  // Set/Get particle type to be biased
  void                         SetParticleBiased( G4ParticleDefinition* );
  const G4ParticleDefinition*  GetParticleBiased( ) const;


 private:
  G4ParticleDefinition* particleToBeBiased;
  const G4int MultiplicationForSecondaries;
};

inline 
 G4bool  G4EvtBiasMechanism::IsApplicable(G4ParticleDefinition* particle) const
{
  return (particle == particleToBeBiased);
}

inline 
 void  G4EvtBiasMechanism::SetParticleBiased(G4ParticleDefinition* particle)
{
  particleToBeBiased = particle;
}

inline 
 const G4ParticleDefinition*  G4EvtBiasMechanism::GetParticleBiased() const 
{
  return particleToBeBiased;
}

#endif

