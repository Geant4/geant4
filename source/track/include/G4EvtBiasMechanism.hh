// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EvtBiasMechanism.hh,v 1.1 1999-01-07 16:14:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
// 
// ------------------------------------------------------------
//  This is a simple example for Event Biasing

#ifndef G4EvtBiasMechanism_h
#define G4EvtBiasMechanism_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEvtBiasMechanism.hh"

class G4VParticleChange;

class G4EvtBiasMechanism :public G4VEvtBiasMechanism
{
 public:
  G4EvtBiasMechanism(const G4String& name = "EBSample", G4int mulFactor=10);
  G4EvtBiasMechanism(const G4EvtBiasMechanism&);

  ~G4EvtBiasMechanism();

  virtual G4VParticleChange* ApplyMath( G4VParticleChange*, const G4Step& );

  void                   SetParticleBiased( G4ParticleDefinition* );
  G4ParticleDefinition*  GetParticleBiased( ) const;
  virtual G4bool IsApplicable(G4ParticleDefinition*) const;

 private:
  G4ParticleDefinition* particleToBeBiased;
  const G4int MultiplicationForSecondaries;
};
inline G4bool  G4EvtBiasMechanism::IsApplicable(G4ParticleDefinition* particle) const
{
  return (particle == particleToBeBiased);
}
inline void  G4EvtBiasMechanism::SetParticleBiased(G4ParticleDefinition* particle)
{
  particleToBeBiased = particle;
}

inline G4ParticleDefinition*  G4EvtBiasMechanism::GetParticleBiased() const 
{
  return particleToBeBiased;
}

#endif

