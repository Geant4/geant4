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
//
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation
//      HPW, 10DEC 98, the decay part originally written by Gunter Folger
//                in his FTF-test-program.
//
//      V. Uzhinsky Nov. 2012 
//      introduced new method PropagateNuclNucl for nucleus-nucleus interactions
//
// -----------------------------------------------------------------------------
//
// Class Description
// Trivial implementation of an intra-nuclear transport. It pworvides coupling
// of high energy generators with pre equilibrium decay models.
// To be used in your physics list in case you need this physics.
// Class Description - End

#ifndef G4GeneratorPrecompoundInterface_h
#define G4GeneratorPrecompoundInterface_h 1

#include "G4VIntraNuclearTransportModel.hh"
#include "G4ReactionProductVector.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "globals.hh"

class G4KineticTrackVector;
class G4V3DNucleus;
class G4ParticleDefinition;

class G4GeneratorPrecompoundInterface : public G4VIntraNuclearTransportModel
{
public:

  G4GeneratorPrecompoundInterface(G4VPreCompoundModel* p = 0);
  virtual ~G4GeneratorPrecompoundInterface();

  virtual G4HadFinalState*
  ApplyYourself(const G4HadProjectile &aTrack, G4Nucleus &targetNucleus );

  virtual G4ReactionProductVector*
  Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus);

  virtual G4ReactionProductVector*
  PropagateNuclNucl(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus,
                                                  G4V3DNucleus* theProjectileNucleus);

  inline void SetCaptureThreshold(G4double);

  inline void SetDeltaM(G4double);
  inline void SetDeltaR(G4double);

  void MakeCoalescence(G4KineticTrackVector* theSecondaries);

  virtual void PropagateModelDescription(std::ostream&) const;

private:

  G4GeneratorPrecompoundInterface(const G4GeneratorPrecompoundInterface& right);
  const G4GeneratorPrecompoundInterface& operator=(const G4GeneratorPrecompoundInterface &right);
  G4bool operator==(G4GeneratorPrecompoundInterface& right) {return (this == &right);}
  G4bool operator!=(G4GeneratorPrecompoundInterface& right) {return (this != &right);}

  G4double CaptureThreshold;
  G4double DeltaM;
  G4double DeltaR;

  const G4ParticleDefinition* proton;
  const G4ParticleDefinition* neutron;
  const G4ParticleDefinition* lambda;

  const G4ParticleDefinition* deuteron;
  const G4ParticleDefinition* triton;
  const G4ParticleDefinition* He3;
  const G4ParticleDefinition* He4;

  const G4ParticleDefinition* ANTIproton;
  const G4ParticleDefinition* ANTIneutron;

  const G4ParticleDefinition* ANTIdeuteron;
  const G4ParticleDefinition* ANTItriton;
  const G4ParticleDefinition* ANTIHe3;
  const G4ParticleDefinition* ANTIHe4;

  G4int secID;  // Creator model ID
};

inline
void G4GeneratorPrecompoundInterface::SetCaptureThreshold(G4double value)
{
  CaptureThreshold = value;
}

inline
void G4GeneratorPrecompoundInterface::SetDeltaM(G4double value)
{ 
  DeltaM = value;
}

inline 
void G4GeneratorPrecompoundInterface::SetDeltaR(G4double value)
{ 
  DeltaR = value;
}

#endif // G4GeneratorPrecompoundInterface_h

