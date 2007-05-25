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
// $Id: G4ProcessWrappers.hh,v 1.1 2007-05-25 19:14:37 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation - wrapper definitions.
//
#ifndef G4PROCESSWRAPPERS_HH
#define G4PROCESSWRAPPERS_HH

#include "G4ForceCondition.hh"
#include "G4Functor.hh"
#include "G4GPILSelection.hh"
#include "G4ProcessLists.hh"
#include "G4TypeList.hh"

class G4FunctorIdentifier;
class G4Step;
class G4Track;
class G4VParticleChange;

namespace G4ProcessWrappers {
    
  // Rest, continuous and discrete "DoIt" methods have same signature
  // I.e : G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&)
  //       G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&)
  //       G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&)
  typedef G4Functor<G4VParticleChange*, G4FunctorIdentifier, G4TypeList_2(G4Track, G4Step)> DoItWrapper;
  
  // DoIt relay signature defined as: 
  //       G4VParticleChange* func(DoItWrapper&, const G4Track&, const G4Step&)
  typedef G4Functor<G4VParticleChange*, G4FunctorIdentifier, G4TypeList_3(DoItWrapper&, G4Track, G4Step)> DoItRelayWrapper;

  // Formal wrapper definitions for rest, continuous and discrete processes, for both "DoIt" and "GPIL" methods
  template <typename List> struct Wrappers {};
  
  template <> struct Wrappers<G4ProcessLists::RestDoIt> 
  {
    typedef DoItWrapper SeedWrapper;
    typedef DoItRelayWrapper RelayWrapper;
  };
  
  template <> struct Wrappers<G4ProcessLists::ContinuousDoIt> 
  {
    typedef DoItWrapper SeedWrapper;
    typedef DoItRelayWrapper RelayWrapper;
  };

  template <> struct Wrappers<G4ProcessLists::DiscreteDoIt> 
  {
    typedef DoItWrapper SeedWrapper;
    typedef DoItRelayWrapper RelayWrapper;
  };
  

  // Rest, continuous and discrete GPIL signatures are all different :
  //        G4double AtRestGPIL(const G4Track& track, G4ForceCondition* condition)
  //        G4double AlongStepGPIL(const G4Track& track,       
  //                               G4double previousStepSize,     
  //                               G4double currentMinimumStep,
  //                               G4double& proposedSafety,    
  //                               G4GPILSelection* selection);
  //        G4double PostStepGPIL(const G4Track& track,
  //                              G4double previousStepSize,
  //                              G4ForceCondition* condition);
  //
  // GPIL relays defined as :
  //       G4double func(SeedWrapper&, SeedArgs) 
  template <> struct Wrappers<G4ProcessLists::RestGPIL> 
  {
    typedef G4TypeList_2(G4Track, G4ForceCondition*) SeedArgs;
    typedef G4Functor<G4double, G4FunctorIdentifier, SeedArgs> SeedWrapper;

    typedef G4TypeList_3(SeedWrapper&, G4Track, G4ForceCondition*) RelayArgs;
    typedef G4Functor<G4double, G4FunctorIdentifier, RelayArgs> RelayWrapper;
  };

  template <> struct Wrappers<G4ProcessLists::ContinuousGPIL> 
  {
    typedef G4TypeList_5(G4Track, G4double, G4double, G4double&, G4GPILSelection*) SeedArgs;
    typedef G4Functor<G4double, G4FunctorIdentifier, SeedArgs> SeedWrapper;

    typedef G4TypeList_6(SeedWrapper&, G4Track, G4double, G4double, G4double&, G4GPILSelection*) RelayArgs;
    typedef G4Functor<G4double, G4FunctorIdentifier, RelayArgs> RelayWrapper;
  };

  template <> struct Wrappers<G4ProcessLists::DiscreteGPIL> 
  {
    typedef G4TypeList_3(G4Track, G4double, G4ForceCondition*) SeedArgs;
    typedef G4Functor<G4double, G4FunctorIdentifier, SeedArgs> SeedWrapper;

    typedef G4TypeList_4(SeedWrapper&, G4Track, G4double, G4ForceCondition*) RelayArgs;
    typedef G4Functor<G4double, G4FunctorIdentifier, RelayArgs> RelayWrapper;
  };
}

// Convenient typedefs
typedef G4ProcessWrappers::Wrappers<G4ProcessLists::RestGPIL>::SeedWrapper G4RestGPILWrapper;
typedef G4ProcessWrappers::Wrappers<G4ProcessLists::ContinuousGPIL>::SeedWrapper G4ContinuousGPILWrapper;
typedef G4ProcessWrappers::Wrappers<G4ProcessLists::DiscreteGPIL>::SeedWrapper G4DiscreteGPILWrapper;

typedef G4ProcessWrappers::Wrappers<G4ProcessLists::RestDoIt>::SeedWrapper G4RestDoItWrapper;
typedef G4ProcessWrappers::Wrappers<G4ProcessLists::ContinuousDoIt>::SeedWrapper G4ContinuousDoItWrapper;
typedef G4ProcessWrappers::Wrappers<G4ProcessLists::DiscreteDoIt>::SeedWrapper G4DiscreteDoItWrapper;

typedef G4ProcessWrappers::DoItWrapper G4DoItWrapper;

#endif
