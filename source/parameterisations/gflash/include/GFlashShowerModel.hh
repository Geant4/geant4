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
// $Id: GFlashShowerModel.hh,v 1.8 2005/11/30 19:29:44 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  GFlashShowerModel
//
//  Class description:
//
//  GFlash parameterisation shower model.

// Authors: E.Barberio & Joanna Weng - 9.11.04
//---------------------------------------------------------------
#ifndef GFlashShowerModel_h
#define GFlashShowerModel_h 1

//G4 Standard
#include "G4VFastSimulationModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ios.hh"

//GFlash
#include "GFlashShowerModelMessenger.hh"
#include "GFlashParticleBounds.hh"
#include "GFlashEnergySpot.hh"
#include "GFlashHitMaker.hh"
#include  <vector>

class GVFlashShowerParameterisation;
class GFlashHomoShowerParameterisation;
class GFlashSamplingShowerParameterisation;

class GFlashShowerModel : public G4VFastSimulationModel
{
  public:  // with description

    GFlashShowerModel (G4String, G4Envelope*);
    GFlashShowerModel (G4String);
    ~GFlashShowerModel ();  
      // Constructors, destructor

    G4bool ModelTrigger(const G4FastTrack &); 
    G4bool IsApplicable(const G4ParticleDefinition&);
    void DoIt(const G4FastTrack&, G4FastStep&);
      // Checks whether conditions of fast parameterisation are fullfilled
  
    // setting

    inline void SetFlagParamType(G4int I)
      { FlagParamType = I; }
    inline void SetFlagParticleContainment(G4int I)
      { FlagParticleContainment = I; }
    inline void SetStepInX0(G4double Lenght)
      { StepInX0=Lenght; } 
    inline void SetParameterisation(GVFlashShowerParameterisation &DP)
      { Parameterisation=&DP;}
    inline void SetHitMaker(GFlashHitMaker &Maker)
      { HMaker=&Maker; }
    inline void SetParticleBounds(GFlashParticleBounds &SpecificBound)
      { PBound =&SpecificBound; }
  
    // getting

    inline G4int GetFlagParamType()
      { return FlagParamType; }
    inline G4int GetFlagParticleContainment()
      { return FlagParticleContainment; }  
    inline G4double GetStepInX0()
      { return StepInX0; }

  public:  // without description

    // Gets ?  
    GFlashParticleBounds  *PBound;
    GVFlashShowerParameterisation *Parameterisation;  

  private:

    void ElectronDoIt(const G4FastTrack&, G4FastStep&);
    //  void GammaDoIt(const G4FastTrack&, G4FastStep&);
    //  void NeutrinoDoIt(const G4FastTrack&, G4FastStep&);
    G4bool CheckParticleDefAndContainment(const G4FastTrack &fastTrack);
    G4bool CheckContainment(const G4FastTrack &fastTrack);
  
  private:

    GFlashHitMaker *HMaker;  
    GFlashShowerModelMessenger* Messenger;
  
    //Control Flags
    G4int FlagParamType;           ///0=no GFlash 1=only em showers parametrized
    G4int FlagParticleContainment; ///0=no check  ///1=only fully contained...
    G4double StepInX0;  
    G4double EnergyStop;
  
};
#endif
