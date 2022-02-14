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
// G4VRestProcess
// 
// Class description:
//
// Abstract class which defines the public behavior of
// physics interactions at rest.

// Authors:
// - 2 December 1995, G.Cosmo - First implementation, based on object model
// - 18 December 1996, H.Kurashige - New Physics scheme
// --------------------------------------------------------------------
#ifndef G4VRestProcess_hh
#define G4VRestProcess_hh 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VProcess.hh"

class G4VRestProcess : public G4VProcess 
{
  public:

    G4VRestProcess(const G4String& aName,
                   G4ProcessType aType = fNotDefined);
    G4VRestProcess(G4VRestProcess&);

    virtual ~G4VRestProcess();

    G4VRestProcess& operator=(const G4VRestProcess&) = delete;

    virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition
                            );

    virtual G4VParticleChange* AtRestDoIt(
                             const G4Track& ,
                             const G4Step&
                            );

    virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
                             G4double  ,
                             G4double  ,
                             G4double& ,
                             G4GPILSelection*
                           ) { return -1.0; }
      // No operation in AlongStepGPIL

    virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& ,
                             G4double   ,
                             G4ForceCondition* 
                            ) { return -1.0; }
      // No operation in PostStepGPIL

    virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& ,
                             const G4Step&
                            ) { return 0; }
      // No operation in PostStepDoIt

    virtual G4VParticleChange* AlongStepDoIt(
                             const G4Track& ,
                             const G4Step& 
                            ) { return 0; }
      // No operation in AlongStepDoIt

  protected:

    virtual G4double GetMeanLifeTime(const G4Track& aTrack,
                                     G4ForceCondition* condition) = 0;
      // Calculates the mean life-time (i.e. for decays) of the
      // particle at rest due to the occurrence of the given process,
      // or converts the probability of interaction (i.e. for
      // annihilation) into the life-time of the particle for the
      // occurrence of the given process

  private:

    G4VRestProcess();
      // Hidden default constructor
};

#endif
