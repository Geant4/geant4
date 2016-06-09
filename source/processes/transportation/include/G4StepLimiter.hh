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
// $Id: G4StepLimiter.hh,v 1.2 2005/03/07 16:21:41 asaim Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// class description
//
// A "process" to be registered to the process manager of each particle,
// in the UserPhysicsList, in order to take into account the MaxAllowedStep
// defined in the G4UserLimits object attached to a logical volume.
//
// ------------------------------------------------------------
//                  23 Jan. 2004  H.Kurashige
// ------------------------------------------------------------
#ifndef G4StepLimiter_h
#define G4StepLimiter_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VProcess.hh"

class G4StepLimiter : public G4VProcess 
{
  public:  // with description     

     G4StepLimiter(const G4String& processName ="StepLimiter" );

     virtual ~G4StepLimiter();

     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );

     virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& ,
                             const G4Step& 
                            );
                            
  public:  // without description 
                                 
     //  no operation in  AtRestGPIL
     virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
                             G4ForceCondition* 
                            ){ return -1.0; };
                            
     //  no operation in  AtRestDoIt      
     virtual G4VParticleChange* AtRestDoIt(
                             const G4Track& ,
                             const G4Step&
                            ){return 0;};

     //  no operation in  AlongStepGPIL
     virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
                             G4double  ,
                             G4double  ,
                             G4double& ,
                             G4GPILSelection*
                            ){ return -1.0; };

     //  no operation in  AlongStepDoIt
     virtual G4VParticleChange* AlongStepDoIt(
                             const G4Track& ,
                             const G4Step& 
                            ) {return 0;};

  private:
  
  // hide assignment operator as private 
      G4StepLimiter(G4StepLimiter&);
      G4StepLimiter& operator=(const G4StepLimiter& right);

};

#endif










