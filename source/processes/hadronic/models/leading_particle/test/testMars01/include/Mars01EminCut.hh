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
// $Id: Mars01EminCut.hh,v 1.1 2001-12-13 14:58:42 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
// ------------------------------------------------------------

#ifndef Mars01EminCut_h
#define Mars01EminCut_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VProcess.hh"


class Mars01EminCut : public G4VProcess 
{
  public:     

     Mars01EminCut(const G4String& processName ="ExN05SpecialCut" );

     virtual ~Mars01EminCut();

     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );

     virtual G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );
			    
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

     void SetEminThreshold(G4double );
     G4double GetEminThreshold() const; 
  private:
  
  // hide assignment operator as private 
     Mars01EminCut& operator=(const Mars01EminCut& right){return *this;};

  G4double theEminThreshold;

};

inline
 void Mars01EminCut::SetEminThreshold(G4double value)
{
  theEminThreshold = value;
}

inline
 G4double Mars01EminCut::GetEminThreshold() const
{
  return  theEminThreshold ;
}

inline 
  G4double Mars01EminCut::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            )
{
   // condition is set to "Not Forced"
  *condition = NotForced;

  G4double     proposedStep = DBL_MAX;
  G4double     eKine = track.GetDynamicParticle()->GetKineticEnergy();
  if (eKine <  theEminThreshold) proposedStep = DBL_MIN;
  return proposedStep;
}

#endif


