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
// $Id: Tst32EminCut.hh,v 1.3 2006-06-29 21:58:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
// ------------------------------------------------------------

#ifndef Tst32EminCut_h
#define Tst32EminCut_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VProcess.hh"


class Tst32EminCut : public G4VProcess 
{
  public:     

     Tst32EminCut(const G4String& processName ="ExN05SpecialCut" );

     virtual ~Tst32EminCut();

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
     Tst32EminCut& operator=(const Tst32EminCut&){return *this;};

  G4double theEminThreshold;

};

inline
 void Tst32EminCut::SetEminThreshold(G4double value)
{
  theEminThreshold = value;
}

inline
 G4double Tst32EminCut::GetEminThreshold() const
{
  return  theEminThreshold ;
}

inline 
  G4double Tst32EminCut::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double, //   previousStepSize,
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


