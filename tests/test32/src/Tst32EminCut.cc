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
// $Id: Tst32EminCut.cc,v 1.1 2002-06-13 12:16:35 jwellisc Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
// ------------------------------------------------------------

#include "Tst32EminCut.hh"
#include "G4VParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"

Tst32EminCut::Tst32EminCut(const G4String& aName)
  : G4VProcess(aName),
    theEminThreshold( 1.0*MeV)
{
}

Tst32EminCut::~Tst32EminCut() 
{                                     
}                                     

G4VParticleChange* Tst32EminCut::PostStepDoIt(
			     const G4Track& aTrack,
			     const G4Step& 
			    )
//
// Stop the current particle 
// 			    			    			    
{
   aParticleChange.Initialize(aTrack);
   aParticleChange.SetLocalEnergyDeposit (aTrack.GetKineticEnergy()) ;
   aParticleChange.SetStatusChange(fStopAndKill);
   return &aParticleChange;
}

