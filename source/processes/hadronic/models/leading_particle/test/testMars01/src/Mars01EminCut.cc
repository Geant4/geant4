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
// $Id: Mars01EminCut.cc,v 1.1 2001-12-13 14:58:46 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
// ------------------------------------------------------------

#include "Mars01EminCut.hh"
#include "G4VParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"

Mars01EminCut::Mars01EminCut(const G4String& aName)
  : G4VProcess(aName),
    theEminThreshold( 1.0*MeV)
{
}

Mars01EminCut::~Mars01EminCut() 
{                                     
}                                     

G4VParticleChange* Mars01EminCut::PostStepDoIt(
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

