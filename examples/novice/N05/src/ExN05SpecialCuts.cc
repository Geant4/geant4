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
// $Id: ExN05SpecialCuts.cc,v 1.9 2006/06/29 17:53:44 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class implementation file 
//
// ------------------------------------------------------------
//                  15 April 1998  M.Maire
// ------------------------------------------------------------

#include "ExN05SpecialCuts.hh"
#include "G4VParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"

ExN05SpecialCuts::ExN05SpecialCuts(const G4String& aName)
  : G4VProcess(aName)
{
   if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

ExN05SpecialCuts::~ExN05SpecialCuts() 
{                                     
}                                     

G4VParticleChange*
ExN05SpecialCuts::PostStepDoIt( const G4Track& aTrack, const G4Step& )
{
  //
  // Stop the current particle, if requested by G4UserLimits 
  // 			    			    			    
  aParticleChange.Initialize(aTrack);
  aParticleChange.ProposeEnergy(0.) ;
  aParticleChange.ProposeLocalEnergyDeposit (aTrack.GetKineticEnergy()) ;
  aParticleChange.ProposeTrackStatus(fStopButAlive);
  return &aParticleChange;
}

G4double
ExN05SpecialCuts::
PostStepGetPhysicalInteractionLength(const G4Track&,G4double,G4ForceCondition*)
{
  return DBL_MAX;
}
