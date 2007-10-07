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
// $Id: G4DNAProcess.hh,v 1.1 2007-10-07 12:48:36 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author:  Maria Grazia Pia (Maria.Grazia.Pia@ge.infn.it)
//
// History:
// -----------
// 29 Sep 2007 MGP   First (incomplete) impleentation
//
// -------------------------------------------------------------------

// Class description:
// Host class for DNA physics policies
// Further documentation available in PAPER REF. TO BE ADDED

// -------------------------------------------------------------------

#ifndef G4DNAPROCESS_HH
#define G4DNAPROCESS_HH 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;

template < class TCrossSection, class TFinalState>
class G4DNAProcess : public G4VDiscreteProcess {
  
public:

  // ---- MGP ---- Note: process name to be replaced with a better identifying mechanism  
  G4DNAProcess(const G4String& processName = "DNAProcess")  { /* nop */; }
  
  ~G4DNAProcess() { /* nop */; }

  //  ---- MGP ---- Dummy initially: process is always applicable
  G4bool IsApplicable(const G4ParticleDefinition&) { return true; } 
  
  //  ---- MGP ---- Dummy initially: no PhysicsTable (usefulness to be verified)
  void BuildPhysicsTable(const G4ParticleDefinition& particle) { /* nop */; }
 
 // For testing purpose only
  G4double DumpMeanFreePath(const G4Track& aTrack, 
			    G4double previousStepSize, 
			    G4ForceCondition* condition) 
  { return GetMeanFreePath(aTrack, previousStepSize, condition); }

  // Implemented in terms of TFinalState
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step)
  {
    aParticleChange.Initialize(track);

    // Interaction product
    const G4FinalStateProduct& product = finalState.GenerateFinalState(track,step);

    // Number of secondary products to be generated
    G4int nSecondaries  = product.NumberOfSecondaries();
    aParticleChange.SetNumberOfSecondaries(nSecondaries);

    // Secondaries
    for (G4int l = 0; l<nSecondaries; l++ )
    {
      G4VDynamicParticle* particle = product.GetSecondaries()[l];
      if (particle != 0) 
	{
	  aParticleChange.AddSecondary(particle);
	}
    }

   // Take care of incident particle to be killed, if necessary; dump its energy deposit locally
    if (product.KillIncidentParticle())
      {
	G4double deposit = product.EnergyDeposit();
	aParticleChange.ProposeLocalEnergyDeposit(deposit);
	aParticleChange.ProposeTrackStatus(fStopAndKill);
	aParticleChange.ProposeEnergy(0.);
	aParticleChange.ProposeMomentumDirection( 0., 0., 0. );
      }

    return G4VDiscreteProcess::PostStepDoIt(track,step );
  }

protected:
  
  // Implemented in terms of TCrossSection
  G4double GetMeanFreePath(const G4Track& track, 
			   G4double previousStepSize, 
			   G4ForceCondition* condition)
  { 
    G4double meanFreePath = DBL_MAX;
    G4double cross = crossSection.CalculateCrossSection(track);
    if (cross > 0.0) meanFreePath = 1.0/cross;
    return meanFreePath;
  }

private:

 // Hide copy constructor and assignment operator as private 
  G4DNAProcess& operator=(const G4DNAProcess& right);
  G4DNAProcess(const G4DNAProcess& );

  // Policy classes
  TCrossSection crossSection;
  TFinalState finalState;

};

#endif



