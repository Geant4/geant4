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
// $Id: G4NeutronKiller.hh 68048 2013-03-13 14:34:07Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutronKiller
//
// Description: The process to kill neutrons to save CPU
//
// Author:      V.Ivanchenko 26/09/00 for HARP software
//
//----------------------------------------------------------------------------
//
// Class description:
//
// G4NeutronKiller allows to remove unwanted neutrons from simulation in 
// order to improve CPU performance. There are two parameters:
//                 
// 1) low energy threshold for neutron kinetic energy (default 0)
// 2) time limit for neutron track (default DBL_MAX) 
// 
// These parameters can be changed via Set methods or by UI commands:
//   /physics_engine/neutron/energyCut
//   /physics_engine/neutron/timeLimit
//
// If a neutron track is killed no energy deposition is added to the step 
//
// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef DNADEV
#ifndef NeutronKiller_h
#define NeutronKiller_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4Track.hh"

class G4TrackKillerMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
class G4TrackKillerManager
{

	G4TrackKillerManager();
	virtual ~G4TrackKillerManager();

	std::map<G4ParticleDefinition*, G4TrackKiller*>
public :
};
*/

class G4TrackKiller : public G4VDiscreteProcess, public G4UImessenger
{
public:

  G4TrackKiller(const G4String& processName = "nKiller",
		  G4ProcessType   aType =  fGeneral );

  virtual ~G4TrackKiller();

  G4bool IsApplicable(const G4ParticleDefinition&);

  void  SetTimeLimit(G4double);
  void  SetMinKinEnergyLimit(G4double);
  void  SetMaxKinEnergyLimit(G4double);

  G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
						 G4double previousStepSize,
						 G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  G4double GetMeanFreePath(const G4Track&, G4double,G4ForceCondition*);

private:

  // hide assignment operator as private
  G4TrackKiller(const G4TrackKiller&);
  G4TrackKiller& operator = (const G4TrackKiller &right);

  G4double fMinKinEnergyThreshold;
  G4double fMaxKinEnergyThreshold;
  G4double fTimeThreshold;
     
  G4TrackKillerMessenger* pMess;
  G4UIcmdWithAString* fp;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4TrackKiller::PostStepGetPhysicalInteractionLength(
				 const G4Track& aTrack,
				 G4double, G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  
  G4double limit = DBL_MAX; 
  if(aTrack.GetGlobalTime() > fTimeThreshold ||
     aTrack.GetKineticEnergy() < fMinKinEnergyThreshold) limit = 0.0;
  return limit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4TrackKiller::GetMeanFreePath(const G4Track&,G4double,
						 G4ForceCondition*)
{
  return DBL_MAX;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4VParticleChange* G4TrackKiller::PostStepDoIt(const G4Track& aTrack,
							const G4Step&)
{
  pParticleChange->Initialize(aTrack);
  pParticleChange->ProposeTrackStatus(fStopAndKill);
  return pParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
#endif

