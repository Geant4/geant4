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
// $Id: G4QDiscProcessMixer.cc,v 1.2 2007/08/31 09:36:57 mkossov Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//      ---------------- G4QDiscProcessMixer class -----------------
//                 by Mikhail Kossov, Aug 2007.
// G4QDiscProcessMixer class of the CHIPS Simulation Branch in GEANT4
// ---------------------------------------------------------------
// ****************************************************************************************
// ********** This CLASS is temporary moved from the "chips/interface" directory *********
// ****************************************************************************************

//#define debug

#include "G4QDiscProcessMixer.hh"

// Constructor
G4QDiscProcessMixer::G4QDiscProcessMixer(const G4String& name,
                                         const G4ParticleDefinition* particle,
                                         G4ProcessType pType):
  G4VDiscreteProcess(name, pType), DPParticle(particle)
{
#ifdef debug
  G4cout<<"G4QDiscProcessMixer::Constructor is called processName="<<name<<G4endl;
#endif
  if (verboseLevel>0) G4cout<<GetProcessName()<<" process is created "<<G4endl;
}

// Destructor
G4QDiscProcessMixer::~G4QDiscProcessMixer()
{
  for_each(theDPVector.begin(), theDPVector.end(), DeleteDiscreteProcess());
}

void G4QDiscProcessMixer::AddDiscreteProcess(G4VDiscreteProcess* DP, G4double ME)
{
  if(ME>theDPVector[theDPVector.size()-1]->second)
  {
    std::pair<G4VDiscreteProcess*, G4double>* QDiscProc =
      new std::pair<G4VDiscreteProcess*, G4double>(DP, ME);
    theDPVector.push_back( QDiscProc );
  }
  else // Wrong Max Energy Order for the new process in the sequence of processes
  {
    G4cerr<<"G4QDiscProcessMixer::AddDiscreteProcess:LastMaxE("<<theDPVector.size()-1<<")="
          <<theDPVector[theDPVector.size()-1]->second<<" >= MaxE="<<ME<<G4endl;
    G4Exception("G4QDiscProcessMixer::AddDiscreteProcess: Wrong Max Energy Order");
  }
}

G4bool G4QDiscProcessMixer::IsApplicable(const G4ParticleDefinition& particle) 
{
  if(particle == *DPParticle) return true;
  return false;
}

G4double G4QDiscProcessMixer::PostStepGetPhysicalInteractionLength(const G4Track& Track,
			                                                         G4double   PrevStSize,
			                                                         G4ForceCondition* F)
{
  G4double kEn=Track.GetDynamicParticle()->GetKineticEnergy(); // Projectile kinetic energy
  G4int maxDP=theDPVector.size();
  if(maxDP) for(G4int i=0; i<maxDP; ++i) if(kEn < theDPVector[i]->second)
    return theDPVector[i]->first->PostStepGetPhysicalInteractionLength(Track,PrevStSize,F);
  return DBL_MAX;
}

G4VParticleChange* G4QDiscProcessMixer::PostStepDoIt(const G4Track& Track,
                                                     const G4Step& Step)
{
  G4double kEn=Track.GetDynamicParticle()->GetKineticEnergy(); // Projectile kinetic energy
  G4int maxDP=theDPVector.size();
  if(maxDP) for(G4int i=0; i<maxDP; ++i) if(kEn < theDPVector[i]->second)
  {
    //EnMomConservation= theDPVector[i]->first->GetEnegryMomentumConservation();
    //nOfNeutrons      = theDPVector[i]->first->GetNumberOfNeutronsInTarget();
    return theDPVector[i]->first->PostStepDoIt(Track, Step);
  }
  return G4VDiscreteProcess::PostStepDoIt(Track, Step);
}

//G4LorentzVector G4QDiscProcessMixer::GetEnegryMomentumConservation()
//                                                              {return EnMomConservation;}

//G4int G4QDiscProcessMixer::GetNumberOfNeutronsInTarget()
//                                                                    {return nOfNeutrons;}
