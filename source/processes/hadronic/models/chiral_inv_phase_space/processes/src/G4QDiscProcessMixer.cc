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
// $Id$
//
//      ---------------- G4QDiscProcessMixer class -------------------
//                 by Mikhail Kossov, Aug 2007.
// G4QDiscProcessMixer class of the CHIPS Simulation Branch in GEANT4
// ------------------------------------------------------------------------
// Short description: universal mixer of processes (NOT models as in GHAD!)
// depending on the application energy region (defined by users).
// ------------------------------------------------------------------------

//#define debug

#include "G4QDiscProcessMixer.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicDeprecate.hh"


// Constructor
G4QDiscProcessMixer::G4QDiscProcessMixer(const G4String& name,
                                         const G4ParticleDefinition* particle,
                                         G4ProcessType pType):
  G4VDiscreteProcess(name, pType), DPParticle(particle)
{
  G4HadronicDeprecate("G4QDiscProcessMixer");

#ifdef debug
  G4cout<<"G4QDiscProcessMixer::Constructor is called processName="<<name<<G4endl;
#endif
  if (verboseLevel>0) G4cout<<GetProcessName()<<" process is created "<<G4endl;
}

// Destructor
G4QDiscProcessMixer::~G4QDiscProcessMixer()
{
  // Now the responsibility of deleting is deligated to the user, who created them
  //for_each(theDPVector.begin(), theDPVector.end(), DeleteDiscreteProcess());
}

void G4QDiscProcessMixer::AddDiscreteProcess(G4VDiscreteProcess* DP, G4double ME)
{
  static const G4double maxEn = 1.E8*megaelectronvolt; // Conditional INF
  if(!theDPVector.size()) // The first process in the DiscreteProcessVector (MaxE=INF)
  {
    std::pair<G4VDiscreteProcess*, G4double>* QDiscProc =
                                   new std::pair<G4VDiscreteProcess*, G4double>(DP, maxEn);
    theDPVector.push_back( QDiscProc );
  }
  else
  {
    if(ME < theDPVector[theDPVector.size()-1]->second)
    {
      std::pair<G4VDiscreteProcess*, G4double>* QDiscProc =
                                      new std::pair<G4VDiscreteProcess*, G4double>(DP, ME);
      theDPVector.push_back( QDiscProc );
    }
    else // Wrong Max Energy Order for the new process in the sequence of processes
    {
      // G4cerr<<"G4QDiscProcessMixer::AddDiscreteProcess:LastMaxE("<<theDPVector.size()-1
      //       <<")="<<theDPVector[theDPVector.size()-1]->second<<" <= MaxE="<<ME<<G4endl;
      // G4Exception("G4QDiscProcessMixer::AddDiscreteProcess: Wrong Max Energy Order");
      G4ExceptionDescription ed;
      ed << " LastMaxE(" << theDPVector.size()-1 << ")="
         << theDPVector[theDPVector.size()-1]->second << " <= MaxE=" << ME << G4endl;
      G4Exception("G4QDiscProcessMixer::AddDiscreteProcess()", "HAD_CHPS_0000",
                  FatalException, ed);
    }
  }
}

G4bool G4QDiscProcessMixer::IsApplicable(const G4ParticleDefinition& particle) 
{
#ifdef debug
  G4cout<<"G4QDiscProcessMixer::IsApplicable: part="<<particle.GetParticleName()<<" = "
        <<DPParticle->GetParticleName()<<G4endl;
#endif
  if(particle == *DPParticle) return true;
  return false;
}

G4double G4QDiscProcessMixer::PostStepGetPhysicalInteractionLength(const G4Track& Track,
                                                                   G4double  PrevStSize,
                                                                   G4ForceCondition* F)
{
  G4double kEn=Track.GetDynamicParticle()->GetKineticEnergy(); // Projectile kinetic energy
  G4int maxDP=theDPVector.size();
  if(maxDP) for(G4int i=maxDP-1; i>-1; i--) if(kEn < theDPVector[i]->second)
  {
#ifdef debug
    G4cout<<"G4QDPMix::PSGetPIL:"<<i<<",kEn="<<kEn<<" < "<<theDPVector[i]->second<<", PIL="
        <<theDPVector[i]->first->PostStepGetPhysicalInteractionLength(Track,PrevStSize,F)
        <<G4endl;
#endif
    return theDPVector[i]->first->PostStepGetPhysicalInteractionLength(Track,PrevStSize,F);
  }
  return DBL_MAX;
}

// compilation fake class (length is calculated in PostStepGetPhysicalInteractionLength)
G4double G4QDiscProcessMixer::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*)
{
  return DBL_MAX;
}

G4VParticleChange* G4QDiscProcessMixer::PostStepDoIt(const G4Track& Track,
                                                     const G4Step& Step)
{
  G4double kEn=Track.GetDynamicParticle()->GetKineticEnergy(); // Projectile kinetic energy
  G4int maxDP=theDPVector.size();
  if(maxDP) for(G4int i=maxDP-1; i>-1; i--) if(kEn < theDPVector[i]->second)
  {
#ifdef debug
    G4cout<<"G4QDPMix::PSDoIt: i="<<i<<",kEn="<<kEn<<" < "<<theDPVector[i]->second<<G4endl;
#endif
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
