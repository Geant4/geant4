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
// $Id: G4ParametrizedHadronicVertex.cc,v 1.1 2003/11/11 19:08:58 hpw Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// --------------------------------------------------------------
#include "G4ParametrizedHadronicVertex.hh"
#include "G4HadFinalState.hh"
#include "G4ParticleChange.hh"

G4VParticleChange * G4ParametrizedHadronicVertex::
ApplyYourself(G4Nucleus & theTarget, const G4Track &thePhoton)
{   
    static G4ParticleChange theTotalResult; 

    theTotalResult.Clear();
    theTotalResult.SetLocalEnergyDeposit(0.);
    theTotalResult.Initialize(thePhoton);
    theTotalResult.SetStatusChange(fAlive);
    G4double theKineticEnergy = thePhoton.GetKineticEnergy();
    G4HadFinalState * aR = 0;
    G4HadProjectile thePro(thePhoton);
    if(RandBit::shootBit())
    {
      if(theKineticEnergy<20*GeV) aR = theLowEPionMinus.ApplyYourself(thePro, theTarget);
      else aR = theHighEPionMinus.ApplyYourself(thePro, theTarget);
    }
    else
    {
      if(theKineticEnergy<20*GeV) aR = theLowEPionPlus.ApplyYourself(thePro, theTarget);
      else aR = theHighEPionPlus.ApplyYourself(thePro, theTarget);
    }
    aR->SetTrafoToLab(thePro.GetTrafoToLab());
    if(aR->GetStatusChange()==stopAndKill)
    {
      theTotalResult.SetStatusChange(fStopAndKill);
      theTotalResult.SetEnergyChange( 0.0 );
    }
    if(aR->GetStatusChange()==suspend)
    {
      theTotalResult.SetStatusChange(fSuspend);
    }
    if(aR->GetStatusChange()!=stopAndKill )
    {
      G4double newWeight = aR->GetWeightChange()*thePhoton.GetWeight();
      theTotalResult.SetWeightChange(newWeight); 
      if(aR->GetEnergyChange()>-.5) theTotalResult.SetEnergyChange(aR->GetEnergyChange());
      G4LorentzVector newDirection(aR->GetMomentumChange().unit(), 1.);
      newDirection*=aR->GetTrafoToLab();
      theTotalResult.SetMomentumDirectionChange(newDirection.vect());
    }

    theTotalResult.SetLocalEnergyDeposit(aR->GetLocalEnergyDeposit());
    theTotalResult.SetNumberOfSecondaries(aR->GetNumberOfSecondaries());

    G4ThreeVector it(0., 0., 1.);
    G4double what = 2.*pi*G4UniformRand();
    for(G4int i=0; i<aR->GetNumberOfSecondaries(); i++)
    {
      G4LorentzVector theM = aR->GetSecondary(i)->GetParticle()->Get4Momentum();
      theM.rotate(what, it);
      theM*=aR->GetTrafoToLab();
      aR->GetSecondary(i)->GetParticle()->Set4Momentum(theM);
      G4double time = aR->GetSecondary(i)->GetTime();
      if(time<0) time = thePhoton.GetGlobalTime();
      G4Track* track = new G4Track(aR->GetSecondary(i)->GetParticle(),
				   thePhoton.GetGlobalTime(),
				   thePhoton.GetPosition());
      G4double newWeight = thePhoton.GetWeight()*aR->GetSecondary(i)->GetWeight();
      track->SetWeight(newWeight);
      theTotalResult.AddSecondary(track);
    }
  return &theTotalResult;
}
