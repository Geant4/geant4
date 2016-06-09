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
// $Id: G4ParametrizedHadronicVertex.cc,v 1.6 2006/06/29 20:57:40 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
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
    theTotalResult.ProposeLocalEnergyDeposit(0.);
    theTotalResult.Initialize(thePhoton);
    theTotalResult.ProposeTrackStatus(fAlive);
    G4double theKineticEnergy = thePhoton.GetKineticEnergy();
    G4HadFinalState * aR = 0;
    G4HadProjectile thePro(thePhoton);
    if(CLHEP::RandBit::shootBit())
    {
      if(theKineticEnergy<20*GeV) aR = theLowEPionMinus->ApplyYourself(thePro, theTarget);
      else aR = theHighEPionMinus->ApplyYourself(thePro, theTarget);
    }
    else
    {
      if(theKineticEnergy<20*GeV) aR = theLowEPionPlus->ApplyYourself(thePro, theTarget);
      else aR = theHighEPionPlus->ApplyYourself(thePro, theTarget);
    }
    aR->SetTrafoToLab(thePro.GetTrafoToLab());
    if(aR->GetStatusChange()==stopAndKill)
    {
      theTotalResult.ProposeTrackStatus(fStopAndKill);
      theTotalResult.ProposeEnergy( 0.0 );
    }
    if(aR->GetStatusChange()==suspend)
    {
      theTotalResult.ProposeTrackStatus(fSuspend);
    }
    if(aR->GetStatusChange()!=stopAndKill )
    {
      G4double newWeight = aR->GetWeightChange()*thePhoton.GetWeight();
      theTotalResult.ProposeParentWeight(newWeight); 
      if(aR->GetEnergyChange()>-.5) theTotalResult.ProposeEnergy(aR->GetEnergyChange());
      G4LorentzVector newDirection(aR->GetMomentumChange().unit(), 1.);
      newDirection*=aR->GetTrafoToLab();
      theTotalResult.ProposeMomentumDirection(newDirection.vect());
    }

    theTotalResult.ProposeLocalEnergyDeposit(aR->GetLocalEnergyDeposit());
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
