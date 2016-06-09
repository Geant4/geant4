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
// $Id: G4ParametrizedHadronicVertex.cc,v 1.7 2009-03-04 19:09:20 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#include "G4ParametrizedHadronicVertex.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"

G4ParametrizedHadronicVertex::G4ParametrizedHadronicVertex()
{
  theLowEPionPlus = new G4LEPionPlusInelastic;
  theLowEPionMinus = new G4LEPionMinusInelastic;
  theHighEPionPlus = new G4HEPionPlusInelastic;
  theHighEPionMinus = new G4HEPionMinusInelastic;
}

G4ParametrizedHadronicVertex::~G4ParametrizedHadronicVertex()
{}

G4VParticleChange* 
G4ParametrizedHadronicVertex::ApplyYourself(G4Nucleus & theTarget, 
					    const G4Track &thePhoton)
{   
  theTotalResult.Clear();
  theTotalResult.Initialize(thePhoton);

  G4double theKineticEnergy = thePhoton.GetKineticEnergy();
  G4HadFinalState * aR = 0;
  G4HadProjectile thePro(thePhoton);
  if(G4UniformRand() > 0.5)
    {
      if(theKineticEnergy<20*GeV) aR = theLowEPionMinus->ApplyYourself(thePro, theTarget);
      else aR = theHighEPionMinus->ApplyYourself(thePro, theTarget);
    }
  else
    {
      if(theKineticEnergy<20*GeV) aR = theLowEPionPlus->ApplyYourself(thePro, theTarget);
      else aR = theHighEPionPlus->ApplyYourself(thePro, theTarget);
    }
  if(!aR) return &theTotalResult;

  aR->SetTrafoToLab(thePro.GetTrafoToLab());

  G4int nsec = aR->GetNumberOfSecondaries();
  theTotalResult.SetNumberOfSecondaries(nsec);

  G4ThreeVector it(0., 0., 1.);
  G4double what = twopi*G4UniformRand();
  for(G4int i=0; i<nsec; i++)
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
