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
// This software was developed by Lawrence Livermore National Laboratory.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// 3. The name of the author may not be used to endorse or promote products
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Copyright (c) 2006 The Regents of the University of California.
// All rights reserved.
// UCRL-CODE-224807
//
//
// $Id: G4FissionLibrary.cc 67966 2013-03-13 09:38:38Z gcosmo $
//
// neutron_hp -- source file
// J.M. Verbeke, Jan-2007
// A low energy neutron-induced fission model.
//

#include "G4FissionLibrary.hh"
#include "G4ParticleHPManager.hh"
#include "G4SystemOfUnits.hh"

G4FissionLibrary::G4FissionLibrary()
  : G4ParticleHPFinalState(), theIsotope(0), targetMass(0.0)
{
  hasXsec = false;
  fe=0;
}

G4FissionLibrary::~G4FissionLibrary()
{}

G4ParticleHPFinalState * G4FissionLibrary::New()
{
  G4FissionLibrary * theNew = new G4FissionLibrary;
  return theNew;
}

//void G4FissionLibrary::Init (G4double A, G4double Z, G4String & dirName, G4String &)
void G4FissionLibrary::Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String &, G4ParticleDefinition*)
{
  G4String tString = "/FS/";
  G4bool dbool;
  theIsotope = static_cast<G4int>(1000*Z+A);
  G4ParticleHPDataUsed aFile = theNames.GetName(static_cast<G4int>(A), static_cast<G4int>(Z), M, dirName, tString, dbool);
  G4String filename = aFile.GetName();

  if(!dbool)
  {
    hasAnyData = false;
    hasFSData = false;
    hasXsec = false;
    return;
  }
  //std::ifstream theData(filename, std::ios::in);
  std::istringstream theData(std::ios::in);
  G4ParticleHPManager::GetInstance()->GetDataStream(filename,theData);

  // here it comes
  G4int infoType, dataType;
  hasFSData = false;
  while (theData >> infoType) // Loop checking, 11.03.2015, T. Koi
  {
    hasFSData = true;
    theData >> dataType;
    switch(infoType)
    {
      case 1:
        if(dataType==4) theNeutronAngularDis.Init(theData);
        if(dataType==5) thePromptNeutronEnDis.Init(theData);
        if(dataType==12) theFinalStatePhotons.InitMean(theData);
        if(dataType==14) theFinalStatePhotons.InitAngular(theData);
        if(dataType==15) theFinalStatePhotons.InitEnergies(theData);
        break;
      case 2:
        if(dataType==1) theFinalStateNeutrons.InitMean(theData);
        break;
      case 3:
        if(dataType==1) theFinalStateNeutrons.InitDelayed(theData);
        if(dataType==5) theDelayedNeutronEnDis.Init(theData);
        break;
      case 4:
        if(dataType==1) theFinalStateNeutrons.InitPrompt(theData);
        break;
      case 5:
        if(dataType==1) theEnergyRelease.Init(theData);
        break;
      default:
        G4cout << "G4FissionLibrary::Init: unknown data type"<<dataType<<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, "G4FissionLibrary::Init: unknown data type");
        break;
    }
  }
  targetMass = theFinalStateNeutrons.GetTargetMass();
  //theData.close();
}

G4HadFinalState* G4FissionLibrary::ApplyYourself(const G4HadProjectile & theTrack)
{  

  if ( theResult.Get() == NULL ) theResult.Put( new G4HadFinalState );
  theResult.Get()->Clear();

  // prepare neutron
  G4double eKinetic = theTrack.GetKineticEnergy();
  const G4HadProjectile* incidentParticle = &theTrack;
  G4ReactionProduct theNeutron(incidentParticle->GetDefinition() );
  theNeutron.SetMomentum(incidentParticle->Get4Momentum().vect() );
  theNeutron.SetKineticEnergy(eKinetic);

  // prepare target
  G4Nucleus aNucleus;
  G4ReactionProduct theTarget; 
  G4ThreeVector neuVelo = (1./incidentParticle->GetDefinition()->GetPDGMass())*theNeutron.GetMomentum();
  theTarget = aNucleus.GetBiasedThermalNucleus( targetMass, neuVelo, theTrack.GetMaterial()->GetTemperature());

  // set neutron and target in the FS classes 
  //theNeutronAngularDis.SetNeutron(theNeutron);
  theNeutronAngularDis.SetProjectileRP(theNeutron);
  theNeutronAngularDis.SetTarget(theTarget);

  // boost to target rest system
  theNeutron.Lorentz(theNeutron, -1*theTarget);

  eKinetic = theNeutron.GetKineticEnergy();    

  // dice neutron and gamma multiplicities, energies and momenta in Lab. @@
  // no energy conservation on an event-to-event basis. we rely on the data to be ok. @@
  // also for mean, we rely on the consistency of the data. @@

  G4int nPrompt=0, gPrompt=0;
  SampleMult(theTrack, &nPrompt, &gPrompt, eKinetic);

  // Build neutrons and add them to dynamic particle vector
  G4double momentum;
  for(G4int i=0; i<nPrompt; i++)
  {
    G4DynamicParticle * it = new G4DynamicParticle;
    it->SetDefinition(G4Neutron::Neutron());
    it->SetKineticEnergy(fe->getNeutronEnergy(i)*MeV);
    momentum = it->GetTotalMomentum();
    G4ThreeVector temp(momentum*fe->getNeutronDircosu(i), 
                       momentum*fe->getNeutronDircosv(i), 
                       momentum*fe->getNeutronDircosw(i));
    it->SetMomentum( temp );
//    it->SetGlobalTime(fe->getNeutronAge(i)*second);
    theResult.Get()->AddSecondary(it);
//    G4cout <<"G4FissionLibrary::ApplyYourself: energy of prompt neutron " << i << " = " << it->GetKineticEnergy()<<G4endl;
  }

  // Build gammas, lorentz transform them, and add them to dynamic particle vector
  for(G4int i=0; i<gPrompt; i++)
  {
    G4ReactionProduct * thePhoton = new G4ReactionProduct;
    thePhoton->SetDefinition(G4Gamma::Gamma());
    thePhoton->SetKineticEnergy(fe->getPhotonEnergy(i)*MeV);
    momentum = thePhoton->GetTotalMomentum();
    G4ThreeVector temp(momentum*fe->getPhotonDircosu(i), 
                       momentum*fe->getPhotonDircosv(i), 
                       momentum*fe->getPhotonDircosw(i));
    thePhoton->SetMomentum( temp );
    thePhoton->Lorentz(*thePhoton, -1.*theTarget);
    
    G4DynamicParticle * it = new G4DynamicParticle;
    it->SetDefinition(thePhoton->GetDefinition());
    it->SetMomentum(thePhoton->GetMomentum());
//    it->SetGlobalTime(fe->getPhotonAge(i)*second);
//    G4cout <<"G4FissionLibrary::ApplyYourself: energy of prompt photon " << i << " = " << it->GetKineticEnergy()<<G4endl;
    theResult.Get()->AddSecondary(it);
    delete thePhoton;  
  }
//  G4cout <<"G4FissionLibrary::ApplyYourself: Number of secondaries = "<<theResult.GetNumberOfSecondaries()<< G4endl;
//  G4cout <<"G4FissionLibrary::ApplyYourself: Number of induced prompt neutron = "<<nPrompt<<G4endl;
//  G4cout <<"G4FissionLibrary::ApplyYourself: Number of induced prompt photons = "<<gPrompt<<G4endl;

  // finally deal with local energy depositions.
  G4double eDepByFragments = theEnergyRelease.GetFragmentKinetic();
  theResult.Get()->SetLocalEnergyDeposit(eDepByFragments);
//   G4cout << "G4FissionLibrary::local energy deposit" << eDepByFragments<<G4endl;
  // clean up the primary neutron
  theResult.Get()->SetStatusChange(stopAndKill);
  return theResult.Get();
}

void G4FissionLibrary::SampleMult(const G4HadProjectile & theTrack, G4int* nPrompt,
                                   G4int* gPrompt, G4double eKinetic)
{
   G4double promptNeutronMulti = 0;
   promptNeutronMulti = theFinalStateNeutrons.GetPrompt(eKinetic); // prompt nubar from Geant
   G4double delayedNeutronMulti = 0;
   delayedNeutronMulti = theFinalStateNeutrons.GetDelayed(eKinetic); // delayed nubar from Geant

   G4double time = theTrack.GetGlobalTime()/second;
   G4double totalNeutronMulti = theFinalStateNeutrons.GetMean(eKinetic);
   if(delayedNeutronMulti==0&&promptNeutronMulti==0) {
     // no data for prompt and delayed neutrons in Geant
     // but there is perhaps data for the total neutron multiplicity, in which case 
     // we use it for prompt neutron emission
     if (fe != 0) delete fe;
     fe = new G4fissionEvent(theIsotope, time, totalNeutronMulti, eKinetic);
   } else {
     // prompt nubar != 0 || delayed nubar != 0
     if (fe != 0) delete fe;
     fe = new G4fissionEvent(theIsotope, time, promptNeutronMulti, eKinetic);
   }
   *nPrompt = fe->getNeutronNu();
   if (*nPrompt == -1) *nPrompt = 0; // the fission library libFission.a has no data for neutrons
   *gPrompt = fe->getPhotonNu();
   if (*gPrompt == -1) *gPrompt = 0; // the fission library libFission.a has no data for gammas
}

