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
// $Id: G4TheoFSGenerator.cc,v 1.12 2003/06/16 17:06:57 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// G4TheoFSGenerator
#include "G4DynamicParticle.hh"
#include "G4TheoFSGenerator.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"

G4TheoFSGenerator::G4TheoFSGenerator()
{
 theParticleChange = new G4ParticleChange;
}

G4TheoFSGenerator::G4TheoFSGenerator(const G4TheoFSGenerator &) : G4HadronicInteraction()
{
}


G4TheoFSGenerator::~G4TheoFSGenerator()
{
  delete theParticleChange;
}


const G4TheoFSGenerator & G4TheoFSGenerator::operator=(const G4TheoFSGenerator &)
{
  G4Exception("G4CrossSectionBase::operator= meant to not be accessable");
  return *this;
}


int G4TheoFSGenerator::operator==(const G4TheoFSGenerator &) const
{
  return 0;
}

int G4TheoFSGenerator::operator!=(const G4TheoFSGenerator &) const
{
  return 1;
}


G4VParticleChange * G4TheoFSGenerator::ApplyYourself(const G4Track & thePrimary, G4Nucleus &theNucleus)
{
  // init particle change
  theParticleChange->Clear();
  theParticleChange->Initialize(thePrimary);
  theParticleChange->SetStatusChange(fStopAndKill);
  
  // check if models have been registered, and use default, in case this is not true @@
  
  // get result from high energy model
  const G4DynamicParticle * aPart = thePrimary.GetDynamicParticle();
  G4KineticTrackVector * theInitialResult =
               theHighEnergyGenerator->Scatter(theNucleus, *aPart);
  
  G4double predecayEnergy = 0;
  for(size_t hpw=0; hpw<theInitialResult->size(); hpw++)
  {
    predecayEnergy += (*theInitialResult)[hpw]->Get4Momentum().t();
  }
  
  G4ReactionProductVector * theTransportResult = NULL;
  G4int hitCount = 0;
  const std::vector<G4Nucleon *> & they = theHighEnergyGenerator->GetWoundedNucleus()->GetNucleons();
  for(size_t them=0; them<they.size(); them++)
  {
    if(they[them]->AreYouHit()) hitCount ++;
  }
  if(hitCount != theHighEnergyGenerator->GetWoundedNucleus()->GetMassNumber() )
  {
    theTransportResult = 
               theTransport->Propagate(theInitialResult, theHighEnergyGenerator->GetWoundedNucleus());
  }
  else
  {
     theTransportResult = theDecay.Propagate(theInitialResult, theHighEnergyGenerator->GetWoundedNucleus());
  }

  // Fill particle change
  unsigned int i;
  theParticleChange->SetNumberOfSecondaries(theTransportResult->size());
  for(i=0; i<theTransportResult->size(); i++)
  {
    G4DynamicParticle * aNew = 
       new G4DynamicParticle(theTransportResult->operator[](i)->GetDefinition(),
                             theTransportResult->operator[](i)->GetTotalEnergy(),
                             theTransportResult->operator[](i)->GetMomentum());
    G4double newTime = theParticleChange->GetGlobalTime(theTransportResult->operator[](i)->GetFormationTime());
    theParticleChange->AddSecondary(aNew, newTime);
    delete theTransportResult->operator[](i);
  }
  
  // some garbage collection
  delete theTransportResult;
  
  // Done
  return theParticleChange;
}

