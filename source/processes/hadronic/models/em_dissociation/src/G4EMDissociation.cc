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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 17191/03/NL/LvH (Aurora Programme).                     *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		G4EMDissociation.cc
//
// Version:		B.1
// Date:		15/04/04
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		17191/03/NL/LvH
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 17 October 2003, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// 15 March 2004, P R Truscott, QinetiQ Ltd, UK
// Beta release
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "G4EMDissociation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4LorentzVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4EMDissociationCrossSection.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4IonTable.hh"
#include "G4DecayProducts.hh"
#include "G4DynamicParticle.hh"
#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"
#include "Randomize.hh"
#include "globals.hh"

G4EMDissociation::G4EMDissociation():G4HadronicInteraction("EMDissociation") {

  // Send message to stdout to advise that the G4EMDissociation model is being
  // used.
  PrintWelcomeMessage();

  // No de-excitation handler has been supplied - define the default handler.
  theExcitationHandler            = new G4ExcitationHandler;
  theExcitationHandler->SetMinEForMultiFrag(5.0*MeV);
  handlerDefinedInternally = true;

  // This EM dissociation model needs access to the cross-sections held in
  // G4EMDissociationCrossSection.
  dissociationCrossSection = new G4EMDissociationCrossSection;
  thePhotonSpectrum = new G4EMDissociationSpectrum;

  // Set the minimum and maximum range for the model (despite nomanclature, this
  // is in energy per nucleon number).    
  SetMinEnergy(100.0*MeV);
  SetMaxEnergy(500.0*GeV);

  // Set the default verbose level to 0 - no output.
  verboseLevel = 0;
}

G4EMDissociation::G4EMDissociation (G4ExcitationHandler *aExcitationHandler)
{
  // Send message to stdout to advise that the G4EMDissociation model is being
  // used.
  PrintWelcomeMessage();
  
  theExcitationHandler     = aExcitationHandler;
  handlerDefinedInternally = false;

  // This EM dissociation model needs access to the cross-sections held in
  // G4EMDissociationCrossSection.
  dissociationCrossSection = new G4EMDissociationCrossSection;
  thePhotonSpectrum = new G4EMDissociationSpectrum;

  // Set the minimum and maximum range for the model (despite nomanclature, this
  // is in energy per nucleon number)    
  SetMinEnergy(100.0*MeV);
  SetMaxEnergy(500.0*GeV);
  verboseLevel = 0;
}


G4EMDissociation::~G4EMDissociation() {
  if (handlerDefinedInternally) delete theExcitationHandler;
  // delete dissociationCrossSection;
  // Cross section deleted by G4CrossSectionRegistry; don't do it here
  // Bug reported by Gong Ding in Bug Report #1339
  delete thePhotonSpectrum;
}


G4HadFinalState *G4EMDissociation::ApplyYourself
  (const G4HadProjectile &theTrack, G4Nucleus &theTarget)
{
  // The secondaries will be returned in G4HadFinalState &theParticleChange -
  // initialise this.

  theParticleChange.Clear();
  theParticleChange.SetStatusChange(stopAndKill);

  // Get relevant information about the projectile and target (A, Z) and
  // energy/nuc, momentum, velocity, Lorentz factor and rest-mass of the
  // projectile.

  const G4ParticleDefinition *definitionP = theTrack.GetDefinition();
  const G4double AP  = definitionP->GetBaryonNumber();
  const G4double ZP  = definitionP->GetPDGCharge();
  G4LorentzVector pP = theTrack.Get4Momentum();
  G4double E         = theTrack.GetKineticEnergy()/AP;
  G4double MP        = theTrack.GetTotalEnergy() - E*AP;
  G4double b         = pP.beta();
  G4double AT        = theTarget.GetA_asInt();
  G4double ZT        = theTarget.GetZ_asInt();
  G4double MT        = G4NucleiProperties::GetNuclearMass(AT,ZT);

  // Depending upon the verbosity level, output the initial information on the
  // projectile and target
  if (verboseLevel >= 2) {
    G4cout.precision(6);
    G4cout <<"########################################"
           <<"########################################"
           <<G4endl;
    G4cout <<"IN G4EMDissociation" <<G4endl;
    G4cout <<"Initial projectile A=" <<AP 
           <<", Z=" <<ZP
           <<G4endl; 
    G4cout <<"Initial target     A=" <<AT
           <<", Z=" <<ZT
           <<G4endl;
    G4cout <<"Projectile momentum and Energy/nuc = " <<pP <<" ," <<E <<G4endl;
  }

  // Initialise the variables which will be used with the phase-space decay and
  // to boost the secondaries from the interaction.
  
  G4ParticleDefinition *typeNucleon  = NULL;
  G4ParticleDefinition *typeDaughter = NULL;
  G4double Eg                        = 0.0;
  G4double mass                      = 0.0;
  G4ThreeVector boost = G4ThreeVector(0.0, 0.0, 0.0);

  // Determine the cross-sections at the giant dipole and giant quadrupole
  // resonance energies for the projectile and then target.  The information is
  // initially provided in the G4PhysicsFreeVector individually for the E1
  // and E2 fields. These are then summed.

  G4double bmin = thePhotonSpectrum->GetClosestApproach(AP, ZP, AT, ZT, b);
  G4PhysicsFreeVector *crossSectionP = dissociationCrossSection->
    GetCrossSectionForProjectile(AP, ZP, AT, ZT, b, bmin);
  G4PhysicsFreeVector *crossSectionT = dissociationCrossSection->
    GetCrossSectionForTarget(AP, ZP, AT, ZT, b, bmin);

  G4double totCrossSectionP = (*crossSectionP)[0]+(*crossSectionP)[1];
  G4double totCrossSectionT = (*crossSectionT)[0]+(*crossSectionT)[1];

  // Now sample whether the interaction involved EM dissociation of the projectile
  // or the target.

  if (G4UniformRand() <
    totCrossSectionP / (totCrossSectionP + totCrossSectionT)) {

    // It was the projectile which underwent EM dissociation.  Define the Lorentz
    // boost to be applied to the secondaries, and sample whether a proton or a
    // neutron was ejected.  Then determine the energy of the virtual gamma ray
    // which passed from the target nucleus ... this will be used to define the
    // excitation of the projectile.

    mass  = MP;
    if (G4UniformRand() < dissociationCrossSection->
      GetWilsonProbabilityForProtonDissociation (AP, ZP))
    {
      if (verboseLevel >= 2)
        G4cout <<"Projectile underwent EM dissociation producing a proton"
               <<G4endl;
      typeNucleon = G4Proton::ProtonDefinition();
      typeDaughter = G4IonTable::GetIonTable()->
      GetIon((G4int) ZP-1, (G4int) AP-1, 0.0);
    }
    else
    {
      if (verboseLevel >= 2)
        G4cout <<"Projectile underwent EM dissociation producing a neutron"
               <<G4endl;
      typeNucleon = G4Neutron::NeutronDefinition();
      typeDaughter = G4IonTable::GetIonTable()->
      GetIon((G4int) ZP, (G4int) AP-1, 0.0);
    }
    if (G4UniformRand() < (*crossSectionP)[0]/totCrossSectionP)
    {
      Eg = crossSectionP->GetLowEdgeEnergy(0);
      if (verboseLevel >= 2)
        G4cout <<"Transition type was E1" <<G4endl;
    }
    else
    {
      Eg = crossSectionP->GetLowEdgeEnergy(1);
      if (verboseLevel >= 2)
        G4cout <<"Transition type was E2" <<G4endl;
    }

    // We need to define a Lorentz vector with the original momentum, but total
    // energy includes the projectile and virtual gamma.  This is then used
    // to calculate the boost required for the secondaries.

    pP.setE(pP.e()+Eg);
    boost = pP.findBoostToCM();
  }
  else
  {
    // It was the target which underwent EM dissociation.  Sample whether a
    // proton or a neutron was ejected.  Then determine the energy of the virtual 
    // gamma ray which passed from the projectile nucleus ... this will be used to
    // define the excitation of the target.

    mass = MT;
    if (G4UniformRand() < dissociationCrossSection->
      GetWilsonProbabilityForProtonDissociation (AT, ZT))
    {
      if (verboseLevel >= 2)
        G4cout <<"Target underwent EM dissociation producing a proton"
               <<G4endl;
      typeNucleon = G4Proton::ProtonDefinition();
      typeDaughter = G4IonTable::GetIonTable()->
      GetIon((G4int) ZT-1, (G4int) AT-1, 0.0);
    }
    else
    {
      if (verboseLevel >= 2)
        G4cout <<"Target underwent EM dissociation producing a neutron"
               <<G4endl;
      typeNucleon = G4Neutron::NeutronDefinition();
      typeDaughter = G4IonTable::GetIonTable()->
      GetIon((G4int) ZT, (G4int) AT-1, 0.0);
    }
    if (G4UniformRand() < (*crossSectionT)[0]/totCrossSectionT)
    {
      Eg = crossSectionT->GetLowEdgeEnergy(0);
      if (verboseLevel >= 2)
        G4cout <<"Transition type was E1" <<G4endl;
    }
    else
    {
      Eg = crossSectionT->GetLowEdgeEnergy(1);
      if (verboseLevel >= 2)
        G4cout <<"Transition type was E2" <<G4endl;
    }

    // Add the projectile to theParticleChange, less the energy of the
    // not-so-virtual gamma-ray.  Not that at the moment, no lateral momentum
    // is transferred between the projectile and target nuclei.

    G4ThreeVector v = pP.vect();
    v.setMag(1.0);
    G4DynamicParticle *changedP = new G4DynamicParticle (definitionP, v, E*AP-Eg);
    theParticleChange.AddSecondary (changedP);
    if (verboseLevel >= 2)
    {
      G4cout <<"Projectile change:" <<G4endl;
      changedP->DumpInfo();
    }
  }

  // Perform a two-body decay based on the restmass energy of the parent and
  // gamma-ray, and the masses of the daughters. In the frame of reference of
  // the nucles, the angular distribution is sampled isotropically, but the
  // the nucleon and secondary nucleus are boosted if they've come from the
  // projectile.

  G4double e  = mass + Eg;
  G4double mass1 = typeNucleon->GetPDGMass();
  G4double mass2 = typeDaughter->GetPDGMass();
  G4double pp = (e+mass1+mass2)*(e+mass1-mass2)*
                (e-mass1+mass2)*(e-mass1-mass2)/(4.0*e*e);
  if (pp < 0.0) {
    pp = 1.0*eV;
//    if (verboseLevel >`= 1)
//    {
//      G4cout <<"IN G4EMDissociation::ApplyYoursef" <<G4endl;
//      G4cout <<"Error in mass of secondaries compared with primary:" <<G4endl;
//      G4cout <<"Rest mass of primary      = " <<mass <<" MeV" <<G4endl;
//      G4cout <<"Virtual gamma energy      = " <<Eg   <<" MeV" <<G4endl;
//      G4cout <<"Rest mass of secondary #1 = " <<mass1   <<" MeV" <<G4endl;
//      G4cout <<"Rest mass of secondary #2 = " <<mass2   <<" MeV" <<G4endl;
//    }
  }
  else
    pp = std::sqrt(pp);
  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
  G4double phi      = 2.0*pi*G4UniformRand()*rad;
  G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),costheta);
  G4DynamicParticle *dynamicNucleon =
    new G4DynamicParticle(typeNucleon, direction*pp);
  dynamicNucleon->Set4Momentum(dynamicNucleon->Get4Momentum().boost(-boost));
  G4DynamicParticle *dynamicDaughter =
    new G4DynamicParticle(typeDaughter, -direction*pp);
  dynamicDaughter->Set4Momentum(dynamicDaughter->Get4Momentum().boost(-boost));

  // The "decay" products have to be transferred to the G4HadFinalState object.
  // Furthermore, the residual nucleus should be de-excited.

  theParticleChange.AddSecondary (dynamicNucleon);
  if (verboseLevel >= 2) {
    G4cout <<"Nucleon from the EMD process:" <<G4endl;
    dynamicNucleon->DumpInfo();
  }

  G4Fragment* theFragment = new
    G4Fragment((G4int) typeDaughter->GetBaryonNumber(),
    (G4int) typeDaughter->GetPDGCharge(), dynamicDaughter->Get4Momentum());

  if (verboseLevel >= 2) {
    G4cout <<"Dynamic properties of the prefragment:" <<G4endl;
    G4cout.precision(6);
    dynamicDaughter->DumpInfo();
    G4cout <<"Nuclear properties of the prefragment:" <<G4endl;
    G4cout <<theFragment <<G4endl;
  }

  G4ReactionProductVector* products =
                      theExcitationHandler->BreakItUp(*theFragment);
  delete theFragment;
  theFragment = NULL;
  
  G4DynamicParticle* secondary = 0;
  G4ReactionProductVector::iterator iter;
  for (iter = products->begin(); iter != products->end(); ++iter) {
    secondary = new G4DynamicParticle((*iter)->GetDefinition(),
    (*iter)->GetTotalEnergy(), (*iter)->GetMomentum());
    theParticleChange.AddSecondary (secondary);
  }
  delete products;

  delete crossSectionP;
  delete crossSectionT;

  if (verboseLevel >= 2)
    G4cout <<"########################################"
           <<"########################################"
           <<G4endl;
 
  return &theParticleChange;
}


void G4EMDissociation::PrintWelcomeMessage ()
{
  G4cout <<G4endl;
  G4cout <<" ****************************************************************"
         <<G4endl;
  G4cout <<" EM dissociation model for nuclear-nuclear interactions activated"
         <<G4endl;
  G4cout <<" (Written by QinetiQ Ltd for the European Space Agency)"
         <<G4endl;
  G4cout <<" ****************************************************************"
         <<G4endl;
  G4cout << G4endl;

  return;
}

