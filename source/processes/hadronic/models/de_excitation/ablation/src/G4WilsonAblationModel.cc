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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the	      *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 17191/03/NL/LvH (Aurora Programme). 		      *
// *                                                                  *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4WilsonAblationModel.cc
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
// 6 October 2003, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// 15 March 2004, P R Truscott, QinetiQ Ltd, UK
// Beta release
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "G4WilsonAblationModel.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Alpha.hh"
#include "G4He3.hh"
#include "G4Triton.hh"
#include "G4Deuteron.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4AlphaEvaporationChannel.hh"
#include "G4He3EvaporationChannel.hh"
#include "G4TritonEvaporationChannel.hh"
#include "G4DeuteronEvaporationChannel.hh"
#include "G4ProtonEvaporationChannel.hh"
#include "G4NeutronEvaporationChannel.hh"
#include "G4LorentzVector.hh"
#include "G4VEvaporationChannel.hh"

#include <iomanip>
#include <numeric>
////////////////////////////////////////////////////////////////////////////////
//
G4WilsonAblationModel::G4WilsonAblationModel()
{
//
//
// Send message to stdout to advise that the G4Abrasion model is being used.
//
  PrintWelcomeMessage();
//
//
// Set the default verbose level to 0 - no output.
//
  verboseLevel = 0;  
//
//
// Set the binding energy per nucleon .... did I mention that this is a crude
// model for nuclear de-excitation?
//
  B = 10.0 * MeV;
//
//
// It is possuble to switch off secondary particle production (other than the
// final nuclear fragment).  The default is on.
//
  produceSecondaries = true;
//
//
// Now we need to define the decay modes.  We're using the G4Evaporation model
// to help determine the kinematics of the decay.
//
  nFragTypes  = 6;
  fragType[0] = G4Alpha::Alpha();
  fragType[1] = G4He3::He3();
  fragType[2] = G4Triton::Triton();
  fragType[3] = G4Deuteron::Deuteron();
  fragType[4] = G4Proton::Proton();
  fragType[5] = G4Neutron::Neutron();
//
//
// Set verboseLevel default to no output.
//
  verboseLevel = 0;
}
////////////////////////////////////////////////////////////////////////////////
//
G4WilsonAblationModel::~G4WilsonAblationModel()
{;}
////////////////////////////////////////////////////////////////////////////////
//
G4FragmentVector *G4WilsonAblationModel::BreakItUp
  (const G4Fragment &theNucleus)
{
//
//
// Initilise the pointer to the G4FragmentVector used to return the information
// about the breakup.
//
  fragmentVector = new G4FragmentVector;
  fragmentVector->clear();
//
//
// Get the A, Z and excitation of the nucleus.
//
  G4int A     = (G4int) theNucleus.GetA();
  G4int Z     = (G4int) theNucleus.GetZ();
  G4double ex = theNucleus.GetExcitationEnergy();
  if (verboseLevel >= 2)
  {
    G4cout <<"oooooooooooooooooooooooooooooooooooooooo"
           <<"oooooooooooooooooooooooooooooooooooooooo"
           <<G4endl;
    G4cout.precision(6);
    G4cout <<"IN G4WilsonAblationModel" <<G4endl;
    G4cout <<"Initial prefragment A=" <<A
           <<", Z=" <<Z
           <<", excitation energy = " <<ex/MeV <<" MeV"
           <<G4endl; 
  }
//
//
// Check that there is a nucleus to speak of.  It's possible there isn't one
// or its just a proton or neutron.  In either case, the excitation energy
// (from the Lorentz vector) is not used.
//
  if (A == 0)
  {
    if (verboseLevel >= 2)
    {
      G4cout <<"No nucleus to decay" <<G4endl;
      G4cout <<"oooooooooooooooooooooooooooooooooooooooo"
             <<"oooooooooooooooooooooooooooooooooooooooo"
             <<G4endl;
    }
    return fragmentVector;
  }
  else if (A == 1)
  {
    G4LorentzVector lorentzVector = theNucleus.GetMomentum();
    lorentzVector.setE(lorentzVector.e()-ex+10.0*eV);
    if (Z == 0)
    {
      G4Fragment *fragment = new G4Fragment(lorentzVector,G4Neutron::Neutron());
      fragmentVector->push_back(fragment);
    }
    else
    {
      G4Fragment *fragment = new G4Fragment(lorentzVector,G4Proton::Proton());
      fragmentVector->push_back(fragment);
    }
    if (verboseLevel >= 2)
    {
      G4cout <<"Final fragment is in fact only a nucleon) :" <<G4endl;
      G4cout <<(*fragmentVector)[0] <<G4endl;
      G4cout <<"oooooooooooooooooooooooooooooooooooooooo"
             <<"oooooooooooooooooooooooooooooooooooooooo"
             <<G4endl;
    }
    return fragmentVector;
  }
//
//
// Then the number of nucleons ablated (either as nucleons or light nuclear
// fragments) is based on a simple argument for the binding energy per nucleon.
//
  G4int DAabl = (G4int) (ex / B);
  if (DAabl > A) DAabl = A;
  if (verboseLevel >= 2)
    G4cout <<"Number of nucleons ejected = " <<DAabl <<G4endl;

//
//
// Determine the nuclear fragment from the ablation process by sampling the
// Rudstam equation.
//
  G4int AF = A - DAabl;
  G4int ZF = 0;
  if (AF > 0)
  {
    G4double AFd = static_cast<G4double>(AF);
    G4double R = 11.8 / std::pow(AFd, 0.45);
    G4int minZ = Z - DAabl;
    if (minZ <= 0) minZ = 1;
//
//
// Here we define an integral probability distribution based on the Rudstam
// equation assuming a constant AF.
//    
    G4double sig[100];
    G4double sum = 0.0;
    for (G4int ii=minZ; ii<= Z; ii++)
    {
      sum   += std::exp(-R*std::pow(std::abs(ii - 0.486*AFd + 3.8E-04*AFd*AFd),1.5));
      sig[ii] = sum;
    }
//
//
// Now sample that distribution to determine a value for ZF.
//
    G4double xi  = G4UniformRand();
    G4int iz     = minZ;
    G4bool found = false;
    while (iz <= Z && !found)
    {
      found = (xi <= sig[iz]/sum);
      if (!found) iz++;
    }
    if (iz > Z)
      ZF = Z;
    else
      ZF = iz;
  }
  G4int DZabl = Z - ZF;
  if (verboseLevel >= 2)
    G4cout <<"Final fragment      A=" <<AF
           <<", Z=" <<ZF
           <<G4endl;
//
//
// Now determine the nucleons or nuclei which have bee ablated.  The preference
// is for the production of alphas, then other nuclei in order of decreasing
// binding energy. The energies assigned to the products of the decay are
// provisional for the moment (the 10eV is just to avoid errors with negative
// excitation energies due to rounding).
//
  G4double totalEpost = 0.0;
  evapType.clear();
  for (G4int ift=0; ift<nFragTypes; ift++)
  {
    G4ParticleDefinition *type = fragType[ift];
    G4double n  = std::floor((G4double) DAabl / type->GetBaryonNumber() + 1.0E-10);
    G4double n1 = 1.0E+10;
    if (fragType[ift]->GetPDGCharge() > 0.0)
      n1 = std::floor((G4double) DZabl / type->GetPDGCharge() + 1.0E-10);
    if (n > n1) n = n1;
    if (n > 0.0)
    {
      G4double mass = type->GetPDGMass();
      for (G4int j=0; j<(G4int) n; j++)
      {
        totalEpost += mass;
        evapType.push_back(type);
      }
      DAabl -= (G4int) (n * type->GetBaryonNumber() + 1.0E-10);
      DZabl -= (G4int) (n * type->GetPDGCharge() + 1.0E-10);
      if (verboseLevel >= 2)
        G4cout <<"Particle type: " <<std::setw(10) <<type->GetParticleName()
               <<", number of particles emitted = " <<n
               <<G4endl;
    }
  }
//
//
// Determine the properties of the final nuclear fragment.
//
  G4double massFinalFrag = 0.0;
  if (AF > 0.0)
    massFinalFrag = G4ParticleTable::GetParticleTable()->GetIonTable()->
      GetIonMass(ZF,AF);
  totalEpost   += massFinalFrag;
//
//
// Add the total energy from the fragment.  Note that the fragment is assumed
// to be de-excited and does not undergo photo-evaporation .... I did mention
// this is a bit of a crude model?
//
  G4double massPreFrag      = theNucleus.GetGroundStateMass();
  G4double totalEpre        = massPreFrag + ex;
  G4double excess           = totalEpre - totalEpost;
//  G4Fragment *resultNucleus(theNucleus);
  G4Fragment *resultNucleus = new G4Fragment(A, Z, theNucleus.GetMomentum());
  G4ThreeVector boost(0.0,0.0,0.0);
  G4int nEvap = 0;
  if (produceSecondaries && evapType.size()>0)
  {
    if (excess > 0.0)
    {
      SelectSecondariesByEvaporation (resultNucleus);
      nEvap = fragmentVector->size();
      boost = resultNucleus->GetMomentum().findBoostToCM();
      if (evapType.size() > 0)
        SelectSecondariesByDefault (boost);
    }
    else
      SelectSecondariesByDefault(G4ThreeVector(0.0,0.0,0.0));
  }
  if (AF > 0)
  {
    G4double mass = G4ParticleTable::GetParticleTable()->GetIonTable()->
      GetIonMass(ZF,AF);
    G4double e    = mass + 10.0*eV;
    G4double p    = std::sqrt(e*e-mass*mass);
    G4ThreeVector direction(0.0,0.0,1.0);
    G4LorentzVector lorentzVector = G4LorentzVector(direction*p, e);
    lorentzVector.boost(-boost);
    *resultNucleus = G4Fragment(AF, ZF, lorentzVector);
    fragmentVector->push_back(resultNucleus);
  }
//
//
// Provide verbose output on the ablation products if requested.
//
  if (verboseLevel >= 2)
  {
    if (nEvap > 0)
    {
      G4cout <<"----------------------" <<G4endl;
      G4cout <<"Evaporated particles :" <<G4endl;
      G4cout <<"----------------------" <<G4endl;
    }
    G4int ie = 0;
    G4FragmentVector::iterator iter;
    for (iter = fragmentVector->begin(); iter != fragmentVector->end(); ++iter)
    {
      if (ie == nEvap)
      {
        G4cout <<*iter  <<G4endl;
        G4cout <<"---------------------------------" <<G4endl;
        G4cout <<"Particles from default emission :" <<G4endl;
        G4cout <<"---------------------------------" <<G4endl;
      }
      G4cout <<*iter <<G4endl;
    }
    G4cout <<"oooooooooooooooooooooooooooooooooooooooo"
           <<"oooooooooooooooooooooooooooooooooooooooo"
           <<G4endl;
  }

  return fragmentVector;    
}
////////////////////////////////////////////////////////////////////////////////
//
void G4WilsonAblationModel::SelectSecondariesByEvaporation
  (G4Fragment *intermediateNucleus)
{
  G4bool evaporate = true;
  while (evaporate && evapType.size() != 0)
  {
//
//
// Here's the cheaky bit.  We're hijacking the G4Evaporation model, in order to
// more accurately sample to kinematics, but the species of the nuclear
// fragments will be the ones of our choosing as above.
//
    std::vector <G4VEvaporationChannel*>  theChannels;
    theChannels.clear();
    VectorOfFragmentTypes::iterator iter;
    std::vector <VectorOfFragmentTypes::iterator> iters;
    iters.clear();
    iter = std::find(evapType.begin(), evapType.end(), G4Alpha::Alpha());
    if (iter != evapType.end())
    {
      theChannels.push_back(new G4AlphaEvaporationChannel);
      iters.push_back(iter);
    }
    iter = std::find(evapType.begin(), evapType.end(), G4He3::He3());
    if (iter != evapType.end())
    {
      theChannels.push_back(new G4He3EvaporationChannel);
      iters.push_back(iter);
    }
    iter = std::find(evapType.begin(), evapType.end(), G4Triton::Triton());
    if (iter != evapType.end())
    {
      theChannels.push_back(new G4TritonEvaporationChannel);
      iters.push_back(iter);
    }
    iter = std::find(evapType.begin(), evapType.end(), G4Deuteron::Deuteron());
    if (iter != evapType.end())
    {
      theChannels.push_back(new G4DeuteronEvaporationChannel);
      iters.push_back(iter);
    }
    iter = std::find(evapType.begin(), evapType.end(), G4Proton::Proton());
    if (iter != evapType.end())
    {
      theChannels.push_back(new G4ProtonEvaporationChannel);
      iters.push_back(iter);
    }
    iter = std::find(evapType.begin(), evapType.end(), G4Neutron::Neutron());
    if (iter != evapType.end())
    {
      theChannels.push_back(new G4NeutronEvaporationChannel);
      iters.push_back(iter);
    }
    G4int nChannels = theChannels.size();

    std::vector<G4VEvaporationChannel*>::iterator iterEv;
    for (iterEv=theChannels.begin(); iterEv!=theChannels.end(); iterEv++)
      (*iterEv)->Initialize(*intermediateNucleus);
    G4double totalProb = std::accumulate(theChannels.begin(),
      theChannels.end(), 0.0, SumProbabilities());
    if (totalProb > 0.0)
    {
//
//
// The emission probability for at least one of the evaporation channels is
// positive, therefore work out which one should be selected and decay
// the nucleus.
//
      G4double totalProb1      = 0.0;
      G4double probEvapType[6] = {0.0};
      for (G4int ich=0; ich<nChannels; ich++)
      {
        totalProb1      += theChannels[ich]->GetEmissionProbability();
        probEvapType[ich]  = totalProb1 / totalProb;
      }
      G4double xi = G4UniformRand();
      G4int i     = 0;
      for (i=0; i<nChannels; i++)
        if (xi < probEvapType[i]) break;
      if (i > nChannels) i = nChannels - 1;
      G4FragmentVector *evaporationResult = theChannels[i]->
        BreakUp(*intermediateNucleus);
      fragmentVector->push_back((*evaporationResult)[0]);
      *intermediateNucleus = *(*evaporationResult)[1];
      delete evaporationResult->back();
      delete evaporationResult;
      evapType.erase(iters[i]);
    }
    else
    {
//
//
// Probability for further evaporation is nil so have to escape from this
// routine and set the energies of the secondaries to 10eV.
//
      evaporate = false;
    }
  }
  
  return;
}
////////////////////////////////////////////////////////////////////////////////
//
void G4WilsonAblationModel::SelectSecondariesByDefault (G4ThreeVector boost)
{
  for (unsigned i=0; i<evapType.size(); i++)
  {
    G4ParticleDefinition *type = fragType[i];
    G4double mass              = type->GetPDGMass();
    G4double e                 = mass + 10.0*eV;
    G4double p                 = std::sqrt(e*e-mass*mass);
    G4double costheta          = 2.0*G4UniformRand() - 1.0;
    G4double sintheta          = std::sqrt((1.0 - costheta)*(1.0 + costheta));
    G4double phi               = twopi * G4UniformRand() * rad;
    G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),costheta);
    G4LorentzVector lorentzVector = G4LorentzVector(direction*p, e);
    lorentzVector.boost(-boost);
    G4Fragment *fragment          = 
      new G4Fragment(lorentzVector, type);
    fragmentVector->push_back(fragment);
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void G4WilsonAblationModel::PrintWelcomeMessage ()
{
  G4cout <<G4endl;
  G4cout <<" *****************************************************************"
         <<G4endl;
  G4cout <<" Nuclear ablation model for nuclear-nuclear interactions activated"
         <<G4endl;
  G4cout <<" (Written by QinetiQ Ltd for the European Space Agency)"
         <<G4endl;
  G4cout <<" *****************************************************************"
         <<G4endl;
  G4cout << G4endl;

  return;
}
////////////////////////////////////////////////////////////////////////////////
//
