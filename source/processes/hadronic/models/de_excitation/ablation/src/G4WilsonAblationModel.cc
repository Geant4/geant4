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
// MODULE:              G4WilsonAblationModel.cc
//
// Version:		1.0
// Date:		08/12/2009
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
// 08 December 2009, P R Truscott, QinetiQ Ltd, UK
// Ver 1.0
// Updated as a result of changes in the G4Evaporation classes.  These changes
// affect mostly SelectSecondariesByEvaporation, and now you have variables
// associated with the evaporation model which can be changed:
//    OPTxs to select the inverse cross-section
//    OPTxs = 0      => Dostrovski's parameterization
//    OPTxs = 1 or 2 => Chatterjee's paramaterization
//    OPTxs = 3 or 4 => Kalbach's parameterization
//    useSICB        => use superimposed Coulomb Barrier for inverse cross
//                      sections
// Other problem found with G4Fragment definition using Lorentz vector and
// **G4ParticleDefinition**.  This does not allow A and Z to be defined for the
// fragment for some reason.  Now the fragment is defined more explicitly:
//    G4Fragment *fragment = new G4Fragment(A, Z, lorentzVector);
// to avoid this quirk.  Bug found in SelectSecondariesByDefault: *type is now
// equated to evapType[i] whereas previously it was equated to fragType[i].
// 
// 06 August 2015, A. Ribon, CERN
// Migrated std::exp and std::pow to the faster G4Exp and G4Pow.
//
// 09 June 2017, C. Mancini Terracciano, INFN
// Fixed bug on the initialization of Photon Evaporation model
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include <iomanip>
#include <numeric>

#include "G4WilsonAblationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
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
#include "G4PhotonEvaporation.hh"
#include "G4LorentzVector.hh"
#include "G4VEvaporationChannel.hh"

#include "G4Exp.hh"
#include "G4Pow.hh"


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
  for(G4int i=0; i<200; ++i) { fSig[i] = 0.0; }
//
//
// Set verboseLevel default to no output.
//
  verboseLevel = 0;
  theChannelFactory = new G4EvaporationFactory(new G4PhotonEvaporation());
  theChannels = theChannelFactory->GetChannel();
//
//
// Set defaults for evaporation classes.  These can be overridden by user
// "set" methods.
//
  OPTxs   = 3;
  useSICB = false;
  fragmentVector = 0;
}
////////////////////////////////////////////////////////////////////////////////
//
G4WilsonAblationModel::~G4WilsonAblationModel()
{}

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
  G4int A     = theNucleus.GetA_asInt();
  G4int Z     = theNucleus.GetZ_asInt();
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
// The following lines are no longer accurate given we now treat the final fragment
//  if (verboseLevel >= 2)
//    G4cout <<"Number of nucleons ejected = " <<DAabl <<G4endl;

//
//
// Determine the nuclear fragment from the ablation process by sampling the
// Rudstam equation.
//
  G4int AF = A - DAabl;
  G4int ZF = 0;
  
  if (AF > 0)
  {
    G4Pow* g4calc = G4Pow::GetInstance(); 
    G4double AFd = (G4double) AF;
    G4double R = 11.8 / g4calc->powZ(AF, 0.45);
    G4int minZ = std::max(1, Z - DAabl);
//
//
// Here we define an integral probability distribution based on the Rudstam
// equation assuming a constant AF.
//    
    G4int zmax = std::min(199, Z);
    G4double sum = 0.0;
    for (ZF=minZ; ZF<=zmax; ++ZF)
    {
      sum += G4Exp(-R*g4calc->powA(std::abs(ZF - 0.486*AFd + 3.8E-04*AFd*AFd),1.5));
      fSig[ZF] = sum;
    }
//
//
// Now sample that distribution to determine a value for ZF.
//
    sum *= G4UniformRand();
    for (ZF=minZ; ZF<=zmax; ++ZF) {
      if(sum <= fSig[ZF]) { break; }
    } 
  }
  G4int DZabl = Z - ZF;
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
    }
  }
//
//
// Determine the properties of the final nuclear fragment.  Note that if
// the final fragment is predicted to have a nucleon number of zero, then
// really it's the particle last in the vector evapType which becomes the
// final fragment.  Therefore delete this from the vector if this is the
// case.
//
  G4double massFinalFrag = 0.0;
  if (AF > 0)
    massFinalFrag = G4ParticleTable::GetParticleTable()->GetIonTable()->
      GetIonMass(ZF,AF);
  else
  {
    G4ParticleDefinition *type = evapType[evapType.size()-1];
    AF                         = type->GetBaryonNumber();
    ZF                         = (G4int) (type->GetPDGCharge() + 1.0E-10);
    evapType.erase(evapType.end()-1);
  }
  totalEpost   += massFinalFrag;
//
//
// Provide verbose output on the nuclear fragment if requested.
//
  if (verboseLevel >= 2)
  {
    G4cout <<"Final fragment      A=" <<AF
           <<", Z=" <<ZF
           <<G4endl;
    for (G4int ift=0; ift<nFragTypes; ift++)
    {
      G4ParticleDefinition *type = fragType[ift];
      G4int n                    = std::count(evapType.begin(),evapType.end(),type);
      if (n > 0) 
        G4cout <<"Particle type: " <<std::setw(10) <<type->GetParticleName()
               <<", number of particles emitted = " <<n <<G4endl;
    }
  }
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
    G4Fragment* frag = new G4Fragment(AF, ZF, lorentzVector);
    fragmentVector->push_back(frag);
  }
  delete resultNucleus;
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
    for (iter = fragmentVector->begin(); iter != fragmentVector->end(); iter++)
    {
      if (ie == nEvap)
      {
//        G4cout <<*iter  <<G4endl;
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
  G4Fragment theResidualNucleus = *intermediateNucleus;
  G4bool evaporate = true;
  // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
  while (evaporate && evapType.size() != 0)
  {
//
//
// Here's the cheaky bit.  We're hijacking the G4Evaporation model, in order to
// more accurately sample to kinematics, but the species of the nuclear
// fragments will be the ones of our choosing as above.
//
    std::vector <G4VEvaporationChannel*>  theChannels1;
    theChannels1.clear();
    std::vector <G4VEvaporationChannel*>::iterator i;
    VectorOfFragmentTypes::iterator iter;
    std::vector <VectorOfFragmentTypes::iterator> iters;
    iters.clear();
    iter = std::find(evapType.begin(), evapType.end(), G4Alpha::Alpha());
    if (iter != evapType.end())
    {
      theChannels1.push_back(new G4AlphaEvaporationChannel);
      i = theChannels1.end() - 1;
      (*i)->SetOPTxs(OPTxs);
      (*i)->UseSICB(useSICB);
//      (*i)->Initialize(theResidualNucleus);
      iters.push_back(iter);
    }
    iter = std::find(evapType.begin(), evapType.end(), G4He3::He3());
    if (iter != evapType.end())
    {
      theChannels1.push_back(new G4He3EvaporationChannel);
      i = theChannels1.end() - 1;
      (*i)->SetOPTxs(OPTxs);
      (*i)->UseSICB(useSICB);
//      (*i)->Initialize(theResidualNucleus);
      iters.push_back(iter);
    }
    iter = std::find(evapType.begin(), evapType.end(), G4Triton::Triton());
    if (iter != evapType.end())
    {
      theChannels1.push_back(new G4TritonEvaporationChannel);
      i = theChannels1.end() - 1;
      (*i)->SetOPTxs(OPTxs);
      (*i)->UseSICB(useSICB);
//      (*i)->Initialize(theResidualNucleus);
      iters.push_back(iter);
    }
    iter = std::find(evapType.begin(), evapType.end(), G4Deuteron::Deuteron());
    if (iter != evapType.end())
    {
      theChannels1.push_back(new G4DeuteronEvaporationChannel);
      i = theChannels1.end() - 1;
      (*i)->SetOPTxs(OPTxs);
      (*i)->UseSICB(useSICB);
//      (*i)->Initialize(theResidualNucleus);
      iters.push_back(iter);
    }
    iter = std::find(evapType.begin(), evapType.end(), G4Proton::Proton());
    if (iter != evapType.end())
    {
      theChannels1.push_back(new G4ProtonEvaporationChannel);
      i = theChannels1.end() - 1;
      (*i)->SetOPTxs(OPTxs);
      (*i)->UseSICB(useSICB);
//      (*i)->Initialize(theResidualNucleus);
      iters.push_back(iter);
    }
    iter = std::find(evapType.begin(), evapType.end(), G4Neutron::Neutron());
    if (iter != evapType.end())
    {
      theChannels1.push_back(new G4NeutronEvaporationChannel);
      i = theChannels1.end() - 1;
      (*i)->SetOPTxs(OPTxs);
      (*i)->UseSICB(useSICB);
//      (*i)->Initialize(theResidualNucleus);
      iters.push_back(iter);
    }
    G4int nChannels = theChannels1.size();

    G4double totalProb = 0.0;
    G4int ich = 0;
    G4double probEvapType[6] = {0.0};
    std::vector<G4VEvaporationChannel*>::iterator iterEv;
    for (iterEv=theChannels1.begin(); iterEv!=theChannels1.end(); iterEv++) {
      totalProb += (*iterEv)->GetEmissionProbability(intermediateNucleus);
      probEvapType[ich] = totalProb;
      ++ich;
    }
    if (totalProb > 0.0) {
//
//
// The emission probability for at least one of the evaporation channels is
// positive, therefore work out which one should be selected and decay
// the nucleus.
//
      G4double xi = totalProb*G4UniformRand();
      G4int ii     = 0;
      for (ii=0; ii<nChannels; ii++) {
        if (xi < probEvapType[ii]) { break; }
      }
      if (ii >= nChannels) { ii = nChannels - 1; }
      G4FragmentVector *evaporationResult = theChannels1[ii]->
        BreakUpFragment(intermediateNucleus);
      fragmentVector->push_back((*evaporationResult)[0]);
      intermediateNucleus = (*evaporationResult)[1];
      delete evaporationResult;
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
    G4ParticleDefinition *type = evapType[i];
    G4double mass              = type->GetPDGMass();
    G4double e                 = mass + 10.0*eV;
    G4double p                 = std::sqrt(e*e-mass*mass);
    G4double costheta          = 2.0*G4UniformRand() - 1.0;
    G4double sintheta          = std::sqrt((1.0 - costheta)*(1.0 + costheta));
    G4double phi               = twopi * G4UniformRand() * rad;
    G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),costheta);
    G4LorentzVector lorentzVector = G4LorentzVector(direction*p, e);
    lorentzVector.boost(-boost);
// Possibility that the following line is not correctly carrying over A and Z
// from particle definition.  Force values.  PRT 03/12/2009.
//    G4Fragment *fragment          = 
//      new G4Fragment(lorentzVector, type);
    G4int A = type->GetBaryonNumber();
    G4int Z = (G4int) (type->GetPDGCharge() + 1.0E-10);
    G4Fragment *fragment          = 
      new G4Fragment(A, Z, lorentzVector);

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
  G4cout <<" !!! WARNING: This model is not well validation and should not be used for accurate simulation !!!"
         <<G4endl;
  G4cout <<" *****************************************************************"
         <<G4endl;
  G4cout << G4endl;

  return;
}
////////////////////////////////////////////////////////////////////////////////
//
