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
// $Id: G4PreCompoundFragmentVector.cc,v 1.3.2.1 2001/06/28 19:13:35 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 

#include "G4PreCompoundFragmentVector.hh"

#include "G4PreCompoundNeutron.hh"
#include "G4PreCompoundProton.hh"
#include "G4PreCompoundDeuteron.hh"
#include "G4PreCompoundTriton.hh"
#include "G4PreCompoundHe3.hh"
#include "G4PreCompoundAlpha.hh"

G4PreCompoundFragmentVector::G4PreCompoundFragmentVector() :
  TotalEmissionProbability(0.0)
{
  //  theChannels.reserve(6)
  // neutron 
    //  theChannels.insert(new G4PreCompoundNeutron());
  theChannels.push_back(new G4PreCompoundNeutron());
  // proton
  //  theChannels.insert(new G4PreCompoundProton());
  theChannels.push_back(new G4PreCompoundProton());
  // deuterium
  //  theChannels.insert(new G4PreCompoundDeuteron());
  theChannels.push_back(new G4PreCompoundDeuteron());
  // triton
  //  theChannels.insert(new G4PreCompoundTriton());
  theChannels.push_back(new G4PreCompoundTriton());
  // helium3
  //  theChannels.insert(new G4PreCompoundHe3());
  theChannels.push_back(new G4PreCompoundHe3());
  // alpha
  //  theChannels.insert(new G4PreCompoundAlpha());
  theChannels.push_back(new G4PreCompoundAlpha());
}


G4PreCompoundFragmentVector::~G4PreCompoundFragmentVector()
{
  //  theChannels.clearAndDestroy();
  for (G4std::vector<G4VPreCompoundFragment*>::iterator i=theChannels.begin(); 
       i != theChannels.end(); i++) delete *i;
  theChannels.clear();
}

const G4PreCompoundFragmentVector & G4PreCompoundFragmentVector::operator=(const G4PreCompoundFragmentVector &right)
{
  G4Exception("G4PreCompoundFragmentVector::operator= meant to not be accessable");
  return *this;
}


G4bool G4PreCompoundFragmentVector::operator==(const G4PreCompoundFragmentVector &right) const
{
  return false;
}

G4bool G4PreCompoundFragmentVector::operator!=(const G4PreCompoundFragmentVector &right) const
{
  return true;
}



G4double G4PreCompoundFragmentVector::CalculateProbabilities(const G4Fragment & aFragment)
{
  TotalEmissionProbability = 0.0;
  G4std::vector<G4VPreCompoundFragment*>::iterator aChannel; 
  for (aChannel=theChannels.begin(); aChannel != theChannels.end(); aChannel++) {
    (*aChannel)->CalcExcitonLevelDensityRatios(aFragment.GetNumberOfExcitons(),
					       aFragment.GetNumberOfParticles());	
    (*aChannel)->CalcCondensationProbability(aFragment.GetA());
    // Calculate emission probailities
    if (aFragment.GetNumberOfParticles() <= (*aChannel)->GetA()-0.01) {
      // if number of particles less than a fragment atomic number 
      // set probability to emit a fragment 0
      (*aChannel)->SetEmissionProbability(0.0);
    } else if (aFragment.GetNumberOfExcitons() <= (*aChannel)->GetA()+0.01 && 
	       aFragment.GetNumberOfExcitons() != 1) {
      (*aChannel)->SetEmissionProbability(0.0);
    } else if (aFragment.GetNumberOfCharged() <= (*aChannel)->GetZ()-0.01) {
      // if number of charged particles (protons) is less than charge of fragment
      // set probability to emit a fragment 0
      (*aChannel)->SetEmissionProbability(0.0);
    } else if ((*aChannel)->GetMaximalKineticEnergy() <= 0.0) {
      // if the energy threshold for emitted fragment is less or equal 0
      // set probability to emit a fragment 0
      (*aChannel)->SetEmissionProbability(0.0);
    } else {
      // Compute total (integrated over kinetic energy) emission 
      // probability of a fragment and
      // Summing channel emission probabilities
      TotalEmissionProbability += (*aChannel)->CalcEmissionProbability(aFragment);
    }
  }
  return TotalEmissionProbability;
}


G4VPreCompoundFragment * G4PreCompoundFragmentVector::ChooseFragment(void)
{
  const G4int NumOfFrags = theChannels.size();
  G4double * running = new G4double[NumOfFrags];
  running[0] = (*theChannels.begin())->GetEmissionProbability();
  //  G4std::vector<G4VPreCompoundFragment*>::iterator aChannel; 
  G4int i;
  for (i = 1; i < NumOfFrags; i++) {
    running[i]=running[i-1]+theChannels[i]->GetEmissionProbability();
  }
	
  // Choose an emission channel
  G4double aChannel = G4UniformRand()*TotalEmissionProbability;
  G4int ChosenChannel = -1;
  for (i = 0; i < NumOfFrags; i++) {
    if (aChannel <= running[i]) {
      ChosenChannel = i;
      break;
    }
  }
  delete [] running;
  if (ChosenChannel < 0) 
    G4Exception("G4PreCompoundFragmentVector::ChooseFragment: I can't determine a channel");
	
  return theChannels[ChosenChannel];
}


