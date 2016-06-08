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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PreCompoundFragmentVector.cc,v 1.9 2002/06/17 15:29:30 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
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
    // neutron 
    theChannels.push_back(new G4PreCompoundNeutron());
    // proton
    theChannels.push_back(new G4PreCompoundProton());
    // deuterium
    theChannels.push_back(new G4PreCompoundDeuteron());
    // triton
    theChannels.push_back(new G4PreCompoundTriton());
    // helium3
    theChannels.push_back(new G4PreCompoundHe3());
    // alpha
    theChannels.push_back(new G4PreCompoundAlpha());
}


G4PreCompoundFragmentVector::~G4PreCompoundFragmentVector()
{
  G4std::for_each(theChannels.begin(), theChannels.end(), 
		  DeleteFragment());
  theChannels.clear();
}

const G4PreCompoundFragmentVector & 
G4PreCompoundFragmentVector::
operator=(const G4PreCompoundFragmentVector &right)
{
    G4Exception("G4PreCompoundFragmentVector::operator= meant to not be accessable");
    return *this;
}


G4bool G4PreCompoundFragmentVector::
operator==(const G4PreCompoundFragmentVector &right) const
{
    return false;
}

G4bool G4PreCompoundFragmentVector::
operator!=(const G4PreCompoundFragmentVector &right) const
{
    return true;
}



G4double G4PreCompoundFragmentVector::
CalculateProbabilities(const G4Fragment & aFragment)
{
  TotalEmissionProbability = 0.0;
  pcfvector::iterator aChannel; 
  for (aChannel=theChannels.begin(); aChannel != theChannels.end(); 
       aChannel++) 
    {
      // Calculate emission probailities
      // Compute total (integrated over kinetic energy) emission 
      // probability of a fragment and
      // Summing channel emission probabilities
      TotalEmissionProbability += (*aChannel)->CalcEmissionProbability(aFragment);
    }
  return TotalEmissionProbability;
}


G4VPreCompoundFragment * G4PreCompoundFragmentVector::
ChooseFragment(void)
{
  const G4int NumOfFrags = theChannels.size();
  G4std::vector<G4double> running;
  running.reserve(NumOfFrags);
  
  pcfvector::iterator i;
  G4double accumulation = 0.0;
  for (i = theChannels.begin(); i != theChannels.end(); ++i) {
    accumulation += (*i)->GetEmissionProbability();
    running.push_back(accumulation);
  }
	
  // Choose an emission channel
  G4double aChannel = G4UniformRand()*TotalEmissionProbability;
  G4int ChosenChannel = -1;
  G4std::vector<G4double>::iterator ich;
  for (ich = running.begin(); ich != running.end(); ++ich) 
    {
      if (aChannel <= *ich) 
	{
#ifdef G4NO_ISO_VECDIST
          G4std::vector<G4double>::difference_type n = 0;
          G4std::distance(running.begin(),ich,n);
          ChosenChannel = n;
#else
	  ChosenChannel = G4std::distance(running.begin(),ich);
#endif
	  break;
	}
    }
  running.clear();
  if (ChosenChannel < 0) 
  {
      G4cerr
	  << "G4PreCompoundFragmentVector::ChooseFragment: I can't determine a channel\n"
	  << "Probabilities: ";
      for (i = theChannels.begin(); i != theChannels.end(); ++i) 
      {
	  G4cout << (*i)->GetEmissionProbability() << "  ";
      }
      G4cout << '\n';
      return 0;
  }
  
  return theChannels[ChosenChannel];
}
