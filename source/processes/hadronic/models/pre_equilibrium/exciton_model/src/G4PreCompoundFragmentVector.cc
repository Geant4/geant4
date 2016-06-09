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
// $Id: G4PreCompoundFragmentVector.cc,v 1.3 2003/11/04 11:36:27 lara Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 

#include "G4PreCompoundFragmentVector.hh"
#include "G4HadronicException.hh"

const G4PreCompoundFragmentVector & 
G4PreCompoundFragmentVector::
operator=(const G4PreCompoundFragmentVector &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4PreCompoundFragmentVector::operator= meant to not be accessable");
    return *this;
}


G4bool G4PreCompoundFragmentVector::
operator==(const G4PreCompoundFragmentVector &) const
{
    return false;
}

G4bool G4PreCompoundFragmentVector::
operator!=(const G4PreCompoundFragmentVector &) const
{
    return true;
}



G4double G4PreCompoundFragmentVector::
CalculateProbabilities(const G4Fragment & aFragment)
{
  TotalEmissionProbability = 0.0;
  pcfvector::iterator aChannel; 
  for (aChannel=theChannels->begin(); aChannel != theChannels->end(); 
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
  const G4int NumOfFrags = theChannels->size();
  std::vector<G4double> running;
  running.reserve(NumOfFrags);
  
  pcfvector::iterator i;
  G4double accumulation = 0.0;
  for (i = theChannels->begin(); i != theChannels->end(); ++i) {
    accumulation += (*i)->GetEmissionProbability();
    running.push_back(accumulation);
  }
	
  // Choose an emission channel
  G4double aChannel = G4UniformRand()*TotalEmissionProbability;
  G4int ChosenChannel = -1;
  std::vector<G4double>::iterator ich;
  for (ich = running.begin(); ich != running.end(); ++ich) 
    {
      if (aChannel <= *ich) 
	{
#ifdef G4NO_ISO_VECDIST
          std::vector<G4double>::difference_type n = 0;
          std::distance(running.begin(),ich,n);
          ChosenChannel = n;
#else
	  ChosenChannel = std::distance(running.begin(),ich);
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
      for (i = theChannels->begin(); i != theChannels->end(); ++i) 
	{
	  G4cout << (*i)->GetEmissionProbability() << "  ";
	}
      G4cout << '\n';
      return 0;
    }
  else
    {
      for (i = theChannels->begin(); i != theChannels->end(); ++i) 
	{
	  (*i)->IncrementStage();
	}
    } 

  return theChannels->operator[](ChosenChannel);
}
