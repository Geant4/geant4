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
// $Id: G4PreCompoundFragmentVector.cc,v 1.7.2.1 2009/03/03 13:17:04 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-01 $
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
