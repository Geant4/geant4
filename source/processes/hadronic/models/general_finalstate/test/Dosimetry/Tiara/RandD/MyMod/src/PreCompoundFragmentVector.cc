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
// $Id: PreCompoundFragmentVector.cc,v 1.1 2003-10-08 12:32:13 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 

#include "PreCompoundFragmentVector.hh"

#include "PreCompoundNeutron.hh"
#include "PreCompoundProton.hh"
#include "PreCompoundDeuteron.hh"
#include "PreCompoundTriton.hh"
#include "PreCompoundHe3.hh"
#include "PreCompoundAlpha.hh"

PreCompoundFragmentVector::PreCompoundFragmentVector() :
  TotalEmissionProbability(0.0)
{
    // neutron 
    theChannels.push_back(new PreCompoundNeutron());
    // proton
    theChannels.push_back(new PreCompoundProton());
    // deuterium
    theChannels.push_back(new PreCompoundDeuteron());
    // triton
    theChannels.push_back(new PreCompoundTriton());
    // helium3
    theChannels.push_back(new PreCompoundHe3());
    // alpha
    theChannels.push_back(new PreCompoundAlpha());
}


PreCompoundFragmentVector::~PreCompoundFragmentVector()
{
  std::for_each(theChannels.begin(), theChannels.end(), 
		  DeleteFragment());
  theChannels.clear();
}

const PreCompoundFragmentVector & PreCompoundFragmentVector::operator=(const PreCompoundFragmentVector &right)
{
    G4Exception("G4PreCompoundFragmentVector::operator= meant to not be accessable");
    return *this;
}


G4bool PreCompoundFragmentVector::operator==(const PreCompoundFragmentVector &right) const
{
    return false;
}

G4bool PreCompoundFragmentVector::operator!=(const PreCompoundFragmentVector &right) const
{
    return true;
}



G4double PreCompoundFragmentVector::CalculateProbabilities(const G4Fragment & aFragment,G4double dLevelDensity)
{
  TotalEmissionProbability = 0.0;
  double dTmp;
  pcfvector::iterator aChannel; 
  for (aChannel=theChannels.begin(); aChannel != theChannels.end(); 
       aChannel++) 
    {
      // Calculate emission probailities
      // Compute total (integrated over kinetic energy) emission 
      // probability of a fragment and
      // Summing channel emission probabilities
      dTmp = (*aChannel)->CalcEmissionProbability(aFragment,dLevelDensity);
      TotalEmissionProbability += (dTmp<0) ? 0 : dTmp;
    }
  return TotalEmissionProbability;
}


VPreCompoundFragment * PreCompoundFragmentVector::ChooseFragment(void)
{ 
  pcfvector::iterator i;
  G4int n=0,ChosenChannel=-1;
  G4double accumulation = 0.0,dTmp;
  G4double aChannel = G4UniformRand()*TotalEmissionProbability;
  for (i = theChannels.begin(); i != theChannels.end(); ++i,n++) {
    accumulation += ((dTmp=(*i)->GetEmissionProbability())<0) ? 0 : dTmp;
    if(aChannel <= accumulation){
      ChosenChannel = n;
      break;
    }
  }
	
  // Choose an emission channel
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
