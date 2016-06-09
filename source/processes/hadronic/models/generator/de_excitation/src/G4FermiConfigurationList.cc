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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4FermiConfigurationList.hh"




G4FermiConfigurationList::G4FermiConfigurationList():
  TotNumOfConfigurations(0)
{
  for (G4int i = 0; i < MaxNumOfFragments; i++) NumOfConfigurations[i] = 0;
}

G4FermiConfigurationList::G4FermiConfigurationList(const G4FermiConfigurationList &)
{
  G4Exception("G4FermiConfigurationList::copy_constructor meant to not be accessable");
}


const G4FermiConfigurationList & G4FermiConfigurationList::operator=(const G4FermiConfigurationList &)
{
  G4Exception("G4FermiConfigurationList::operator= meant to not be accessable");
  return *this;
}


G4bool G4FermiConfigurationList::operator==(const G4FermiConfigurationList &) const
{
  return false;
}

G4bool G4FermiConfigurationList::operator!=(const G4FermiConfigurationList &) const
{
  return true;
}



G4bool G4FermiConfigurationList::Initialize(const G4int A, const G4int Z, const G4double TotalEnergyRF)
{
  //
  // let's split nucleus into k = 2,...,A fragments
  //
  Configurations.clear();
  NormalizedWeights.clear();
  G4FermiConfiguration aConfiguration;
  std::vector<G4double> NOTNormalizedWeights;
  G4double NormStatWeight = 0.0;
  for (G4int k = 2; k <= A; k++) 
    {
      // Initialize Configuration for k fragments
      aConfiguration.Initialize(k);
      G4bool SplitSuccesed;
      do 
	{
	  // Splits the nucleus into k fragments
	  SplitSuccesed = aConfiguration.SplitNucleus(A,Z);
	  if (SplitSuccesed) 
	    {
	      TotNumOfConfigurations++;
	      NumOfConfigurations[k-1]++;
	    
	      // Non-Normalized statistical weight for given channel with k fragments
	      G4double StatWeight = aConfiguration.DecayProbability(A,TotalEnergyRF);
	      NormStatWeight += StatWeight;
	      // Statistical weights (it will be normalized...)
	      NOTNormalizedWeights.push_back(StatWeight);	
	    
	      // Store configuration
	      Configurations.push_back(aConfiguration);
	    }
	  // Repeat splitting into k fragments (it may be several posibilities for a choosen K)
	} 
      while (SplitSuccesed);
    }

  if (NormStatWeight > 0.0) 
    {
      // Let's normalize statistical weights of channels
      for (G4int i = 0; i < TotNumOfConfigurations; i++) 
	{
	  NormalizedWeights.push_back(NOTNormalizedWeights[i]/NormStatWeight);
	}
      return true;
    }
  else return false;

}



G4FermiConfiguration G4FermiConfigurationList::ChooseConfiguration(void)
{
  G4double RandomWeight =  G4UniformRand();
  G4double AcumWeight = 0.0;
  G4int thisConfig = 0;
  do 
    {
      AcumWeight += NormalizedWeights[thisConfig];  // We are adding the prob. of each configuration
      thisConfig++;
    } 
  while ((thisConfig <= TotNumOfConfigurations) && (AcumWeight < RandomWeight));

  return Configurations[thisConfig - 1];

}
