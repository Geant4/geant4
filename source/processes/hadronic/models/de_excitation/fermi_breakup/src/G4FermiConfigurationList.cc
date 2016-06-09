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
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4FermiConfigurationList.hh"
#include "G4FermiSplitter.hh"
#include "G4HadronicException.hh"

G4FermiFragmentsPool & G4FermiConfigurationList::GetFragmentsPoolInstance()
{
  static G4FermiFragmentsPool theFragmentsPool;
  return theFragmentsPool;
}

G4FermiConfigurationList::G4FermiConfigurationList()
{
}

G4FermiConfigurationList::G4FermiConfigurationList(const G4FermiConfigurationList &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4FermiConfigurationList::copy_constructor meant to not be accessable");
}


const G4FermiConfigurationList & G4FermiConfigurationList::
operator=(const G4FermiConfigurationList &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4FermiConfigurationList::operator= meant to not be accessable");
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



G4bool G4FermiConfigurationList::
Initialize(const G4int A, const G4int Z, const G4double TotalEnergyRF)
{
  //
  // let's split nucleus into k = 2,...,A fragments
  //
  Configurations.clear();
  NormalizedWeights.clear();
  G4FermiSplitter aSplitter(&GetFragmentsPoolInstance());
  G4double NormStatWeight = 0.0;
  for (G4int k = 2; k <= A; k++) 
    {
      // Initialize Configuration for k fragments
      aSplitter.Initialize(A,Z,k);
      for (G4int i = 0; i < aSplitter.GetNumberOfSplits(); i++)
	{
	  // Create a configuration from a split
	  G4FermiConfiguration * aConfiguration = new G4FermiConfiguration(aSplitter.GetSplit(i));

	  // Store configuration
	  Configurations.push_back(aConfiguration);
	  	    
	  // Non-Normalized statistical weight for given channel with k fragments
	  G4double StatWeight = aConfiguration->DecayProbability(A,TotalEnergyRF);
	  NormStatWeight += StatWeight;
	  // Statistical weights (it will be normalized...)
	  NormalizedWeights.push_back(StatWeight);	
	} 
    }
  
  if (NormStatWeight > 0.0) 
    {
      // Let's normalize statistical weights of channels
      std::transform(NormalizedWeights.begin(), NormalizedWeights.end(), 
		     NormalizedWeights.begin(),
		     std::bind2nd(std::divides<G4double>(),NormStatWeight));
      return true;
    }

  return false;
}



G4FermiConfiguration G4FermiConfigurationList::ChooseConfiguration(void)
{
  G4double RandomWeight =  G4UniformRand();
  G4double AcumWeight = 0.0;
  std::vector<G4double>::iterator thisConfig;
  for (thisConfig = NormalizedWeights.begin(); thisConfig != NormalizedWeights.end(); ++thisConfig)
    {
      // We are adding the prob. of each configuration
      AcumWeight += *thisConfig; 
      if (AcumWeight >= RandomWeight) break;
    } 
  if (thisConfig == NormalizedWeights.end())
    {
      throw G4HadronicException(__FILE__, __LINE__, "G4FermiConfigurationList::ChooseConfigration: Cannot choose a configuration");
      G4FermiConfiguration dummy;
      return dummy;
    }
  else 
    {
#ifdef G4NO_ISO_VECDIST
      std::vector<G4double>::difference_type n = 0;
      std::distance(NormalizedWeights.begin(),thisConfig,n);
      return *(Configurations[n]);
#else
      return *(Configurations[std::distance(NormalizedWeights.begin(),thisConfig)]);
#endif
    }
}
