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
// $Id: G4StatMFMicroManager.cc 91834 2015-08-07 07:24:22Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4StatMFMicroManager.hh"
#include "G4HadronicException.hh"

// Copy constructor
G4StatMFMicroManager::G4StatMFMicroManager(const G4StatMFMicroManager & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMicroManager::copy_constructor meant to not be accessable");
}

// Operators

G4StatMFMicroManager & G4StatMFMicroManager::
operator=(const G4StatMFMicroManager & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StatMFMicroManager::operator= meant to not be accessable");
    return *this;
}


G4bool G4StatMFMicroManager::operator==(const G4StatMFMicroManager & ) const
{
    return false;
}
 

G4bool G4StatMFMicroManager::operator!=(const G4StatMFMicroManager & ) const
{
    return true;
}

// constructor
G4StatMFMicroManager::G4StatMFMicroManager(const G4Fragment & theFragment, 
					   G4int multiplicity,
					   G4double FreeIntE, G4double SCompNuc) : 
  _Normalization(0.0)
{
  // Perform class initialization
  Initialize(theFragment,multiplicity,FreeIntE,SCompNuc);
}

// destructor
G4StatMFMicroManager::~G4StatMFMicroManager() 
{
  if (!_Partition.empty()) 
    {
      std::for_each(_Partition.begin(),_Partition.end(),
		      DeleteFragment());
    }
}

void G4StatMFMicroManager::Initialize(const G4Fragment & theFragment, G4int im, 
				      G4double FreeIntE, G4double SCompNuc) 
{
  G4int i;

  G4double U = theFragment.GetExcitationEnergy();

  G4int A = theFragment.GetA_asInt();
  G4int Z = theFragment.GetZ_asInt();
	
  // Statistical weights
  _WW = 0.0;

  // Mean breakup multiplicity
  _MeanMultiplicity = 0.0;

  // Mean channel temperature
  _MeanTemperature = 0.0;

  // Mean channel entropy
  _MeanEntropy = 0.0;	
	
  // Keep fragment atomic numbers
  // 	G4int * FragmentAtomicNumbers = new G4int(static_cast<G4int>(A+0.5));
  //	G4int * FragmentAtomicNumbers = new G4int(m);
  G4int FragmentAtomicNumbers[4];
	
  // We distribute A nucleons between m fragments mantaining the order
  // FragmentAtomicNumbers[m-1]>FragmentAtomicNumbers[m-2]>...>FragmentAtomicNumbers[0]
  // Our initial distribution is 
  // FragmentAtomicNumbers[m-1]=A, FragmentAtomicNumbers[m-2]=0, ..., FragmentAtomicNumbers[0]=0
  FragmentAtomicNumbers[im-1] = A;
  for (i = 0; i <  (im - 1); i++) FragmentAtomicNumbers[i] = 0;

  // We try to distribute A nucleons in partitions of m fragments
  // MakePartition return true if it is possible 
  // and false if it is not	

  // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
  while (MakePartition(im,FragmentAtomicNumbers)) {
    // Allowed partitions are stored and its probability calculated
			
    G4StatMFMicroPartition * aPartition = new G4StatMFMicroPartition(A,Z);
    G4double PartitionProbability = 0.0;
			
    for (i = im-1; i >= 0; i--) aPartition->SetPartitionFragment(FragmentAtomicNumbers[i]);
    PartitionProbability = aPartition->CalcPartitionProbability(U,FreeIntE,SCompNuc);
    _Partition.push_back(aPartition);
			
    _WW += PartitionProbability;
    _MeanMultiplicity += im*PartitionProbability;
    _MeanTemperature += aPartition->GetTemperature() * PartitionProbability;
    if (PartitionProbability > 0.0) 
      _MeanEntropy += PartitionProbability * aPartition->GetEntropy();
  }
}

G4bool G4StatMFMicroManager::MakePartition(G4int k, G4int * ANumbers)
// Distributes A nucleons between k fragments
// mantaining the order ANumbers[k-1] > ANumbers[k-2] > ... > ANumbers[0]
// If it is possible returns true. In other case returns false
{
  G4int l = 1;
  // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
  while (l < k) {
    G4int tmp = ANumbers[l-1] + ANumbers[k-1];
    ANumbers[l-1] += 1;
    ANumbers[k-1] -= 1;
    if (ANumbers[l-1] > ANumbers[l] || ANumbers[k-2] > ANumbers[k-1]) {
      ANumbers[l-1] = 1;
      ANumbers[k-1] = tmp - 1;
      l++;
    } else return true;
  }
  return false;
}

void G4StatMFMicroManager::Normalize(G4double Norm)
{
  _Normalization = Norm;
  _WW /= Norm;
  _MeanMultiplicity /= Norm;
  _MeanTemperature /= Norm;
  _MeanEntropy /= Norm; 
	
  return;
}

G4StatMFChannel* 
G4StatMFMicroManager::ChooseChannel(G4int A0, G4int Z0, G4double MeanT)
{
  G4double RandNumber = _Normalization * _WW * G4UniformRand();
  G4double AccumWeight = 0.0;
	
  for (std::vector<G4StatMFMicroPartition*>::iterator i = _Partition.begin();
       i != _Partition.end(); ++i)
    {
	AccumWeight += (*i)->GetProbability();
	if (RandNumber < AccumWeight)
	  return (*i)->ChooseZ(A0,Z0,MeanT);
    }

  throw G4HadronicException(__FILE__, __LINE__, 
			    "G4StatMFMicroCanonical::ChooseChannel: Couldn't find a channel.");
  return 0;
}
