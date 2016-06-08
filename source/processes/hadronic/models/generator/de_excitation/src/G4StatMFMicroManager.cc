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
// $Id: G4StatMFMicroManager.cc,v 1.5 2001/10/05 16:13:44 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#include "G4StatMFMicroManager.hh"



// Copy constructor
G4StatMFMicroManager::G4StatMFMicroManager(const G4StatMFMicroManager & right)
{
    G4Exception("G4StatMFMicroManager::copy_constructor meant to not be accessable");
}

// Operators

G4StatMFMicroManager & G4StatMFMicroManager::
operator=(const G4StatMFMicroManager & right)
{
    G4Exception("G4StatMFMicroManager::operator= meant to not be accessable");
    return *this;
}


G4bool G4StatMFMicroManager::operator==(const G4StatMFMicroManager & right) const
{
    return false;
}
 

G4bool G4StatMFMicroManager::operator!=(const G4StatMFMicroManager & right) const
{
    return true;
}



// constructor
G4StatMFMicroManager::G4StatMFMicroManager(const G4Fragment & theFragment, const G4int multiplicity,
					   const G4double FreeIntE, const G4double SCompNuc) : 
    _Normalization(0.0)
{
    // Perform class initialization
    Initialize(theFragment,multiplicity,FreeIntE,SCompNuc);
}


// destructor
G4StatMFMicroManager::~G4StatMFMicroManager() 
{
    while (!_Partition.empty()) {
	delete _Partition.back();
	_Partition.pop_back();
    }
}



// Initialization method

void G4StatMFMicroManager::Initialize(const G4Fragment & theFragment, const G4int m, 
				      const G4double FreeIntE, const G4double SCompNuc) 
{
    G4int i;

    G4double U = theFragment.GetExcitationEnergy();

    G4double A = theFragment.GetA();
    G4double Z = theFragment.GetZ();
	
    // Statistical weights
    _WW = 0.0;

    // Mean breakup multiplicity
    _MeanMultiplicity = 0.0;

    // Mean channel temperature
    _MeanTemperature = 0.0;

    // Mean channel entropy
    _MeanEntropy = 0.0;	
	
    // Keep fragment atomic numbers
// 	G4int * FragmentAtomicNumbers = new G4int(G4int(A+0.5));
//	G4int * FragmentAtomicNumbers = new G4int(m);
    G4int FragmentAtomicNumbers[4];
	
    // We distribute A nucleons between m fragments mantaining the order
    // FragmentAtomicNumbers[m-1]>FragmentAtomicNumbers[m-2]>...>FragmentAtomicNumbers[0]
    // Our initial distribution is 
    // FragmentAtomicNumbers[m-1]=A, FragmentAtomicNumbers[m-2]=0, ..., FragmentAtomicNumbers[0]=0
    FragmentAtomicNumbers[m-1] = G4int(A);
    for (i = 0; i <  (m - 1); i++) FragmentAtomicNumbers[i] = 0;

    // We try to distribute A nucleons in partitions of m fragments
    // MakePartition return true if it is possible 
    // and false if it is not	
    while (MakePartition(m,FragmentAtomicNumbers)) {
	// Allowed partitions are stored and its probability calculated
			
	G4StatMFMicroPartition * aPartition = new G4StatMFMicroPartition(A,Z);
	G4double PartitionProbability = 0.0;
			
	for (i = m-1; i >= 0; i--) aPartition->SetPartitionFragment(FragmentAtomicNumbers[i]);
	PartitionProbability = aPartition->CalcPartitionProbability(U,FreeIntE,SCompNuc);
	_Partition.push_back(aPartition);
			
	_WW += PartitionProbability;
	_MeanMultiplicity += m*PartitionProbability;
	_MeanTemperature += aPartition->GetTemperature() * PartitionProbability;
	if (PartitionProbability > 0.0) 
	    _MeanEntropy += PartitionProbability * aPartition->GetEntropy();
			
    }
		
	
    // garbage collection
// 	delete [] FragmentAtomicNumbers;
	
}


G4bool G4StatMFMicroManager::MakePartition(const G4int k, G4int * ANumbers)
    // Distributes A nucleons between k fragments
    // mantaining the order ANumbers[k-1] > ANumbers[k-2] > ... > ANumbers[0]
    // If it is possible returns true. In other case returns false
{
    G4int l = 1;
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



void G4StatMFMicroManager::Normalize(const G4double Norm)
{
    _Normalization = Norm;
    _WW /= Norm;
    _MeanMultiplicity /= Norm;
    _MeanTemperature /= Norm;
    _MeanEntropy /= Norm; 
	
// 	for (G4int i = 0; i < _Partition.entries(); i++) _Partition(i)->Normalize(Norm);
	
    return;
}

G4StatMFChannel * G4StatMFMicroManager::ChooseChannel(const G4double A0, const G4double Z0, 
						      const G4double MeanT)
{
    G4double RandNumber = _Normalization * _WW * G4UniformRand();
    G4double AccumWeight = 0.0;
	
    for (unsigned int i = 0; i < _Partition.size(); i++) {
	AccumWeight += _Partition[i]->GetProbability();
	if (RandNumber < AccumWeight) 
	    return _Partition[i]->ChooseZ(A0,Z0,MeanT);
    }

    G4Exception
	("G4StatMFMicroCanonical::ChooseChannel: Couldn't find a channel.");
    return 0;
}
