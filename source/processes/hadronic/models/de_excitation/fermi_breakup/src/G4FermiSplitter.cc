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
// by V. Lara

#include "G4FermiSplitter.hh"
#include "G4FermiIntegerPartition.hh"
#include "G4HadronicException.hh"
#include <strstream>


G4int G4FermiSplitter::Initialize(const G4int a, const G4int z, const G4int n)
{
  // Check argument correctness
  if (a < 0 || z < 0 || n < 0) 
    {
      std::ostrstream errOs;
      errOs << "G4FermiSplitter::Initialize() Error: Non valid arguments A = "
	    << a << " Z = " << z << " #fragments = " << n;
      throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    }
  if (z > a || n > a)
    {
      std::ostrstream errOs;
      errOs << "G4FermiSplitter::Initialize() Error: Non physical arguments = "
	    << a << " Z = " << z << " #fragments = " << n;
      throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    }
  A = a;
  Z = z;
  K = n;

  G4FermiIntegerPartition PartitionA;
  PartitionA.Initialize(A,K);
  do // for each partition of A
    {
      G4FermiIntegerPartition PartitionZ;
      PartitionZ.EnableNull();
      PartitionZ.Initialize(Z,K);
      std::vector<G4int> partA = PartitionA.GetPartition();
      std::vector<G4int> multiplicities;
      do  // for each partition of Z
	{
	  std::vector<G4int> partZ = PartitionZ.GetPartition();
	  // Get multiplicities
	  // provided that for each a,z pair there could be several 
	  // posibilities, we have to count them in order to create 
	  // all possible combinations.
	  multiplicities.clear();
	  G4int num_rows = 1;
	  std::pair<G4int,G4int> az_pair;
	  for (G4int i = 0; i < K; i++)
	    {
	      az_pair = std::make_pair(partA[i],partZ[i]);
	      G4int m = theFragmentsPool->Count(az_pair);
	      if (m > 0)
		{
		  multiplicities.push_back(m);
		  num_rows *= m;
		}
	      else
		{
		  break;
		}
	    }
	  // if size of multiplicities is equal to K, it
	  // means that all combinations a,z exist in the 
	  // fragments pool
	  if (static_cast<G4int>(multiplicities.size()) == K)
	    {
	      splits.clear();
	      splits.resize(num_rows);
	      for (G4int i = 0; i < K; i++)
		{
		  az_pair = std::make_pair(partA[i],partZ[i]);
		  for (G4int j = 0; j < num_rows/multiplicities[i]; j++)
		    { // number of times that we have to introduce the same data
		      G4int k=0;
		      std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment*,
			std::less<const std::pair<G4int,G4int> > >::iterator pos;
		      for (pos = theFragmentsPool->LowerBound(az_pair);
			   pos != theFragmentsPool->UpperBound(az_pair); ++pos)
			{
			  G4int tmp = j*multiplicities[i]+k;
			  splits[tmp].push_back(pos->second);
			  k++;
			}
		    }
		}
	    }
	}
      while (PartitionZ.Next());
    } // partition of A
  while (PartitionA.Next());

  return splits.size();
}
