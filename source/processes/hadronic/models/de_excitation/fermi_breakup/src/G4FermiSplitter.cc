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
// by V. Lara

#include "G4FermiSplitter.hh"
#include "G4FermiIntegerPartition.hh"
#include "G4HadronicException.hh"
#include <sstream>


G4int G4FermiSplitter::Initialize(const G4int a, const G4int z, const G4int n)
{
  // Check argument correctness
  if (a < 0 || z < 0 || n < 0) 
    {
      std::ostringstream errOs;
      errOs << "G4FermiSplitter::Initialize() Error: Non valid arguments A = "
	    << a << " Z = " << z << " #fragments = " << n;
      throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    }
  if (z > a || n > a)
    {
      std::ostringstream errOs;
      errOs << "G4FermiSplitter::Initialize() Error: Non physical arguments = "
	    << a << " Z = " << z << " #fragments = " << n;
      throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    }
  A = a;
  Z = z;
  K = n;
  splits.clear();


  // Form all possible partition by combination
  // of A partitions and Z partitions (Z partitions include null parts)
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
              std::vector<std::vector<const G4VFermiFragment*> > tsplits;
              tsplits.clear();
              tsplits.resize(num_rows);
              G4int group_size = num_rows;
	      for (G4int i = 0; i < K; i++)
                {
		  az_pair = std::make_pair(partA[i],partZ[i]);
                  group_size /= multiplicities[i];
                  std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment*,
                    std::less<const std::pair<G4int,G4int> > >::iterator pos;
                  pos = theFragmentsPool->LowerBound(az_pair);
                  for (G4int k = 0; k < num_rows/group_size; k++)
                    {
                      if (pos == theFragmentsPool->UpperBound(az_pair))
                        {
                          pos = theFragmentsPool->LowerBound(az_pair);
                        }
                      for (G4int l = 0; l < group_size; l++)
                        {
                          tsplits[k*group_size+l].push_back(pos->second);
                        }
                      pos++;
                    }
                }
              // Remove wrong splits
              for (std::vector<std::vector<const G4VFermiFragment*> >::iterator 
                     itsplits1 = tsplits.begin(); itsplits1 != tsplits.end(); itsplits1++)
                {
                  std::sort((itsplits1)->begin(), (itsplits1)->end(),
                            std::greater<const G4VFermiFragment*>());
                }
              // add splits  (eliminating a few of them that are repeated)
              std::vector<std::vector<const G4VFermiFragment*> >::iterator 
                itlastsplit =  tsplits.begin();
              splits.push_back((*itlastsplit));
              for (std::vector<std::vector<const G4VFermiFragment*> >::iterator 
                     itsplits2 = itlastsplit+1; itsplits2 != tsplits.end(); itsplits2++)
                {
                  if ( (*itsplits2) != (*itlastsplit)) splits.push_back((*itsplits2));
                  itlastsplit++;
                }
	    }
	}
      while (PartitionZ.Next());
    } // partition of A
  while (PartitionA.Next());

  return splits.size();
}
