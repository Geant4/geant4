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


#include "G4FermiIntegerPartition.hh"

void G4FermiIntegerPartition::Initialize(const G4int N, const G4int size)
{
  total = N;
  int S(size);
  if (!enableNull)
    {
      S = std::min(size,N);
    }
  partition.clear();
  partition.reserve(S);
  G4int ini(1);
  if (enableNull)
    {
      ini = 0;
      partition.push_back(total);
    }
  else
    {
      partition.push_back(total-S+1);
    }
  for (G4int i = 1; i < S; i++) partition.push_back(ini);

#ifdef G4FermiIntegerPartition_debug
  this->TestPartition();
#endif

  return;
}


G4bool G4FermiIntegerPartition::Next()
{
  std::vector<G4int>::iterator first = partition.begin();
  std::vector<G4int>::iterator i = first+1;

  while ( i != partition.end() && (*first) <= (*i)+1 ) i++; 

  if ( i == partition.end() )
    {
      return false;
    }
  else
    {
      (*i)++;
      G4int d = -1;
      for (std::vector<G4int>::iterator j = i-1; j != first; j--)
	{
	  d += (*j) - (*i);
	  (*j) = (*i);
	}
      (*first) += d;
    }

#ifdef G4FermiIntegerPartition_debug
  this->TestPartition();
#endif

  return true;
}

#ifdef G4FermiIntegerPartition_debug

void G4FermiIntegerPartition::TestPartition()
{
  G4int tmp = this->GetSum();
  if (total != tmp)
    {
      std::cerr << "G4FermiIntegerPartition Error: " 
		  << "Partition of " << total << " into " << partition.size()
		  << " is [";
      for (std::vector<G4int>::iterator it = partition.begin();
	   it != partition.end(); ++it)
	{
	  std::cerr << *it << ',';
	}
      std::cerr << "\b\b] = " << tmp << '\n';
    }
  return;
}

#endif
