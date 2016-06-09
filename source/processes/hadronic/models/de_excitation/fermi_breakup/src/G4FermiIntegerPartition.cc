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
