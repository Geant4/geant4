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


#ifndef G4FermiSplitter_hh
#define G4FermiSplitter_hh


#include "G4VFermiFragment.hh"
#include "G4FermiFragmentsPool.hh"
#include <vector>

class G4FermiSplitter
{
public:
  inline G4FermiSplitter(G4FermiFragmentsPool*);
  inline ~G4FermiSplitter();

  inline G4FermiSplitter(const G4FermiSplitter&);

  inline const G4FermiSplitter& operator=(const G4FermiSplitter&);
  inline G4bool operator==(const G4FermiSplitter&);
  inline G4bool operator!=(const G4FermiSplitter&);

  G4int Initialize(const G4int a, const G4int z, const G4int n);

  inline G4int GetNumberOfSplits() const;
  inline std::vector<const G4VFermiFragment*> GetSplit(const G4int i);

private:
  inline G4FermiSplitter();

private:

  G4FermiFragmentsPool * theFragmentsPool;
  G4int A;
  G4int Z;
  G4int K;
  std::vector<std::vector<const G4VFermiFragment*> > splits;
};

#include "G4FermiSplitter.icc"

#endif
