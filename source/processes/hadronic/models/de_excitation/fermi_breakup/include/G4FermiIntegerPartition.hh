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


#ifndef G4FermiIntegerPartition_hh
#define G4FermiIntegerPartition_hh

#define G4FermiIntegerPartition_debug 1
#include "globals.hh"
#include <vector>

#ifdef G4FermiIntegerPartition_debug
#include <numeric>
#endif


class G4FermiIntegerPartition
{
public:
  inline G4FermiIntegerPartition();
  inline ~G4FermiIntegerPartition();
  inline G4FermiIntegerPartition(const G4FermiIntegerPartition&);

  inline const G4FermiIntegerPartition & operator=(const G4FermiIntegerPartition&);
  G4bool operator==(const G4FermiIntegerPartition&);
  G4bool operator!=(const G4FermiIntegerPartition&);

  inline void EnableNull(const G4bool v=true);
  void Initialize(const G4int, const G4int);
  G4bool Next();
  inline std::vector<G4int> GetPartition() const;

private:
#ifdef G4FermiIntegerPartition_debug
  inline G4int GetSum();
  void TestPartition();
#endif

private:
  G4int total;
  G4bool enableNull;
  std::vector<G4int> partition;
};

#include "G4FermiIntegerPartition.icc"

#endif
