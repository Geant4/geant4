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
