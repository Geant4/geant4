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
