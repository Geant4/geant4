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
// $Id: G4FermiFragmentsPool.hh,v 1.5 2006-06-29 20:13:13 gunter Exp $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko - more clean usage of static


#ifndef G4FermiFragmentsPool_hh 
#define G4FermiFragmentsPool_hh 1

#include "globals.hh"
#include "G4VFermiFragment.hh"
#include "G4FermiConfiguration.hh"
#include "G4FermiPhaseSpaceDecay.hh"
#include <vector>

class G4FermiFragmentsPool
{
public:

  static G4FermiFragmentsPool* Instance();

  ~G4FermiFragmentsPool();

  G4bool IsApplicable(G4int Z, G4int A, G4double mass) const;
 
  const std::vector<const G4FermiConfiguration*>* 
  GetConfigurationList(G4int Z, G4int A, G4double mass) const;

  const G4VFermiFragment* GetFragment(G4int Z, G4int A) const;

  const G4FermiPhaseSpaceDecay* GetFermiPhaseSpaceDecay() const;

  G4int GetMaxZ() const;

  G4int GetMaxA() const;

  void DumpFragment(const G4VFermiFragment* f) const;

  void Dump() const;
  
private:

  G4FermiFragmentsPool();

  void Initialise();

  G4bool IsExist(G4int Z, G4int A, std::vector<const G4VFermiFragment*>&) const;

  G4bool IsAvailable(G4int Z, G4int A) const;

  static G4FermiFragmentsPool* theInstance;

  std::vector<const G4VFermiFragment*> fragment_pool;

  G4int maxZ;
  G4int maxA;
  G4int verbose;
 
  // list of configuration sorted by A for 1, 2, 3, 4 final fragments
  std::vector<const G4FermiConfiguration*> list1[17]; 
  std::vector<const G4FermiConfiguration*> list2[17]; 
  std::vector<const G4FermiConfiguration*> list3[17];
  std::vector<const G4FermiConfiguration*> list4[17];

  G4FermiPhaseSpaceDecay thePhaseSpace;
};

#endif

