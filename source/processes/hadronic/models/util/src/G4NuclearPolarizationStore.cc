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
// $Id: G4NuclearPolarizationStore.cc 97302 2016-06-01 09:30:11Z gcosmo $
//
// 23-Jan-2009 V.Ivanchenko make the class to be a singleton
// 17-Aug-2012 V.Ivanchenko added hadronic model factories

#include "G4NuclearPolarizationStore.hh"
#include "G4SystemOfUnits.hh"

G4ThreadLocal G4NuclearPolarizationStore* 
G4NuclearPolarizationStore::instance = nullptr;

G4NuclearPolarizationStore* G4NuclearPolarizationStore::GetInstance()
{
  if(nullptr == instance) {
    static G4ThreadLocalSingleton<G4NuclearPolarizationStore> inst;
    instance = inst.Instance();
  }
  return instance;
}

G4NuclearPolarizationStore::G4NuclearPolarizationStore()
{
  for(G4int i=0; i<maxNumStates; ++i) { nuclist[i] = nullptr; }
  oldIdx = 0;
}

G4NuclearPolarizationStore::~G4NuclearPolarizationStore()
{
  for(G4int i=0; i<maxNumStates; ++i) { 
    delete nuclist[i];
    nuclist[i] = nullptr; 
  }
}

void G4NuclearPolarizationStore::Register(G4NuclearPolarization* ptr)
{
  G4int idx = -1;
  for(G4int i=0; i<maxNumStates; ++i) { 
    if(ptr == nuclist[i]) { return; }
    if(nullptr == nuclist[i]) { idx = i; }
  }
  if(idx >= 0) {
    nuclist[idx] = ptr;
    return;
  }
  // delete oldest object
  delete nuclist[oldIdx];
  nuclist[oldIdx] = ptr;
  // redefine oldIdx
  ++oldIdx;
  if(oldIdx >= maxNumStates) { oldIdx = 0; }
}

G4NuclearPolarization* 
G4NuclearPolarizationStore::FindOrBuild(G4int Z, G4int A, G4double Eexc)
{
  static const G4double tolerance = 10.*CLHEP::eV;
  for(G4int i=0; i<maxNumStates; ++i) { 
    auto nucp = nuclist[i];
    if(nucp && Z == nucp->GetZ() && A == nucp->GetA() && 
       std::abs(Eexc - nucp->GetExcitationEnergy()) < tolerance) { 
      return nucp; 
    }
  }
  G4NuclearPolarization* ptr = new G4NuclearPolarization(Z, A, Eexc);
  Register(ptr);
  return ptr;   
}

void G4NuclearPolarizationStore::RemoveMe(G4NuclearPolarization* ptr)
{
  for(G4int i=0; i<maxNumStates; ++i) { 
    if(ptr == nuclist[i]) { 
      delete ptr;
      nuclist[i] = nullptr;
      // do we need redefine oldIdx?
      if(i == oldIdx) {
	for(G4int j=0; j<maxNumStates; ++j) { 
	  if(j != i && nullptr != nuclist[j]) { 
	    oldIdx = j;
            return;
	  }
	}
	oldIdx = i;   
      }
      return;
    }
  }
}
