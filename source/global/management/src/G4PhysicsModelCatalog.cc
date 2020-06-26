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
// G4PhysicsModelCatalog class implementation
//
// Author: M.Asai (SLAC), 26 September 2013
// --------------------------------------------------------------------

#include "G4PhysicsModelCatalog.hh"

#include "G4Threading.hh"

std::vector<G4String>* G4PhysicsModelCatalog::theCatalog = nullptr;

// --------------------------------------------------------------------
G4PhysicsModelCatalog::G4PhysicsModelCatalog()
{
  if(theCatalog == nullptr)
  {
    static std::vector<G4String> catal;
    theCatalog = &catal;
  }
}

// --------------------------------------------------------------------
G4PhysicsModelCatalog::~G4PhysicsModelCatalog() {}

// --------------------------------------------------------------------
G4int G4PhysicsModelCatalog::Register(const G4String& name)
{
  G4PhysicsModelCatalog();
  G4int idx = GetIndex(name);
  if(idx >= 0)
    return idx;
#ifdef G4MULTITHREADED
  if(G4Threading::IsWorkerThread())
    return -1;
#endif
  theCatalog->push_back(name);
  return theCatalog->size() - 1;
}

// --------------------------------------------------------------------
const G4String& G4PhysicsModelCatalog::GetModelName(G4int idx)
{
  static const G4String undef = "Undefined";
  if(idx >= 0 && idx < Entries())
    return (*theCatalog)[idx];
  return undef;
}

// --------------------------------------------------------------------
G4int G4PhysicsModelCatalog::GetIndex(const G4String& name)
{
  for(G4int idx = 0; idx < Entries(); ++idx)
  {
    if((*theCatalog)[idx] == name)
      return idx;
  }
  return -1;
}

// --------------------------------------------------------------------
G4int G4PhysicsModelCatalog::Entries()
{
  return (theCatalog != nullptr) ? G4int(theCatalog->size()) : -1;
}

// --------------------------------------------------------------------
void G4PhysicsModelCatalog::Destroy() {}
