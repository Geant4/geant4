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
// $Id: G4PhysListFactoryAlt.cc 83151 2014-08-06 16:22:30Z klg $
//
//---------------------------------------------------------------------------
//
// ClassName:  g4alt::G4PhysListFactory
//
// Author:  R. Hatcher  2014-10-15
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4PhysListFactoryAlt.hh"
#include "G4PhysListRegistry.hh"

namespace g4alt {

G4PhysListFactory::G4PhysListFactory() 
  : defName("FTFP_BERT"), verbose(1)
{
}

G4PhysListFactory::~G4PhysListFactory()
{}

G4VModularPhysicsList* 
G4PhysListFactory::ReferencePhysList()
{
  // instantiate PhysList by environment variable "PHYSLIST"
  G4String name = "";
  char* path = getenv("PHYSLIST");
  if (path) {
    name = G4String(path);
  } else {
    name = defName;
    G4cout << "### g4alt::G4PhysListFactory WARNING: "
	   << " environment variable PHYSLIST is not defined"
	   << G4endl
	   << "    Default Physics Lists " << name 
	   << " is instantiated" 
	   << G4endl;
  }
  return GetReferencePhysList(name);
}

G4VModularPhysicsList* 
G4PhysListFactory::GetReferencePhysList(const G4String& name)
{
  G4PhysListRegistry* registry = G4PhysListRegistry::Instance();
  G4int oldRegVerbose = registry->GetVerbose();
  registry->SetVerbose(verbose);
  G4VModularPhysicsList* pl = registry->GetModularPhysicsList(name);
  registry->SetVerbose(oldRegVerbose);
  return pl;
}
  
G4bool G4PhysListFactory::IsReferencePhysList(const G4String& name)
{
  return G4PhysListRegistry::Instance()->IsReferencePhysList(name);
}

const std::vector<G4String>& 
G4PhysListFactory::AvailablePhysLists() const
{
  return G4PhysListRegistry::Instance()->AvailablePhysLists();
}

const std::vector<G4String>& 
G4PhysListFactory::AvailablePhysListsEM() const
{
  return G4PhysListRegistry::Instance()->AvailablePhysListsEM();
}

 void G4PhysListFactory::PrintAvailablePhysLists() const
{
  G4PhysListRegistry::Instance()->PrintAvailablePhysLists();
}

} // end-of-space 'g4alt'
