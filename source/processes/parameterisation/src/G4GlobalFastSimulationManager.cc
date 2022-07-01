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
//
//  
//---------------------------------------------------------------
//
//  G4GlobalFastSimulationManager.cc
//
//  Description:
//    A singleton class which manages the Fast Simulation managers 
//    attached to envelopes. Implementation.
//
//  History:
//    June 98: Verderi && MoraDeFreitas - "G4ParallelWorld" becomes
//             "G4FlavoredParallelWorld"; some method name changes;
//             GetFlavoredWorldForThis now returns a 
//             G4FlavoredParallelWorld pointer.
//    Feb 98: Verderi && MoraDeFreitas - First Implementation.
//    March 98: correction to instanciate dynamically the manager
//    May 07: Move to parallel world scheme
//
//---------------------------------------------------------------

#include "G4GlobalFastSimulationManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4TransportationManager.hh"
#include "G4FastSimulationMessenger.hh"
#include "G4RegionStore.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicalVolumeStore.hh"


// ------------------------------------------
// -- static instance pointer initialisation:
// ------------------------------------------
G4ThreadLocal G4GlobalFastSimulationManager* G4GlobalFastSimulationManager::fGlobalFastSimulationManager = 0;

// --------------------------------------------------
// -- static methods to retrieve the manager pointer:
// --------------------------------------------------
G4GlobalFastSimulationManager* G4GlobalFastSimulationManager::GetGlobalFastSimulationManager()
{
  if(!fGlobalFastSimulationManager)
    fGlobalFastSimulationManager = new G4GlobalFastSimulationManager;
  
  return fGlobalFastSimulationManager;
}


G4GlobalFastSimulationManager* G4GlobalFastSimulationManager::GetInstance()
{
  return G4GlobalFastSimulationManager::GetGlobalFastSimulationManager();
}

// ---------------
// -- constructor
// ---------------
G4GlobalFastSimulationManager::G4GlobalFastSimulationManager()
{
  fTheFastSimulationMessenger = new G4FastSimulationMessenger(this);
}

// -------------
// -- destructor
// -------------
G4GlobalFastSimulationManager::~G4GlobalFastSimulationManager()
{
  delete fTheFastSimulationMessenger;
  fTheFastSimulationMessenger = 0;
}

// ----------------------
// -- management methods:
// ----------------------
void G4GlobalFastSimulationManager::
AddFastSimulationManager(G4FastSimulationManager* fsmanager)
{
  ManagedManagers.push_back(fsmanager);
}

void G4GlobalFastSimulationManager::
RemoveFastSimulationManager(G4FastSimulationManager* fsmanager)
{
  ManagedManagers.remove(fsmanager);
}

void G4GlobalFastSimulationManager::AddFSMP(G4FastSimulationManagerProcess* fp)
{
  fFSMPVector.push_back(fp);
}

void G4GlobalFastSimulationManager::RemoveFSMP(G4FastSimulationManagerProcess* fp)
{
  fFSMPVector.remove(fp);
}

void G4GlobalFastSimulationManager::ActivateFastSimulationModel(const G4String& aName)
{
  G4bool result = false;
  for (size_t ifsm=0; ifsm<ManagedManagers.size(); ifsm++)
    result = result || ManagedManagers[ifsm]->
                       ActivateFastSimulationModel(aName);
  if(result) 
    G4cout << "Model " << aName << " activated.";
  else
    G4cout << "Model " << aName << " not found.";
  G4cout << G4endl;
}

void G4GlobalFastSimulationManager::InActivateFastSimulationModel(const G4String& aName)
{
  G4bool result = false;
  for (size_t ifsm=0; ifsm<ManagedManagers.size(); ifsm++)
    result = result || ManagedManagers[ifsm]->
                       InActivateFastSimulationModel(aName);
  if (result) G4cout << "Model " << aName << " inactivated.";
  else        G4cout << "Model " << aName << " not found.";
  G4cout << G4endl;
}

void G4GlobalFastSimulationManager::Flush()
{
  // loop over all models (that need flushing?) and flush
  for (size_t ifsm=0; ifsm<ManagedManagers.size(); ifsm++)
    ManagedManagers[ifsm]->FlushModels();
}

// ---------------------------------
// -- display fast simulation setup:
// ---------------------------------
void G4GlobalFastSimulationManager::ShowSetup()
{
  std::vector<G4VPhysicalVolume*> worldDone;
  G4VPhysicalVolume* world;
  G4RegionStore* regions = G4RegionStore::GetInstance();
  // ----------------------------------------------------
  // -- loop on regions to get the list of world volumes:
  // ----------------------------------------------------
  G4cout << "\nFast simulation setup:" << G4endl;
  for (size_t i=0; i<regions->size(); i++)
    {
      world = (*regions)[i]->GetWorldPhysical();
      if (world == nullptr) // region does not belong to any (existing) world
        {
          continue;
        }
      G4bool newWorld = true;
      for (size_t ii=0; ii<worldDone.size(); ii++) if (worldDone[ii] == world) {newWorld = false; break;}
      if (newWorld)
	{
	  worldDone.push_back(world);
	  G4Region* worldRegion = world->GetLogicalVolume()->GetRegion();
	  // -- preambule: print physical volume and region names...
	  if (world == G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume())
	    G4cout << "\n    * Mass Geometry with ";
	  else
	    G4cout << "\n    * Parallel Geometry with ";
	  G4cout << "world volume: `" << world->GetName() << "' [region : `" << worldRegion->GetName() << "']"  << G4endl;
	  // -- ... and print G4FSMP(s) attached to this world volume:
	  G4bool findG4FSMP(false);
	  // -- show to what particles this G4FSMP is attached to:
	  std::vector<G4ParticleDefinition*> particlesKnown;
	  for (size_t ip=0; ip<fFSMPVector.size(); ip++)
	    if (fFSMPVector[ip]->GetWorldVolume() == world)
	      {
		G4cout << "        o G4FastSimulationProcess: '" << fFSMPVector[ip]->GetProcessName() << "'" << G4endl;
		G4cout << "                      Attached to:";
		G4ParticleTable* particles = G4ParticleTable::GetParticleTable();
		for (G4int iParticle=0; iParticle<particles->entries(); iParticle++)
		  {
		    G4ParticleDefinition*  particle = particles->GetParticle(iParticle);
		    G4ProcessVector*      processes = particle->GetProcessManager()->GetProcessList();
		    if (processes->contains(fFSMPVector[ip])) {G4cout << " " << particle->GetParticleName(); findG4FSMP = true; particlesKnown.push_back(particle);}
		  }
		G4cout << G4endl;
	      }
	  if (!findG4FSMP) G4cout << "        o G4FastSimulationProcess: (none)" << G4endl;
	  // -- now display the regions in this world volume, with mother<->daughter link shown by indentation:
	  G4cout << "        o Region(s) and model(s) setup:" << G4endl;
	  DisplayRegion(worldRegion, 1, particlesKnown);
	}
    }
}


void G4GlobalFastSimulationManager::DisplayRegion(G4Region* region, G4int depth, std::vector<G4ParticleDefinition*>& particlesKnown) const
{
  G4String indent = "        ";
  for (G4int I=0; I<depth; I++) indent += "    ";
  G4cout << indent << "Region: `" << region->GetName() <<"'"  << G4endl;
  G4FastSimulationManager* fastSimManager = region->GetFastSimulationManager();
  if (fastSimManager)
    {
      indent += "    ";
      G4cout << indent << "Model(s):" << G4endl;
      indent += "    ";
      for (size_t im=0; im<fastSimManager->GetFastSimulationModelList().size(); im++)
	{
	  G4cout << indent << "`" << (fastSimManager->GetFastSimulationModelList())[im]->GetName() << "'";
	  G4cout << " ; applicable to:";
	  G4ParticleTable* particles = G4ParticleTable::GetParticleTable();
	  for (G4int iParticle=0; iParticle<particles->entries(); iParticle++)
	    {
	      if ((fastSimManager->GetFastSimulationModelList())[im]->IsApplicable(*(particles->GetParticle(iParticle))))
		{
		  G4cout << " " << particles->GetParticle(iParticle)->GetParticleName();
		  G4bool known(false);
		  for (size_t l=0; l<particlesKnown.size();l++) if(particlesKnown[l] == particles->GetParticle(iParticle)) {known = true; break;}
		  if (!known) G4cout << "[!!]";
		}
	    }
	  G4cout << G4endl;
	}
    }
  
  // -- all that to check mothership of "region"
  G4PhysicalVolumeStore* physVolStore = G4PhysicalVolumeStore::GetInstance();
  for (size_t ip=0; ip<physVolStore->size(); ip++)
    {
      G4VPhysicalVolume* physVol = (*physVolStore)[ip];
      if (physVol->GetLogicalVolume()->IsRootRegion())
	if (physVol->GetMotherLogical())
	  {
	    G4Region* thisVolMotherRegion = physVol->GetMotherLogical()->GetRegion();
	    if (thisVolMotherRegion == region)
	      DisplayRegion(physVol->GetLogicalVolume()->GetRegion(), depth+1, particlesKnown);
	  }
    }
}


// ----------------------------
// -- management methods : list
// ----------------------------

void G4GlobalFastSimulationManager::ListEnvelopes(const G4String&   aName,
						  listType        theType)
{
  if (theType == ISAPPLICABLE)
    {
      for (size_t ifsm=0; ifsm<ManagedManagers.size(); ifsm++) ManagedManagers[ifsm]->ListModels(aName);
      return;
    }
  
  if(aName == "all")
    {
      G4int titled = 0;
      for (size_t ifsm=0; ifsm<ManagedManagers.size(); ifsm++)
	{
	  if(theType == NAMES_ONLY)
	    {
	      if(!(titled++))
		G4cout << "Current Envelopes for Fast Simulation:\n";
	      G4cout << "   "; 
	      ManagedManagers[ifsm]->ListTitle();
	      G4cout << G4endl;
	    }
	  else ManagedManagers[ifsm]->ListModels();
	}
    }
  else
    {
      for (size_t ifsm=0; ifsm<ManagedManagers.size(); ifsm++)
	if(aName == ManagedManagers[ifsm]-> GetEnvelope()->GetName())
	  {
	    ManagedManagers[ifsm]->ListModels();
	    break;
	  }
    }
}

void G4GlobalFastSimulationManager::ListEnvelopes(const G4ParticleDefinition* aPD)
{
  for (size_t ifsm=0; ifsm<ManagedManagers.size(); ifsm++)
    ManagedManagers[ifsm]->ListModels(aPD);
}


G4VFastSimulationModel* G4GlobalFastSimulationManager::GetFastSimulationModel(const G4String&                   modelName,
									      const G4VFastSimulationModel* previousFound) const
{
  G4VFastSimulationModel* model = 0;
  // -- flag used to navigate accross the various managers;
  bool foundPrevious(false);
  for (size_t ifsm=0; ifsm<ManagedManagers.size(); ifsm++)
    {
      model = ManagedManagers[ifsm]->
              GetFastSimulationModel(modelName, previousFound, foundPrevious);
      if (model) break;
    }
  return model;
}
