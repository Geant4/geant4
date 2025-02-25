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
// G4GeometrySampler
// ----------------------------------------------------------------------

#include "G4GeometrySampler.hh"

#include "G4VIStore.hh"
#include "G4WeightWindowStore.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ImportanceConfigurator.hh"
#include "G4WeightWindowConfigurator.hh"
#include "G4WeightCutOffConfigurator.hh"
#include "G4TransportationManager.hh"

 G4GeometrySampler::
 G4GeometrySampler(G4VPhysicalVolume *world, const G4String& particlename)
  : fParticleName(particlename),
    fWorld(world)
{
}

 G4GeometrySampler::
 G4GeometrySampler(const G4String& worldName, const G4String &particlename)
  : fParticleName(particlename),
    fWorldName(worldName)
{
  fWorld = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
}

G4GeometrySampler::~G4GeometrySampler()
{
  // ClearSampling();
}

void G4GeometrySampler::ClearSampling()
{
  delete fImportanceConfigurator; fImportanceConfigurator = nullptr;
  delete fWeightWindowConfigurator; fWeightWindowConfigurator = nullptr;
  delete fWeightCutOffConfigurator; fWeightCutOffConfigurator = nullptr;
  fIStore = nullptr;
  fConfigurators.clear();
  fIsConfigured = false;
}

G4bool G4GeometrySampler::IsConfigured() const
{
  G4bool isconf = false;
  if (fIsConfigured)
  {
   G4cout << "WARNING - G4GeometrySampler::IsConfigured()"
          << "          Some initialization exists, use ClearSampling()"
          << "          before a new initialization !" << G4endl;
   isconf = true;
  }
  return isconf;
}

// void G4GeometrySampler::PrepareScoring(G4VScorer *scorer)
// {
//   G4cout << " preparing scoring configurator " << G4endl;
//   G4cout << G4endl;
//   G4cout << G4endl;
//   G4cout << G4endl;
//   G4cout << " new fWorld Name: " << fWorld->GetName() << G4endl;
//   G4cout << G4endl;
//   G4cout << G4endl;
//   G4cout << G4endl;
//   fScoreConfigurator = new G4ScoreConfigurator(fWorld, fParticleName, *scorer, paraflag);
//   G4cout << " configured scoring " << G4endl;
//   if (!fScoreConfigurator)
//   {
//     G4Exception("G4GeometrySampler::PrepareScoring()",
//                 "FatalError", FatalException,
//                 "Failed allocation of G4ScoreConfigurator !");
//   }
// }

void
G4GeometrySampler::PrepareImportanceSampling(G4VIStore* istore,
                                           const G4VImportanceAlgorithm  *ialg)
{
  G4cout << "G4GeometrySampler:: preparing importance sampling WorldName is " << fWorldName << G4endl;
  fIStore = istore;
  //  G4cout << "G4GeometrySampler:: creating istore, worldVolume: " << fWorld->GetName() << G4endl;

  fImportanceConfigurator =
    new G4ImportanceConfigurator(&istore->GetWorldVolume(), fParticleName, *fIStore, ialg, paraflag);
  //    new G4ImportanceConfigurator(fWorld, fParticleName, *fIStore, ialg, paraflag);
  fImportanceConfigurator->SetWorldName(fWorldName);

  if (!fImportanceConfigurator)
  {
    G4Exception("G4GeometrySampler::PrepareImportanceSampling()",
                "FatalError", FatalException,
                "Failed allocation of G4ImportanceConfigurator !");
  }
}

void
G4GeometrySampler::PrepareWeightRoulett(G4double wsurvive, 
                                        G4double wlimit,
                                        G4double isource)
{
  //  fGCellFinder = new G4GCellFinder(fWorld);
  G4cout << "G4GeometrySampler:: preparing weight roulette" << G4endl;
  //  fGCellFinder = new G4GCellFinder();
//   if (!fGCellFinder)
//   {
//     G4Exception("G4GeometrySampler::PrepareWeightRoulett()",
//                 "FatalError", FatalException,
//                 "Failed allocation of G4GCellFinder !");
//   }
  
  fWeightCutOffConfigurator = 
    new G4WeightCutOffConfigurator(fWorld, fParticleName,
                                   wsurvive,
                                   wlimit,
                                   isource,
                                   fIStore,
				   paraflag);
				   //*fGCellFinder, paraflag);
  if (!fWeightCutOffConfigurator)
  {
    G4Exception("G4GeometrySampler::PrepareWeightRoulett()",
                "FatalError", FatalException,
                "Failed allocation of G4WeightCutOffConfigurator !");
  }
}

void
G4GeometrySampler::PrepareWeightWindow(G4VWeightWindowStore *wwstore,
                                       G4VWeightWindowAlgorithm *wwAlg,
                                       G4PlaceOfAction placeOfAction)
{

  G4cout << "G4GeometrySampler:: preparing weight window" << G4endl;

  fWWStore = wwstore;
  
  fWeightWindowConfigurator =
    new G4WeightWindowConfigurator(&wwstore->GetWorldVolume(), fParticleName,
                                    *fWWStore,
                                    wwAlg,
                                    placeOfAction, paraflag);
}

void G4GeometrySampler::Configure()
{
  if (!IsConfigured())
  {
    fIsConfigured = true;

    if (fImportanceConfigurator)
    {
      fConfigurators.push_back(fImportanceConfigurator);
    }
    if (fWeightWindowConfigurator)
    {
      fConfigurators.push_back(fWeightWindowConfigurator);
    }
  }

#ifdef G4MULTITHREADED
  G4cout << " make sure AddProcess() is invoked for biasing!!! " << G4endl;
#else
  AddProcess();
#endif

  return;
}

void G4GeometrySampler::AddProcess()
{

  G4VSamplerConfigurator *preConf = nullptr;
  for (auto it = fConfigurators.cbegin();
       it != fConfigurators.cend(); ++it)
    {
      G4VSamplerConfigurator *currConf =*it;
      currConf->Configure(preConf);
      preConf = *it;
    }
  if (fWeightCutOffConfigurator != nullptr)
    {
      fWeightCutOffConfigurator->Configure(nullptr);
    }

  return;
}

void G4GeometrySampler::SetParallel(G4bool para)
{
  paraflag = para;
}

void G4GeometrySampler::SetWorld(const G4VPhysicalVolume* World)
{
  fWorld = World;
}

void G4GeometrySampler::SetParticle(const G4String& particlename)
{
  fParticleName = particlename;
}
