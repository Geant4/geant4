////////////////////////////////////////////////////////////////////////////////
//
#include "MLPhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "MLGeometryConstruction.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ios.hh"
#include <iomanip>

#include "MLPhysicsListMessenger.hh"
#include "MLRunManager.hh"

#include "MLGeneralPhysics.hh"
#include "MLEMPhysics.hh"
#include "MLHEMPhysics.hh"
#include "MLLEEMPhysics.hh"
#include "MLMuonPhysics.hh"
#include "MLHadronPhysics.hh"
#include "MLIonPhysics.hh"
////////////////////////////////////////////////////////////////////////////////
//
MLPhysicsList::MLPhysicsList(): MLVModularPhysicsList(),
  scenario("hadron-em-ln")
{
  physListMessenger = new MLPhysicsListMessenger(this);
//
//
// The default cut value is set to 1.0mm.
//
  defaultCutValue  = 0.1*mm;
  cutForGamma      = defaultCutValue;
  cutForElectron   = defaultCutValue;
  cutForPositron   = defaultCutValue;

  // get the geometry, needed for region definition and cuts
  //
  G4RunManager *runManager = G4RunManager::GetRunManager();
  geometry = (MLGeometryConstruction*)(runManager->GetUserDetectorConstruction());
  //
  SetVerboseLevel(0);
}
////////////////////////////////////////////////////////////////////////////////
//
MLPhysicsList::~MLPhysicsList()
{
  delete physListMessenger;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::BuildList()
{
//
//
// Clean the register and register general physics processes.
//
  CleanAllPhysics();
  RegisterPhysics(new MLGeneralPhysics("general"));

  if (scenario == "em" || scenario == "EM") {
    RegisterPhysics(new MLEMPhysics("standard EM"));
    RegisterPhysics(new MLHEMPhysics("standard EM"));
    RegisterPhysics(new MLIonPhysics("Ion"));
  } else if (scenario == "leem" || scenario == "LEEM") {
    RegisterPhysics(new MLLEEMPhysics("l.e. EM"));
    RegisterPhysics(new MLHEMPhysics("l.e. EM"));
    RegisterPhysics(new MLIonPhysics("l.e. Ion"));
  } else if (scenario == "hadron-em-ln" || scenario == "HADRON-EM-LN") {
//    hadronic without EM and HP_neutron
    cutForGamma         = 100.*cm;
    cutForElectron      = 100.*cm;
    cutForPositron      = 100.*cm;
    RegisterPhysics(new MLEMPhysics("standard EM"));
    RegisterPhysics(new MLHEMPhysics("standard EM"));    
    RegisterPhysics(new MLHadronPhysics("Hadron-n"));
    RegisterPhysics(new MLIonPhysics("Ion"));
  } else if (scenario == "hadron-leem-ln" || scenario == "HADRON-LEEM-LN") {
    cutForGamma         = 100.*cm;
    cutForElectron      = 100.*cm;
    cutForPositron      = 100.*cm;
    RegisterPhysics(new MLLEEMPhysics("l.e. EM"));
    RegisterPhysics(new MLHEMPhysics("l.e. EM"));
    RegisterPhysics(new MLHadronPhysics("Hadron-n"));
    RegisterPhysics(new MLIonPhysics("l.e. Ion"));
  } else if (scenario == "hadron-em+ln" || scenario == "HADRON-EM+LN") {
    // hadronic without EM but with HP_neutron
    cutForGamma         = 100.*cm;
    cutForElectron      = 100.*cm;
    cutForPositron      = 100.*cm;
    RegisterPhysics(new MLEMPhysics("standard EM"));
    RegisterPhysics(new MLHEMPhysics("standard EM"));
    RegisterPhysics(new MLHadronPhysics("Hadron"));
    RegisterPhysics(new MLIonPhysics("Ion"));
  } else if (scenario == "hadron-leem+ln" || scenario == "HADRON-LEEM+LN") {
    // hadronic without EM but with HP_neutron
    cutForGamma         = 100.*cm;
    cutForElectron      = 100.*cm;
    cutForPositron      = 100.*cm;
    RegisterPhysics(new MLLEEMPhysics("l.e. EM"));
    RegisterPhysics(new MLHEMPhysics("l.e. EM"));
    RegisterPhysics(new MLHadronPhysics("Hadron+le"));
    RegisterPhysics(new MLIonPhysics("l.e. Ion"));
  } else if (scenario == "hadron+em+ln" || scenario == "HADRON+EM+LN") {
    RegisterPhysics(new MLEMPhysics("standard EM"));
    RegisterPhysics(new MLHEMPhysics("standard EM"));
    RegisterPhysics(new MLHadronPhysics("Hadron"));
    RegisterPhysics(new MLMuonPhysics("Muon"));
    RegisterPhysics(new MLIonPhysics("Ion"));
  } else if (scenario == "hadron+em-ln" || scenario == "HADRON+EM-LN") {
    RegisterPhysics(new MLEMPhysics("standard EM"));
    RegisterPhysics(new MLHEMPhysics("standard EM"));
    RegisterPhysics(new MLHadronPhysics("Hadron-n"));
    RegisterPhysics(new MLMuonPhysics("Muon"));
    RegisterPhysics(new MLIonPhysics("Ion"));
  } else if (scenario == "hadron+leem+ln" || scenario == "HADRON+LEEM+LN") {
    RegisterPhysics(new MLLEEMPhysics("l.e. EM"));
    RegisterPhysics(new MLHEMPhysics("l.e. EM"));
    RegisterPhysics(new MLHadronPhysics("Hadron"));
    RegisterPhysics(new MLMuonPhysics("Muon"));
    RegisterPhysics(new MLIonPhysics("l.e. Ion"));
  } else if (scenario == "hadron+leem-ln" || scenario == "HADRON+LEEM-LN") {
    RegisterPhysics(new MLLEEMPhysics("l.e. EM"));
    RegisterPhysics(new MLHEMPhysics("l.e. EM"));
    RegisterPhysics(new MLHadronPhysics("Hadron-n"));
    RegisterPhysics(new MLMuonPhysics("Muon"));
    RegisterPhysics(new MLIonPhysics("l.e. Ion"));
  } else if (scenario == "mars5gev" || scenario == "MARS5GEV") {
    RegisterPhysics(new MLEMPhysics("standard EM"));
    RegisterPhysics(new MLHEMPhysics("standard EM"));
    RegisterPhysics(new MLHadronPhysics("Mars5GeV"));
    RegisterPhysics(new MLMuonPhysics("Muon"));
    RegisterPhysics(new MLIonPhysics("Ion"));
  } else if (scenario == "mars5gev-ln" || scenario == "MARS5GEV-LN") {
    RegisterPhysics(new MLEMPhysics("standard EM"));
    RegisterPhysics(new MLHEMPhysics("standard EM"));
    RegisterPhysics(new MLHadronPhysics("Mars5GeV-ln"));
    RegisterPhysics(new MLMuonPhysics("Muon"));
    RegisterPhysics(new MLIonPhysics("Ion"));
  } else if (scenario == "mars5gev+pc" || scenario == "MARS5GEV+PC") {
    RegisterPhysics(new MLEMPhysics("standard EM"));
    RegisterPhysics(new MLHEMPhysics("standard EM"));
    RegisterPhysics(new MLHadronPhysics("Mars5GeV+pc"));
    RegisterPhysics(new MLMuonPhysics("Muon"));
    RegisterPhysics(new MLIonPhysics("Ion"));
  } else if (scenario == "binary" || scenario == "BINARY") {
    RegisterPhysics(new MLEMPhysics("standard EM"));
    RegisterPhysics(new MLHEMPhysics("standard EM"));
    RegisterPhysics(new MLHadronPhysics("Binary"));
    RegisterPhysics(new MLMuonPhysics("Muon"));
    RegisterPhysics(new MLIonPhysics("Ion"));
  } else if (scenario == "classic" || scenario == "CLASSIC") {
    RegisterPhysics(new MLEMPhysics("standard EM"));
    RegisterPhysics(new MLHEMPhysics("standard EM"));
    RegisterPhysics(new MLHadronPhysics("Classic"));
    RegisterPhysics(new MLMuonPhysics("Muon"));
    RegisterPhysics(new MLIonPhysics("Ion"));
  }

  if (verboseLevel > 1) {
    G4cout <<"The registered physics modules are: " <<G4endl;
    G4PhysConstVector::iterator itr;
    for (itr = physicsVector->begin(); itr!= physicsVector->end(); ++itr) {
      G4cout <<"       " <<(*itr)->GetPhysicsName() <<G4endl;
    }
    //    theParticleIterator->reset();
    // while ((*theParticleIterator)()) {
    //  G4ParticleDefinition* particle = theParticleIterator->value();
    //   G4cout <<particle->GetParticleName() <<G4endl;
    //  G4ProcessManager* pmanager = particle->GetProcessManager();
    //  G4ProcessVector* alist     = pmanager->GetProcessList();
    //   G4cout <<"Number of process = " <<alist->size() <<G4endl;
//      for (G4int i = 0; i < alist->size(); i++) {
//        G4cout << "  " << (alist[i])->GetProcessName() <<G4endl;
//      }
//    } 
  }  
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::SetScenario(G4String s)
{
  scenario = s;

  BuildList();
  MLRunManager* runManager =  static_cast <MLRunManager*>
    (MLRunManager::GetRunManager());
  //  runManager->CutOffHasNotBeenChanged(false);
  runManager->PhysicsListHasNotBeenChanged(false);
  runManager->Initialize();
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::SetGlobalCuts()
{
  if (verboseLevel>0) {
    G4cout <<"MLPhysicsList::SetGlobalCuts:" <<G4endl;
    G4cout <<"  GlobalCutLength :   "
	   <<G4BestUnit(defaultCutValue,"Length") <<G4endl;
    G4cout <<"  GammaCutLength :    "
	   <<G4BestUnit(cutForGamma,"Length") <<G4endl;
    G4cout <<"  ElectronCutLength : "
	   <<G4BestUnit(cutForElectron,"Length") <<G4endl;
    G4cout <<"  PositronCutLength :   "
	   <<G4BestUnit(cutForPositron,"Length") <<G4endl;
  } 
  //
  //
  // Set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma.
  //
  //  if (cutForPositron == cutForElectron == cutForGamma == defaultCutValue) {
  //  SetCutsWithDefault();
  // } else {
  SetCutValue(cutForGamma,"gamma");
  SetCutValue(cutForElectron,"e-");
  SetCutValue(cutForPositron,"e+");
  // }
  //  if (verboseLevel>0) DumpCutValuesTable();
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::SetCutsByRegion()
{
  /*  if (verboseLevel>0) {
    G4cout <<"MLPhysicsList::SetCutsByRegion:"  <<G4endl;
    G4cout <<"  GlobalCutLength   : "
           <<G4BestUnit(defaultCutValue,"Length") <<G4endl;
    G4cout <<"  GammaCutLength    : "
           <<G4BestUnit(cutForGamma,"Length") <<G4endl;
    G4cout <<"  ElectronCutLength : "
           <<G4BestUnit(cutForElectron,"Length") <<G4endl;
    G4cout <<"  ProtonCutLength   : "
           <<G4BestUnit(cutForProton,"Length") <<G4endl;
  }
  */
//
  if (verboseLevel>0) {
    G4RegionStore* aRegStore = G4RegionStore::GetInstance();
    G4cout <<"MLPhysicsList::SetCutsByRegion:"  <<G4endl;
    for ( std::vector<G4Region*>::iterator i = aRegStore->begin(); i != aRegStore->end(); i++)
      {
	//    if ((*i)->GetName() == name)
	G4cout <<" Region: " << (*i)->GetName() << G4endl;
	G4cout <<"  GammaCutLength    : "
           <<G4BestUnit((*i)->GetProductionCuts()->GetProductionCut("gamma"),"Length") << G4endl;
	G4cout <<"  e-CutLength    : "
           <<G4BestUnit((*i)->GetProductionCuts()->GetProductionCut("e-"),"Length") << G4endl;
	G4cout <<"  e+CutLength    : "
           <<G4BestUnit((*i)->GetProductionCuts()->GetProductionCut("e+"),"Length") << G4endl;
      }
  }
  // the cuts were set when the region was created
  //  if (verboseLevel>0) DumpCutValuesTable();
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::SetCutInRangeForRegion( G4double cut, G4String mName)
{
  if (cut > 0. && G4RegionStore::GetInstance()->GetRegion(mName) != 0) {
    G4ProductionCuts* cuts = new G4ProductionCuts();
    cuts->SetProductionCut(cut) ;
    G4RegionStore::GetInstance()->GetRegion(mName)->SetProductionCuts(cuts);
    if (verboseLevel > 1)
      G4cout <<"MLPhysicsList::SetCutInRangeForRegion: "
	     <<" cut : " <<cut <<" for region: " <<mName <<G4endl;

  } else {
    G4cout <<"MLPhysicsList::SetCutInRangeForRegion: "
           <<" cut : " <<cut <<" for region:" <<mName  <<G4endl;
    G4cout <<"--> Command rejected. Check the region name or cut." <<G4endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::SetCutInRangeForRegionForParticle( G4double cut,
  G4String mName, G4String part )
{
  if (cut > 0. && G4RegionStore::GetInstance()->GetRegion(mName) != 0 &&
    theParticleTable->FindParticle(part) != 0) {
    G4ProductionCuts* cuts = G4RegionStore::GetInstance()->GetRegion(mName)->GetProductionCuts();
    cuts->SetProductionCut(cut,part) ;
    G4RegionStore::GetInstance()->GetRegion(mName)->SetProductionCuts(cuts);    
    if (verboseLevel > 1) 
      G4cout <<"MLPhysicsList::SetCutInRangeForRegionForParticle: "
	     <<" cut : " <<cut <<" for region " <<mName <<" for " <<part <<G4endl;

  } else {
    G4cout <<"MLPhysicsList::SetCutInRangeForRegionForParticle: "
           <<" cut : " <<cut <<" for region " <<mName  <<" for " <<part <<G4endl;
    G4cout <<"--> Command rejected. Check the region and particle names and the cut! " <<G4endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::AddARegion(G4String s)
{
  G4Region* aReg = new G4Region(s);
  // set default G4ProductionCuts to defaultCutValue
  G4ProductionCuts* cuts = new G4ProductionCuts();
  cuts->SetProductionCut(defaultCutValue) ;
  aReg->SetProductionCuts(cuts);
  // G4RegionStore::GetInstance()->Register(aReg);  // why we dont need to register?
}

////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::DeleteARegion(G4String s)
{
  if(s == "") {
    G4RegionStore::GetInstance()->Clean();
  } else {
    G4RegionStore::GetInstance()->DeRegister( G4RegionStore::GetInstance()->GetRegion(s));
  }
}

////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::ListRegions()
{
  G4RegionStore* aRegStore = G4RegionStore::GetInstance();
  G4cout << " The defined regions are: " << G4endl;
  G4int k = 0;
  for ( std::vector<G4Region*>::iterator i = aRegStore->begin(); i != aRegStore->end(); i++) {
    k++;
    G4cout <<"   Region  " << k << ":" << (*i)->GetName() << G4endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::AddALayerToARegion(G4String s, G4int l)
{
  if (l < 1 || l > geometry->GetNbOfLayers()) {
    G4cout << " MLPhysicsList::AddALayerToARegion " << G4endl;
    G4cout << "   The layer " << l << " dose not exist. " << G4endl;
    return;
  }

  if (!G4RegionStore::GetInstance()->GetRegion(s))  {
    G4cout << " MLPhysicsList::AddALayerToARegion " << G4endl;
    G4cout << "   The region " << s << " dose not exist. " << G4endl;
    return;
  }
  G4Region* aReg = const_cast <G4Region*> (G4RegionStore::GetInstance()->GetRegion(s));
  G4LogicalVolume* aLog = const_cast <G4LogicalVolume*> ( geometry->GetLogicalLayer(l-1) ) ;
  aLog->SetRegion(aReg);
  aReg->AddRootLogicalVolume(aLog);

  //  G4cout << " MLPhysicsList::AddALayerToARegion: Number of maetial = " <<aReg->GetNumberOfMaterials() << G4endl;

}

////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::DeleteALayerFromARegion(G4String s, G4int l)
{
  if (l < 1 || l > geometry->GetNbOfLayers()) {
    G4cout << " MLPhysicsList::DeleteALayerFromARegion " << G4endl;
    G4cout << "   The layer " << l << " dose not exist. " << G4endl;
    return;
  }

  if (!G4RegionStore::GetInstance()->GetRegion(s))  {
    G4cout << " MLPhysicsList::DeleteALayerFromARegion " << G4endl;
    G4cout << "   The region " << s << " dose not exist. " << G4endl;
    return;
  }
  G4Region* aReg = const_cast <G4Region*> (G4RegionStore::GetInstance()->GetRegion(s));
  G4LogicalVolume* aLog = const_cast <G4LogicalVolume*> (geometry->GetLogicalLayer(l-1));
  aReg->RemoveRootLogicalVolume(aLog);
  
}

////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::ListTheLayersInARegion(G4String s)
{
  G4Region* aReg = G4RegionStore::GetInstance()->GetRegion(s);
  if (!aReg)  {
    G4cout << " MLPhysicsList::ListARegion " << G4endl;
    G4cout << "   The region " << s << " dose not exists. " << G4endl;
    return;
  }
  std::vector<G4LogicalVolume*>::iterator itr = aReg->GetRootLogicalVolumeIterator();
  G4cout << " The layers in region " << s <<" are:" << G4endl;
  while (*itr) {  
    G4cout <<(*itr)->GetName() <<G4endl;
    itr++;
  }
}

////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsList::SetCuts()
{

  // set global values first
  SetGlobalCuts() ;
  // deal with regions which have specified cuts
  SetCutsByRegion();
  //
  if (verboseLevel>0)   DumpCutValuesTable();
}

#ifndef USEHBOOK
////////////////////////////////////////////////////////////////////////////////
//
RPTofstream & operator << (RPTofstream &RPTFile, const MLPhysicsList &q)
{
  RPTFile <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
  RPTFile <<"Physics Definition:" <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
//
//  G4RunManager* runManager = G4RunManager::GetRunManager();
//  MLPhysicsList* physList = dynamic_cast <const MLPhysicsList*>
//    (runManager->GetUserPhysicsList());
//  MLPhysicsList *physList = &q;
  RPTFile <<"Simulation scenario : " <<q.GetScenario() <<G4endl;
  RPTFile <<"Applied default cuts: "
         <<G4BestUnit(q.GetDefaultCutValue(),"Length") <<G4endl;

  // now dumple the cuts, need to reproduce the following
  // RPTFile <<  G4ProductionCutsTable::GetProductionCutsTable()->DumpCouples();

  G4ProductionCutsTable* fCutsTable = G4ProductionCutsTable::GetProductionCutsTable();
  RPTFile <<"Production Cuts By Region:"  <<G4endl;
  G4RegionStore* aRegStore = G4RegionStore::GetInstance();
  for ( std::vector<G4Region*>::iterator i = aRegStore->begin(); i != aRegStore->end(); i++)
    {
      //    if ((*i)->GetName() == name)
      RPTFile <<" Region: " << (*i)->GetName() << G4endl;
      RPTFile <<"  GammaCutLength    : "
	     <<G4BestUnit((*i)->GetProductionCuts()->GetProductionCut("gamma"),"Length") << G4endl;
      RPTFile <<"  ElectronCutLength    : "
	     <<G4BestUnit((*i)->GetProductionCuts()->GetProductionCut("e-"),"Length") << G4endl;
      RPTFile <<"  PositronCutLength    : "
	     <<G4BestUnit((*i)->GetProductionCuts()->GetProductionCut("e+"),"Length") << G4endl;

      // now energy cuts
      RPTFile <<"  Cuts in Energy    : " << G4endl;
      G4Region* aRegion = (*i);
      //      std::vector<G4LogicalVolume*>::iterator itr = aRegion->GetRootLogicalVolumeIterator();
      //  RPTFile <<"   The layers in the region are:" << G4endl;
      // while (*itr) {  
      //	RPTFile <<"       " <<(*itr)->GetName() <<G4endl;
      //	itr++;
      //}
      std::vector<G4Material*>::const_iterator mItr = aRegion->GetMaterialIterator();
      size_t nMaterial = aRegion->GetNumberOfMaterials();
      // RPTFile <<    "    Number of material = " << nMaterial << G4endl;
      for(size_t iMate=0;iMate<nMaterial;iMate++, mItr++){
	RPTFile << "    in " << (*mItr)->GetName() << ": " << G4endl;
	G4int cIdx = fCutsTable->GetCoupleIndex((*mItr),aRegion->GetProductionCuts()); 
	
	RPTFile << "        gamma " << G4BestUnit((*(fCutsTable->GetEnergyCutsVector(0)))[cIdx],"Energy")
	       << "    e- " << G4BestUnit((*(fCutsTable->GetEnergyCutsVector(1)))[cIdx],"Energy")
	       << "    e+ " << G4BestUnit((*(fCutsTable->GetEnergyCutsVector(2)))[cIdx],"Energy");
	RPTFile << G4endl;
      }
    }

  //    G4cout << " gamma " << G4BestUnit((*(energyCutTable[0]))[aCouple->GetIndex()],"Energy")

  /*
  G4VPhysicalVolume *tempPV      = NULL;
  G4PhysicalVolumeStore *PVStore = NULL;
  PVStore                        = G4PhysicalVolumeStore::GetInstance();

  std::vector <G4Material*> mList ;
  G4Material* aMat = NULL;
  G4int k = 0;
  G4int i; 
  for ( i=0; i<G4int(PVStore->size());i++) {
    tempPV = (*PVStore)[i];
    aMat   = tempPV->GetLogicalVolume()->GetMaterial();
    k      = count(mList.begin(),mList.end(),aMat);
    if (k==0) mList.push_back(aMat);
  }

  G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();

  G4String name1[3] = {"gammas", "electrons", "positrons"};
  G4String name2[3] = {"gamma",  "e-", "e+"};
  G4ParticleDefinition* particle = NULL;
  G4double* theKineticEnergyCuts = NULL;
  G4double* theRangeCuts         = NULL;
  i                              = 0;

  for (G4int j = 0; j < 3; j++) {
    RPTFile <<" Production cuts for " <<name1[j] <<" : " <<G4endl;
    particle             = theParticleTable->FindParticle(name2[j]);
    theKineticEnergyCuts = particle->GetEnergyCuts();
    theRangeCuts         = particle->GetLengthCuts();
    i                    = 0;
    while (i<G4int(mList.size())) {
      RPTFile <<"   in material: " <<mList[i]->GetName();
      for (size_t idx=0; idx<materialTable->size(); idx++) {
        if ((*materialTable)[idx]->GetName() == mList[i]->GetName()) {
          RPTFile <<"; cut in range = "
                 <<G4BestUnit(theRangeCuts[idx],"Length")
                 << "; cut in energy = " 
                 <<G4BestUnit(theKineticEnergyCuts[idx],"Energy")
                 <<G4endl;
        }
      }
      i++;
    }
  }
  */
  return RPTFile;
}
#endif

////////////////////////////////////////////////////////////////////////////////
