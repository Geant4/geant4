#include "Hall.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "HallMessenger.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"



//#include "stdafx.h"
/* globalni danni za procesa:
   1 - Importance store - za syhranjavane na obemite na paralelnata geometrija
   2 - masiv sys Scorer-i - za broene na chasticite.
   Importance store se promenja pri promjana na sloevete. Vseki pyt toj se izstriva i se poluchava nanovo.
   World volume-to se pravi predvaritelno: Physical volume-to se iska, no sled syzdavaneto na geometrijata
   Za vseki sloj se syzdavat po opredelen broj sloeve ot paralelnata geometrija. Maj njama da mojem da izpolzwame
   edni i syshti physical volumi tyj che te se razlichavat. No Logical volumite i solidite gi pravim edni i syshti
*/
//M #include "G4Scorer.hh"
//M #include "G4ScoreTable.hh"
//M #include "G4ParallelImportanceScoreSampler.hh"

#include "G4ParallelImportanceSampler.hh"


//M #include "G4ParallelScoreSampler.hh"
#include "G4IStore.hh"

//M enum __PARTICLES_FOR_SCORER {NEUTRON,PROTON,DEUTERON,TRITON,HELIUM_3,ALPHA};

static G4VSolid * pPWorldSolid;
static G4VSolid** pPTubesSolid;
static G4VSolid*  pPDetectorSolid;
static G4LogicalVolume * pPWorldLog;
static G4LogicalVolume** pPTubesLog;
static G4LogicalVolume*  pPDetectorLog;
static G4PVPlacement* pPWorld;
static G4PVPlacement ** pPTubes;
static G4PVPlacement*  pPDetector;
static unsigned nTubesNumber=0;
//M static bool pIsPartScorer[6];
//M static G4Scorer * pPartScorer[6];
//M static G4ParallelImportanceScoreSampler* pNeutronSampler;
static G4ParallelImportanceSampler* pNeutronSampler;
//M static G4ParallelScoreSampler* pRestSamplers[6];
static bool bInit=false;
static G4IStore* pPGeom;
static unsigned nParallelLayers=1;
static char* pArrNames[6] = {"neutron",
                             "proton",
                             "deuteron",
                             "triton",
                             "He3",
                             "alpha"};

static void ConstructParallelBase();
void Hall::ConstructParallelLayers()
{
  unsigned nLayers = m_pHall->pShield->nLayerNum,nCurrLayer,nCurrLayer1;
  double dThick,dTotalThick=0;
  if(nLayers==0) return;
  pPTubesSolid = new G4VSolid*[nLayers*nParallelLayers];
  pPTubesLog = new G4LogicalVolume*[nLayers*nParallelLayers];
  pPTubes = new G4PVPlacement*[nLayers*nParallelLayers];
  for(nCurrLayer=0;nCurrLayer<nLayers;nCurrLayer++){
    dThick = m_pHall->pShield->pLayersThick[nCurrLayer];
    dThick /= nParallelLayers;
    for(nCurrLayer1=0;nCurrLayer1<nParallelLayers;nCurrLayer1++){
      pPTubesSolid[nCurrLayer*nParallelLayers+nCurrLayer1] = new G4Tubs(G4String("Parallel layer ")+G4String((char)(nCurrLayer*nParallelLayers+nCurrLayer1)+'0'),0,110*cm,dThick/2.,0,2*pi);
      pPTubesLog[nCurrLayer*nParallelLayers+nCurrLayer1] = new G4LogicalVolume(pPTubesSolid[nCurrLayer*nParallelLayers+nCurrLayer1],G4Material::GetMaterial("Concrete"),G4String("Parallel logical layer ")+G4String((char)(nCurrLayer*nParallelLayers + nCurrLayer1)+'0'));
      pPTubes[nCurrLayer*nParallelLayers+nCurrLayer1] = new G4PVPlacement(0,G4ThreeVector(0,0,76*cm+dTotalThick + dThick/2.),G4String("Parallel physical layer ")+G4String((char)(nCurrLayer*nParallelLayers+nCurrLayer1)+'0'),pPTubesLog[nCurrLayer*nParallelLayers+nCurrLayer1],pPWorld,false,0);
      pPGeom->AddImportanceRegion(1L<<((unsigned)(nCurrLayer*nParallelLayers+nCurrLayer1)),*(pPTubes[nCurrLayer*nParallelLayers+nCurrLayer1]));

      nTubesNumber++;
      G4cout<<76*cm+dTotalThick+dThick/2.<<" mm is the origin of "<<G4String("Parallel layer ")+G4String((char)(nCurrLayer*nParallelLayers+nCurrLayer1)+'0')<<G4endl;
      dTotalThick += dThick;
    }
  }
  pPDetectorSolid = new G4Tubs("Parallel detector solid",0,110*cm,6.35*cm,0,2*pi);
  pPDetectorLog = new G4LogicalVolume(pPDetectorSolid,G4Material::GetMaterial("Air"),G4String("Parallel Detector Logical"));
  pPDetector = new G4PVPlacement(0,G4ThreeVector(0,0,82.35*cm+dTotalThick),"Parallel detector",pPDetectorLog,pPWorld,false,0);
  pPGeom->AddImportanceRegion(1L<<((unsigned)(nCurrLayer*nParallelLayers)),*pPDetector);
  Update();
}

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

static void InitializeParallelGeometry()
{
  return; // hpw, to get it to run for mars
  //Tuk ako e imalo scorer toj trjabva da se izstrie
  //M  pPartScorer[0] = new G4Scorer;
  //M  pNeutronSampler = new G4ParallelImportanceScoreSampler(*pPGeom,*pPartScorer[0],"neutron");
  pNeutronSampler = new G4ParallelImportanceSampler(*pPGeom, "neutron");
  pNeutronSampler->Initialize();
  //M  for(unsigned i=1;i<6;i++){
  //M    pPartScorer[i] = new G4Scorer;
  //M    pRestSamplers[i] = new G4ParallelScoreSampler(*pPWorld,pArrNames[i],*pPartScorer[i]);
    //pRestSamplers[i]->Initialize();
  //M  }
  bInit = true;
}
static void DeleteParallelGeometry()
{
  return; // hpw to avoid crash
  if(pNeutronSampler) delete pNeutronSampler;
  pNeutronSampler=NULL;
  unsigned i;
  //M  for(i=1;i<6;i++){
  //M    if(pRestSamplers[i]) delete pRestSamplers[i];
  //M    pRestSamplers[i]=NULL;
  //M  }
  bInit = false;
  G4ParticleDefinition* pPart;
  G4ProcessManager* pMan;
  G4ProcessVector* pList;
  unsigned index,MaxInd;
  for(i=1;i<6;i++){
    pPart = G4ParticleTable::GetParticleTable()->FindParticle(pArrNames[i]);
    if(pPart!=NULL){
      pMan = pPart->GetProcessManager();
      pList = pMan->GetProcessList();
      MaxInd = pList->entries();
      for(index=0;index<MaxInd;index++){
        if((pList->operator()(index))->GetProcessName()==G4String("ParallelTransport")){
          pMan->RemoveProcess(index);
          break;
        }
      }

    }
  }
}

void Hall::PostLayerInit()

{
  if(nTubesNumber==0){
    G4cout<<"WARNING: No layers are added"<<G4endl;
    return;
  }
  InitializeParallelGeometry();
}

static void ClearImportance(bool ClearSolids = true)
{
  unsigned i;
  if(bInit)
    DeleteParallelGeometry();
  else{
    if(pPGeom) delete pPGeom;
    if(ClearSolids){
      for(i=0;i<nTubesNumber;i++){
        if(pPDetectorSolid) delete pPDetectorSolid;
        pPDetectorSolid = NULL;
        if(pPDetectorLog) delete pPDetectorLog;
        pPDetectorLog = NULL;
        if(pPDetector) delete pPDetector;
        pPDetector = NULL;
        if(pPTubesSolid && pPTubesSolid[i]) delete pPTubesSolid[i];
        if(pPTubesLog && pPTubesLog[i]) delete pPTubesLog[i];
        if(pPTubes && pPTubes[i]) delete pPTubes[i];
      }
    }
    if(pPTubes) delete pPTubes;
    if(ClearSolids)  if(pPWorld) delete pPWorld;
    if(ClearSolids)  if(pPWorldLog) delete pPWorldLog;
    if(pPTubesLog) delete pPTubesLog;
    if(ClearSolids)  if(pPWorldSolid) delete pPWorldSolid;
    if(pPTubesSolid) delete pPTubesSolid;
  }
  //M  for(i=0;i<6;i++){
  //M    if(pIsPartScorer[i]) delete pPartScorer[i];
  //M    pPartScorer[i] = NULL;

  //M  }
  pNeutronSampler = NULL;
  pPTubes = NULL;
  pPWorld = NULL;
  pPTubesLog = NULL;
  pPWorldLog = NULL;
  pPWorldSolid = NULL;
  pPTubesSolid = NULL;
  pPGeom=NULL;
}

Hall::Hall()
{
  m_pHall = NULL;
  m_pLogicalHall = NULL;
  m_pPhysHall = NULL;
  pPWorldSolid = NULL;
  pPTubesSolid = NULL;
  pPWorldLog = NULL;
  pPTubesLog = NULL;
  pPWorld = NULL;
  pPTubes = NULL;
  pPGeom = NULL;
  pPDetectorSolid = NULL;
  pPDetectorLog = NULL;
  pPDetector = NULL;

  pNeutronSampler = NULL;
  //M  for(unsigned i=0;i<6;i++){
  //M    pIsPartScorer[i] = false;
  //M    pPartScorer[i] = NULL;
  //M    pRestSamplers[i] = NULL;
  //M  }
  m_pMessenger = new HallMessenger(this);
}
Hall::~Hall()
{
  G4int i,j;
  ClearImportance(false);
  if(m_pPhysHall){
    delete m_pPhysHall->pMother;
    delete m_pPhysHall->pTarget;
    delete m_pPhysHall->pDeadConcreteA;
    delete m_pPhysHall->pDeadConcreteB;
    delete m_pPhysHall->pIronShieldA;
    delete m_pPhysHall->pIronShieldB;
    delete m_pPhysHall->pDetector;
    delete m_pPhysHall->pDetector20;
    delete m_pPhysHall->pDetector40;
    delete m_pPhysHall->pGhostDetector;
    if(m_pPhysHall->pShield){
      j = m_pPhysHall->pShield->nLayerNum;
      for(i=0;i<j;i++)
        delete m_pPhysHall->pShield->pLayers[i];
      delete m_pPhysHall->pShield->pVLayer;
    }
    delete m_pPhysHall;
  }
  if(m_pLogicalHall){
    delete m_pLogicalHall->pMother;
    delete m_pLogicalHall->pTarget;
    delete m_pLogicalHall->pDeadConcreteA;
    delete m_pLogicalHall->pDeadConcreteB;
    delete m_pLogicalHall->pIronShieldA;
    delete m_pLogicalHall->pIronShieldB;
    delete m_pLogicalHall->pDetector;
    delete m_pLogicalHall->pDetector20;
    delete m_pLogicalHall->pDetector40;
    delete m_pLogicalHall->pGhostDetector;
    if(m_pLogicalHall->pShield){
      j = m_pLogicalHall->pShield->nLayerNum;
      for(i=0;i<j;i++)
        delete m_pLogicalHall->pShield->pLayers[i];
      delete m_pLogicalHall->pShield->pVLayer;
      delete m_pLogicalHall->pShield->pLayerMat;
    }
    delete m_pLogicalHall;
  }
  if(m_pHall){
    delete m_pHall->pMother;
    delete m_pHall->pTarget;
    delete m_pHall->pDetector;
    delete m_pHall->pDetector20;
    delete m_pHall->pDetector40;
    delete m_pHall->pGhostDetector;
    delete m_pHall->pDeadConcreteA;
    delete m_pHall->pDeadConcreteB;
    delete m_pHall->pIronShieldA;
    delete m_pHall->pIronShieldB;
    if(m_pHall->pShield){
      j = m_pHall->pShield->nLayerNum;
      for(i=0;i<j;i++)
        delete m_pHall->pShield->pLayers[i];
      delete m_pHall->pShield->pVLayer;
      delete m_pHall->pShield->pLayersThick;
      delete m_pHall->pShield->pLayersColimator;
    }
    delete m_pHall;
  }
  delete m_pMessenger;
}
    
void Hall::CreateLayers(G4int nLayerNum)
{
  G4int i;
  size_t minSize;
  void** pNew;
  if(nLayerNum<1) return;
  if(m_pHall->pShield){
    if(m_pHall->pShield->nLayerNum==nLayerNum) return;
  }
  else
    {
      m_pHall->pShield = new SHIELDTYPE;
      m_pHall->pShield->nLayerNum = 0;
      m_pHall->pShield->pVLayer = NULL;
      m_pLogicalHall->pShield = new LOGICALSHIELD;
      m_pLogicalHall->pShield->nLayerNum = 0;
      m_pLogicalHall->pShield->pVLayer = NULL;
      m_pPhysHall->pShield = new PHYSICALSHIELD;
      m_pPhysHall->pShield->nLayerNum = 0;
      m_pPhysHall->pShield->pVLayer = NULL;
    }

  minSize = (m_pHall->pShield->nLayerNum < nLayerNum) ? 

    m_pHall->pShield->nLayerNum : nLayerNum;

  if(minSize)
    DeleteLayersDaughters();
  ClearImportance();
  //Compousing Layer's solid
  pNew = new (void*)[nLayerNum];
  if(minSize){
    for(unsigned int j=0;j<minSize;j++)
      ((G4Tubs**)pNew)[j] = m_pHall->pShield->pLayers[j];
  };
  if(minSize==(unsigned)nLayerNum){
    for(i=minSize;i<m_pHall->pShield->nLayerNum;i++)
      delete m_pHall->pShield->pLayers[i];
  }
  if(m_pHall->pShield->nLayerNum!=0)
    delete m_pHall->pShield->pLayers;

  m_pHall->pShield->pLayers = (G4Tubs**)pNew;
  (void*)pNew = new(G4double)[nLayerNum];
  if(minSize){
    for(unsigned j=0;j<minSize;j++)
      ((G4double*)pNew)[j] = m_pHall->pShield->pLayersThick[j];
    delete m_pHall->pShield->pLayersThick;
  }
  m_pHall->pShield->pLayersThick = (G4double*)pNew;
  (void*)pNew = new (G4bool)[nLayerNum];
  if(minSize){
    for(unsigned j=0;j<minSize;j++)
      ((G4bool*)pNew)[j] = m_pHall->pShield->pLayersColimator[j];
    delete m_pHall->pShield->pLayersColimator;
  }
  m_pHall->pShield->pLayersColimator = (G4bool*)pNew;
  CompouseSolidLayer(minSize,nLayerNum);
  m_pHall->pShield->nLayerNum = nLayerNum;

  //Compousing Layer's LogicalVolumes;
  pNew = new (void*)[nLayerNum];
  if(minSize){
    for(unsigned j=0;j<minSize;j++)
      ((G4Material**)pNew)[j] = m_pLogicalHall->pShield->pLayerMat[j];
    delete m_pLogicalHall->pShield->pLayerMat;
  }

  m_pLogicalHall->pShield->pLayerMat = (G4Material**)pNew;
  pNew = new (void*)[nLayerNum];
  if(minSize)
    for(unsigned j=0;j<minSize;j++)
      ((G4LogicalVolume**)pNew)[j] = m_pLogicalHall->pShield->pLayers[j];
  if(minSize==(unsigned)nLayerNum){
    for(i=minSize;i<m_pLogicalHall->pShield->nLayerNum;i++)
      delete m_pLogicalHall->pShield->pLayers[i];
  }
  if(minSize)
    delete m_pLogicalHall->pShield->pLayers;
  m_pLogicalHall->pShield->pLayers = (G4LogicalVolume**)pNew;
  CompouseLogicalLayer(minSize,nLayerNum);
  m_pLogicalHall->pShield->nLayerNum = nLayerNum;

  //Compousing Physical Layer
  pNew = new (void*)[nLayerNum];
  if(minSize)
    for(unsigned j=0;j<minSize;j++)
      ((G4VPhysicalVolume**)pNew)[j] = m_pPhysHall->pShield->pLayers[j];
  if(minSize == (unsigned)nLayerNum){
    for(i=minSize+1;i<m_pPhysHall->pShield->nLayerNum;i++)
      delete m_pPhysHall->pShield->pLayers[i];
  }
  if(minSize)
    delete m_pPhysHall->pShield->pLayers;
  m_pPhysHall->pShield->pLayers = (G4VPhysicalVolume**)pNew;
  CompousePhysicalLayer(minSize,nLayerNum);

  m_pPhysHall->pShield->nLayerNum = nLayerNum;

  DumpLayers();
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->DefineWorldVolume(m_pPhysHall->pMother);

  ConstructParallelBase();
  ConstructParallelLayers();
  PostLayerInit();
  Update();
}
void Hall::DeleteLayersDaughters()
{
  G4int i = m_pLogicalHall->pShield->pVLayer->GetNoDaughters();
  for(G4int j=0;j<i;j++)
    m_pLogicalHall->pShield->pVLayer->RemoveDaughter(m_pLogicalHall->pShield->pVLayer->GetDaughter(j));
  m_pLogicalHall->pMother->RemoveDaughter(m_pPhysHall->pShield->pVLayer);
}
void Hall::CompouseSolidLayer(G4int from,G4int size)
{
  G4int i;
  G4double dRealSize=0;
  for(i=0;i<size;i++){
    if(i<from)
      dRealSize+=m_pHall->pShield->pLayersThick[i];
    else{
      dRealSize +=25.*cm;
      m_pHall->pShield->pLayersThick[i]=25.*cm;
      m_pHall->pShield->pLayers[i] = new G4Tubs(G4String("Layer")+G4String((char)i+'0'),
                                                0.,59.*cm,12.5*cm,0.,2*pi);
      m_pHall->pShield->pLayersColimator[i] = false;
    }
  }
  if(dRealSize>150.*cm)
    G4cerr<<"Warning: Total layer thickness is larger than is allowed"<<G4endl;
  if(m_pHall->pShield->pVLayer)
    delete m_pHall->pShield->pVLayer;
  m_pHall->pShield->pVLayer = new G4Tubs(G4String("Layer"),0.,59.*cm,
                                         dRealSize/2.,0,2*pi);
}
void Hall::CompouseLogicalLayer(G4int from,G4int size)

{
  G4int i;
  static G4VisAttributes vis1;
  static G4VisAttributes vis2;
  if(m_pLogicalHall->pShield->pVLayer)
    delete m_pLogicalHall->pShield->pVLayer;
  m_pLogicalHall->pShield->pVLayer = new G4LogicalVolume(m_pHall->pShield->pVLayer,
                                                         G4Material::GetMaterial(G4String("Vacuum")),
                                                         "Locigal Layer");
  vis1.SetVisibility(false);
  vis1.SetDaughtersInvisible(false);
  m_pLogicalHall->pShield->pVLayer->SetVisAttributes(vis1);

  vis2.SetVisibility(true);
  vis2.SetColour(92./255.,92./255.,92./255.);
  for(i=from;i<size;i++){
    m_pLogicalHall->pShield->pLayerMat[i] = G4Material::GetMaterial("Concrete");
    m_pLogicalHall->pShield->pLayers[i] = new G4LogicalVolume(m_pHall->pShield->pLayers[i],
                                                              m_pLogicalHall->pShield->pLayerMat[i],
                                                              G4String("Logical Layer")+G4String((char)i+'0'));
    m_pLogicalHall->pShield->pLayers[i]->SetVisAttributes(vis2);
  }
}
void Hall::CompousePhysicalLayer(G4int from,G4int size)
{
  G4double ztransl=0.;
  G4double TotalThick=0.;
  G4int i;
  if(m_pPhysHall->pShield->pVLayer)
    delete m_pPhysHall->pShield->pVLayer;
  for(i=0;i<m_pHall->pShield->nLayerNum;i++) TotalThick += m_pHall->pShield->pLayersThick[i];
  m_pPhysHall->pShield->pVLayer = new G4PVPlacement(0,G4ThreeVector(0,0,76.*cm+TotalThick/2.),
                                                    G4String("Physical Layer"),
                                                    m_pLogicalHall->pShield->pVLayer,
                                                    m_pPhysHall->pMother,
                                                    false,
                                                    0);
  m_pPhysHall->pDetector->SetTranslation(G4ThreeVector(0,0,82.35*cm+TotalThick));
  m_pPhysHall->pDetector20->SetTranslation(G4ThreeVector(20.*cm,0,82.35*cm+TotalThick));
  m_pPhysHall->pDetector40->SetTranslation(G4ThreeVector(40.*cm,0.,82.35*cm+TotalThick));
  ztransl = -TotalThick/2.;
  for(i=0;i<size;i++){
    if(i<from){
      m_pLogicalHall->pShield->pVLayer->AddDaughter(m_pPhysHall->pShield->pLayers[i]);
      m_pPhysHall->pShield->pLayers[i]->SetTranslation(G4ThreeVector(0,0,ztransl+m_pHall->pShield->pLayersThick[i]/2.));
      ztransl += m_pHall->pShield->pLayersThick[i];
      continue;
    }
    m_pPhysHall->pShield->pLayers[i] = new G4PVPlacement(0,
                                                         G4ThreeVector(0,0,m_pHall->pShield->pLayersThick[i]/2. +ztransl),
                                                         G4String("Physical Layer ")+G4String((char)i+'0'),
                                                         m_pLogicalHall->pShield->pLayers[i],
                                                         m_pPhysHall->pShield->pVLayer,
                                                         false,
                                                         0);
    ztransl+=m_pHall->pShield->pLayersThick[i];
  }
  if(TotalThick>150.*cm)
    G4cerr<<"Warning: Sum of layers thick is larger than 150 cm: "<<TotalThick/cm<<G4endl;
}

void Hall::SetAsColimator(int iLayer)
{
  if(m_pHall->pShield->nLayerNum < iLayer) return;

  if(m_pHall->pShield->pLayersColimator[iLayer-1]) return;
  m_pHall->pShield->pLayersColimator[iLayer-1] = true;
  delete m_pHall->pShield->pLayers[iLayer-1];
  m_pHall->pShield->pLayers[iLayer-1] = new G4Tubs(G4String("Layer")+G4String(iLayer+'0'),
                                                 5.45*cm,59.*cm,
                                                 m_pHall->pShield->pLayersThick[iLayer-1]/2.,
                                                 0.,2.*pi);
  m_pLogicalHall->pShield->pLayers[iLayer-1]->SetSolid(m_pHall->pShield->pLayers[iLayer-1]);
  Update();
}

void Hall::SetAsShield(int iLayer)
{
  if(m_pHall->pShield->nLayerNum < iLayer) return;
  if(!m_pHall->pShield->pLayersColimator[iLayer-1]) return;
  m_pHall->pShield->pLayersColimator[iLayer-1] = false;
  delete m_pHall->pShield->pLayers[iLayer-1];
  m_pHall->pShield->pLayers[iLayer-1] = new G4Tubs(G4String("Layer")+G4String(iLayer+'0'),
                                                 0.,59.*cm,

                                                 m_pHall->pShield->pLayersThick[iLayer-1]/2.,
                                                 0.,2.*pi);
  m_pLogicalHall->pShield->pLayers[iLayer-1]->SetSolid(m_pHall->pShield->pLayers[iLayer-1]);
  Update();
}

void Hall::LayerMat(G4String szMat,int iLayer)
{
  static G4VisAttributes vis;
  G4Material* pMat = G4Material::GetMaterial(szMat);
  if(!pMat){
    G4cerr<<"Error: Material: "<<szMat<<" does not exist"<<G4endl;
    return;
  }
  if(iLayer>=m_pLogicalHall->pShield->nLayerNum){
    G4cerr<<"Error: Layer number "<<iLayer<<" does not exist"<<G4endl;
    return;
  }
  m_pLogicalHall->pShield->pLayerMat[iLayer-1] = pMat;
  m_pLogicalHall->pShield->pLayers[iLayer-1]->SetMaterial(pMat);
  if(szMat=="Iron")
    vis.SetColour(192./255.,192./255.,192./255.);
  else
    vis.SetColour(92./255.,92./255.,92./255.);
  vis.SetVisibility(true);
  m_pLogicalHall->pShield->pLayers[iLayer]->SetVisAttributes(vis);
  Update();
}

void Hall::LayerDepth(G4double d,G4int iLayer)
{
  if(m_pHall->pShield->nLayerNum < iLayer){
    G4cerr<<"Error: Layer number "<<iLayer<<" does not exist"<<G4endl;
    return;
  }
  ClearImportance();
  m_pHall->pShield->pLayersThick[iLayer-1] = d*cm;
  
  //Refreshing layers
  G4double ztransl=0;
  G4double TotalThick=0;
  G4int i;
  for(i=0;i<m_pHall->pShield->nLayerNum;i++) TotalThick+=m_pHall->pShield->pLayersThick[i];
  ztransl = -TotalThick/2.;

  for(i=0;i<m_pPhysHall->pShield->nLayerNum;i++){
    if(i==iLayer-1){
      m_pHall->pShield->pLayers[i]->SetZHalfLength(d/2.*cm);

      m_pPhysHall->pShield->pLayers[i]->SetTranslation(G4ThreeVector(0,0,ztransl+m_pHall->pShield->pLayersThick[i]/2.));
    }
    else if(i>iLayer-1)
      m_pPhysHall->pShield->pLayers[i]->SetTranslation(G4ThreeVector(0,0,ztransl+m_pHall->pShield->pLayersThick[i]/2.));
    ztransl += m_pHall->pShield->pLayersThick[i];
  }
  m_pHall->pShield->pVLayer->SetZHalfLength(TotalThick/2.);
  m_pPhysHall->pShield->pVLayer->SetTranslation(G4ThreeVector(0,0,76.*cm+TotalThick/2.));
  m_pPhysHall->pDetector->SetTranslation(G4ThreeVector(0,0,82.35*cm+TotalThick));
  m_pPhysHall->pDetector20->SetTranslation(G4ThreeVector(20.*cm,0,82.35*cm+TotalThick));
  m_pPhysHall->pDetector40->SetTranslation(G4ThreeVector(40.*cm,0.,82.35*cm+TotalThick));
  ConstructParallelBase();
  ConstructParallelLayers();
  PostLayerInit();
  Update();
}
void Hall::CreateMaterials()
{
  ELEMENTTYPE Els;
  Els.pElementH = new G4Element("Hydro","H",1.,1.0079*g/mole);
  Els.pElementO = new G4Element("Oxy","O",8.,15.999*g/mole);
  Els.pElementN = new G4Element("Nitro","N",7.,14.006*g/mole);
  Els.pElementFe = new G4Element("Ferrum","Fe",26.,55.845*g/mole);
  Els.pElementNa = new G4Element("Natri","Na",11.,22.98977*g/mole);
  Els.pElementMg = new G4Element("Magn","Mg",12.,24.305*g/mole);
  Els.pElementAl = new G4Element("Alu","Al",13.,26.981538*g/mole);
  Els.pElementSi = new G4Element("Sili","Si",14.,28.0854*g/mole);
  Els.pElementK = new G4Element("Pota","K",19.,39.0983*g/mole);
  Els.pElementCa = new G4Element("Calc","Ca",20.,40.078*g/mole);
  Els.pElementLi7 = new G4Element("Lith-7","Li7",3.,7.*g/mole);
  Els.pElementLi = new G4Element("Lith","Li",3.,6.941*g/mole);

  MATERIALTYPE Mat;
  Mat.pMaterialFe = new G4Material("Iron",7.87*g/cm3,1);
  Mat.pMaterialFe->AddElement(Els.pElementFe,1);
  Mat.pMaterialAir = new G4Material("Air",1.29*mg/cm3,2);
  Mat.pMaterialAir->AddElement(Els.pElementN,0.75);
  Mat.pMaterialAir->AddElement(Els.pElementN,0.25);
  Mat.pMaterialTarget = new G4Material("Enrich Lithium",535*kg/m3,2);
  Mat.pMaterialTarget->AddElement(Els.pElementLi7,0.999);
  Mat.pMaterialTarget->AddElement(Els.pElementLi,0.001);
  Mat.pMaterialVacuum = new G4Material("Vacuum",universe_mean_density,1,
                                       kStateGas, 2.73*kelvin,3.e-18*pascal);
  Mat.pMaterialVacuum->AddElement(Els.pElementH,1);
  Mat.pMaterialConcrete = new G4Material("Concrete",2.31*g/cm3,9);
  Mat.pMaterialConcrete->AddElement(Els.pElementH,/*0.18957226*/0.010853422);
  Mat.pMaterialConcrete->AddElement(Els.pElementO, /*0.52999241*/0.48165594);
  Mat.pMaterialConcrete->AddElement(Els.pElementNa,/*0.01556568*/0.020327181);
  Mat.pMaterialConcrete->AddElement(Els.pElementMg,/*0.0078461149*/0.010832401);
  Mat.pMaterialConcrete->AddElement(Els.pElementAl,/*0.039483675*/0.060514396);
  Mat.pMaterialConcrete->AddElement(Els.pElementSi,/*0.14047077*/0.22409956);
  Mat.pMaterialConcrete->AddElement(Els.pElementK, /*0.0048089091*/0.010680188);
  Mat.pMaterialConcrete->AddElement(Els.pElementCa,/*0.054416603*/0.12388306);
  Mat.pMaterialConcrete->AddElement(Els.pElementFe,/*0.017843584*/0.056603178);
}

void Hall::CreateHall()
{
  if(!m_pHall){
    m_pHall = new HALLTYPE;
    m_pHall->pShield = NULL;
    m_pHall->dTargetThick = 3.6*mm;
  }
  if(!m_pLogicalHall){
    m_pLogicalHall = new LOGICALHALL;
    m_pLogicalHall->pShield = NULL;
  }
  if(!m_pPhysHall){
    m_pPhysHall = new PHYSICALHALL;
    m_pPhysHall->pShield = NULL;
  }

  static G4VisAttributes visMoth;
  static G4VisAttributes visConcrete;
  static G4VisAttributes visIron;

  static G4VisAttributes visTarget;
  static G4VisAttributes visDetectors;
  m_pHall->pMother = new G4Tubs("Mother Solid",0.,100.*cm,350.*cm,0,2*pi);
  m_pHall->pDeadConcreteA = new G4Tubs("Concrete Shield A Solid",26.*cm,100.*cm,110.*cm,0.,2.*pi);
  m_pHall->pDeadConcreteB = new G4Tubs("Concrete Shield B Solid",59.*cm,100.*cm,60.*cm,0.,2.*pi);
  m_pHall->pIronShieldA = new G4Tubs("Iron Shield A Solid",5.45*cm,26.*cm,110.*cm,0.,2.*pi);
  m_pHall->pIronShieldB = new G4Tubs("Iron Shield B Solid",5.45*cm,59.*cm,2.5*cm,0.,2.*pi);
  m_pHall->pTarget = new G4Tubs("Target Solid",0.,5.*cm,1.8*mm,0.,2.*pi);
  m_pHall->pDetector = new G4Tubs("Detector Solid",0,6.35*cm,6.35*cm,0.,2.*pi);
  m_pHall->pDetector20 = new G4Tubs("Detector Solid 20 cm",0.,6.35*cm,6.35*cm,0.,2.*pi);
  m_pHall->pDetector40 = new G4Tubs("Detector Solid 40 cm",0.,6.35*cm,6.35*cm,0.,2.*pi);
  m_pHall->pGhostDetector = new G4Tubs("Ghost Solid",0.,5.*cm,112.5*cm,0.,2.*pi);

  m_pLogicalHall->pMother = new G4LogicalVolume(m_pHall->pMother,
                                                G4Material::GetMaterial("Air"),
                                                "Mother Logical");
  visMoth.SetVisibility(true);
  visMoth.SetDaughtersInvisible(false);
  visMoth.SetColour(245./255.,245./255.,1.,0.01);

  m_pLogicalHall->pMother->SetVisAttributes(visMoth);

  m_pLogicalHall->pDeadConcreteA = new G4LogicalVolume(m_pHall->pDeadConcreteA,
                                                      G4Material::GetMaterial("Concrete"),
                                                      "Dead Concrete Logical A");
  m_pLogicalHall->pDeadConcreteB = new G4LogicalVolume(m_pHall->pDeadConcreteB,
                                                       G4Material::GetMaterial("Concrete"),
                                                       "Dead Concrete Logical B");
  visConcrete.SetColour(92./255.,92./255.,92./255.);
  visConcrete.SetVisibility(true);
  visConcrete.SetForceSolid(true);
  visConcrete.SetForceWireframe(false);
  m_pLogicalHall->pDeadConcreteA->SetVisAttributes(visConcrete);
  m_pLogicalHall->pDeadConcreteB->SetVisAttributes(visConcrete);

  m_pLogicalHall->pIronShieldA = new G4LogicalVolume(m_pHall->pIronShieldA,
                                                    G4Material::GetMaterial("Iron"),
                                                    "Iron Shield Logical A");
  m_pLogicalHall->pIronShieldB = new G4LogicalVolume(m_pHall->pIronShieldB,
                                                     G4Material::GetMaterial("Iron"),
                                                     "Iron Shield Logical B");
  visIron.SetColour(192./255.,192./255.,192./255.);
  visIron.SetVisibility(true);
  visIron.SetForceSolid(true);
  visIron.SetForceWireframe(false);
  m_pLogicalHall->pIronShieldA->SetVisAttributes(visIron);
  m_pLogicalHall->pIronShieldB->SetVisAttributes(visIron);

  m_pLogicalHall->pTarget = new G4LogicalVolume(m_pHall->pTarget,
                                                G4Material::GetMaterial("Enrich Lithium"),
                                                "Target Logical");
  visTarget.SetColour(1.,1.,1.);
  visTarget.SetVisibility(true);
  m_pLogicalHall->pTarget->SetVisAttributes(visTarget);

  m_pLogicalHall->pDetector = new G4LogicalVolume(m_pHall->pDetector,
                                                  G4Material::GetMaterial("Vacuum"),
                                                  "Detector Logical");
  m_pLogicalHall->pDetector20 = new G4LogicalVolume(m_pHall->pDetector20,
                                                    G4Material::GetMaterial("Vacuum"),
                                                    "Detector 20 cm Logical");
  m_pLogicalHall->pDetector40 = new G4LogicalVolume(m_pHall->pDetector40,
                                                    G4Material::GetMaterial("Vacuum"),
                                                    "Detector 40 cm Logical");
  visDetectors.SetColour(0.,0.,1.);
  visDetectors.SetVisibility(true);
  visDetectors.SetForceSolid(true);
  visDetectors.SetForceWireframe(false);
  m_pLogicalHall->pDetector->SetVisAttributes(visDetectors);
  m_pLogicalHall->pDetector20->SetVisAttributes(visDetectors);
  m_pLogicalHall->pDetector40->SetVisAttributes(visDetectors);

  m_pLogicalHall->pGhostDetector = new G4LogicalVolume(m_pHall->pGhostDetector,
                                                       G4Material::GetMaterial("Vacuum"),
                                                       "Ghost Logical");

  m_pLogicalHall->pGhostDetector->SetVisAttributes(visDetectors);


  m_pPhysHall->pMother = new G4PVPlacement(0,G4ThreeVector(),
                                           "Mother",
                                           m_pLogicalHall->pMother,
                                           NULL,
                                           false,
                                           0);
  m_pPhysHall->pDeadConcreteA = new G4PVPlacement(0,G4ThreeVector(0,0,-39.*cm),
                                                  "Dead Concrete A",
                                                  m_pLogicalHall->pDeadConcreteA,
                                                  m_pPhysHall->pMother,
                                                  false,
                                                  0);
  m_pPhysHall->pDeadConcreteB = new G4PVPlacement(0,G4ThreeVector(0,0,131*cm),
                                                  "Dead Concrete B",
                                                  m_pLogicalHall->pDeadConcreteB,
                                                  m_pPhysHall->pMother,
                                                  false,
                                                  0);
  m_pPhysHall->pIronShieldA = new G4PVPlacement(0,G4ThreeVector(0,0,-39*cm),
                                                "Iron Shield A",
                                                m_pLogicalHall->pIronShieldA,

                                                m_pPhysHall->pMother,
                                                false,
                                                0);
  m_pPhysHall->pIronShieldB = new G4PVPlacement(0,G4ThreeVector(0,0,73.5*cm),
                                                "Iron Shield B",
                                                m_pLogicalHall->pIronShieldB,
                                                m_pPhysHall->pMother,
                                                false,
                                                0);
  m_pPhysHall->pTarget = new G4PVPlacement(0,G4ThreeVector(0,0,-325.*cm-
                                                           m_pHall->dTargetThick/2.),

                                           "Target",
                                           m_pLogicalHall->pTarget,
                                           m_pPhysHall->pMother,
                                           false,
                                           0);
  m_pPhysHall->pDetector = new G4PVPlacement(0,G4ThreeVector(0,0,232.35*cm),
                                             "Detector",
                                             m_pLogicalHall->pDetector,
                                             m_pPhysHall->pMother,
                                             false,
                                             0);
  m_pPhysHall->pDetector20 = new G4PVPlacement(0,G4ThreeVector(20.*cm,0,232.35*cm),
                                               "Detector 20 cm",
                                               m_pLogicalHall->pDetector20,
                                               m_pPhysHall->pMother,
                                               false,
                                               0);
  m_pPhysHall->pDetector40 = new G4PVPlacement(0,G4ThreeVector(40.*cm,0,232.35*cm),
                                               "Detector 40 cm",
                                               m_pLogicalHall->pDetector40,
                                               m_pPhysHall->pMother,
                                               false,
                                               0);
  m_pPhysHall->pGhostDetector = new G4PVPlacement(0,G4ThreeVector(0,0,-36.5*cm),
                                                  "Ghost",
                                                  m_pLogicalHall->pGhostDetector,
                                                  m_pPhysHall->pMother,
                                                  false,
                                                  0);

  //Ostatyca se otnasja za parallelnata geometrija
  ConstructParallelBase();
}
static void ConstructParallelBase()
{
  pPWorldSolid = new G4Tubs("Parallel world solid",0,110*cm,360*cm,0,2*pi);
  pPWorldLog = new G4LogicalVolume(pPWorldSolid,G4Material::GetMaterial("Air"),"Parallel world logical");
  pPWorld = new G4PVPlacement(0,G4ThreeVector(),"Parallel World",pPWorldLog,NULL,false,0);
  pPGeom = new G4IStore(*pPWorld);
  pPGeom->AddImportanceRegion(1,*pPWorld,-1);
}

void Hall::CreateTarget(G4double Thick)
{
  if(m_pHall->pTarget)
    delete m_pHall->pTarget;
  m_pHall->pTarget = new G4Tubs("Target",0,5.*cm,Thick/2.,0.,2.*pi);
  m_pLogicalHall->pTarget->SetSolid(m_pHall->pTarget);
  m_pPhysHall->pTarget->SetTranslation(G4ThreeVector(0.,0.,-325.*cm-m_pHall->dTargetThick/2.));
  G4RunManager::GetRunManager()->DefineWorldVolume(m_pPhysHall->pMother);
}
void Hall::Update()
{
  DumpLayers();
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->DefineWorldVolume(m_pPhysHall->pMother);
}

G4VPhysicalVolume* Hall::Construct()
{
  CreateMaterials();
  CreateHall();
  return m_pPhysHall->pMother;
}

void Hall::DumpLayers()
{
#ifdef G4DEBUG
  G4int i,j;
  if(!m_pHall->pShield) return;
  G4Tubs** ppSolids = m_pHall->pShield->pLayers;
  G4VPhysicalVolume* pLayer;
  j = m_pHall->pShield->nLayerNum;
  G4cout<<"Number of solids: "<<j<<G4endl;
  for(i=0;i<j;i++){
    G4cout<<"Layer solid: "<<ppSolids[i]->GetName()<<" Thick: "<<m_pHall->pShield->pLayersThick[i]<<G4endl;
    pLayer = m_pPhysHall->pShield->pLayers[i];
    G4cout<<"Layer: "<<pLayer->GetName()<<" at: "<<pLayer->GetTranslation()<<G4endl;
  }
#endif
}
void Hall::ResetGeometry()
{
  G4int i,j = m_pHall->pShield->nLayerNum;
  for(i=0;i<j;i++){
    m_pLogicalHall->pShield->pVLayer->RemoveDaughter(m_pPhysHall->pShield->pLayers[i]);
    delete m_pHall->pShield->pLayers[i];
    delete m_pLogicalHall->pShield->pLayers[i];
    delete m_pPhysHall->pShield->pLayers[i];
  }
  m_pLogicalHall->pMother->RemoveDaughter(m_pPhysHall->pShield->pVLayer);
  delete m_pHall->pShield->pLayers;
  delete m_pHall->pShield->pLayersThick;
  delete m_pHall->pShield->pLayersColimator;
  delete m_pHall->pShield->pVLayer;
  delete m_pHall->pShield;
  m_pHall->pShield = NULL;

  delete m_pLogicalHall->pShield->pLayerMat;
  delete m_pLogicalHall->pShield->pLayers;
  delete m_pLogicalHall->pShield->pVLayer;
  delete m_pLogicalHall->pShield;
  m_pLogicalHall->pShield = NULL;

  delete m_pPhysHall->pShield->pLayers;
  delete m_pPhysHall->pShield->pVLayer;
  delete m_pPhysHall->pShield;
  ClearImportance();
  m_pPhysHall->pShield = NULL;
  m_pPhysHall->pDetector->SetTranslation(G4ThreeVector(0,0,232.35*cm));
  m_pPhysHall->pDetector20->SetTranslation(G4ThreeVector(20.*cm,0,232.35*cm));
  m_pPhysHall->pDetector40->SetTranslation(G4ThreeVector(40.*cm,0,232.35*cm));
  Update();
}

G4double Hall::GetColimPos()
{
  G4double dRes = 401*cm;
  if(m_pHall->pShield){
    for(G4int j=0;j<m_pHall->pShield->nLayerNum;j++){
      if(!m_pHall->pShield->pLayersColimator[j]) break;
      dRes += m_pHall->pShield->pLayersThick[j];
    }
  }
  return dRes;
}

void Hall::DumpImportance(G4String fileName)
{
  /* M
  std::ofstream of(fileName);
  if(pNeutronSampler){
    of << "-- results for importanct sampled particle: ---" 
       <<  "neutron" << G4endl;
    G4ScoreTable sp(pPGeom);
    sp.PrintHeader(&of);
    sp.PrintTable(pPartScorer[0]->GetMapPtkTallys(),  &of);
    of << G4endl;
    of << "----------------------------------------------- " << G4endl;
  }
  for(unsigned i=1;i<6;i++){
    if(pRestSamplers[i]){
      of << "-- results for scored particle: ---------------" 
         << pArrNames[i] << G4endl;
      G4ScorePrinter sp;
      sp.PrintHeader(&of);
      sp.PrintTable(pPartScorer[i]->GetMapPtkTallys(),  &of);
      of << G4endl;
      of << "----------------------------------------------- " << G4endl;

    }
  }
  */ 
}

void Hall::ResetParallelGeometry()
{
  DeleteParallelGeometry();
  InitializeParallelGeometry();
}

void Hall::UseImportance()
{
  if(bInit){
    G4cout<<"Importance processes are already in use"<<G4endl;
    return;

  }
  InitializeParallelGeometry();
}

void Hall::DontUseImportance()
{
  if(!bInit){
    G4cout<<"Importance processes are not in use"<<G4endl;
    return;
  }
  DeleteParallelGeometry();
}

void Hall::SetBlocksPerLayer(G4int nBl)
{
  if(nBl<1) return;

  if(nBl!=nParallelLayers){
    ClearImportance();
    nParallelLayers=nBl;
    ConstructParallelBase();
    ConstructParallelLayers();
    PostLayerInit();
  }
}
