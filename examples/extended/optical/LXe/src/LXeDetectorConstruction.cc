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
#include "LXeDetectorConstruction.hh"
#include "LXePMTSD.hh"
#include "LXeScintSD.hh"
#include "LXeDetectorMessenger.hh"
#include "LXeMainVolume.hh"
#include "LXeWLSSlab.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryManager.hh"
#include "G4UImanager.hh"

G4bool LXeDetectorConstruction::sphereOn = true;

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeDetectorConstruction::LXeDetectorConstruction()
: LXe_mt(NULL), MPTPStyrene(NULL)
{
  SetDefaults();
  detectorMessenger = new LXeDetectorMessenger(this);
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeDetectorConstruction::~LXeDetectorConstruction(){
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeDetectorConstruction::DefineMaterials(){
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;

  G4int polyPMMA = 1;
  G4int nC_PMMA = 3+2*polyPMMA;
  G4int nH_PMMA = 6+2*polyPMMA;

  G4int polyeth = 1;
  G4int nC_eth = 2*polyeth;
  G4int nH_eth = 4*polyeth;

  //***Elements
  H = new G4Element("H", "H", z=1., a=1.01*g/mole);
  C = new G4Element("C", "C", z=6., a=12.01*g/mole);
  N = new G4Element("N", "N", z=7., a= 14.01*g/mole);
  O = new G4Element("O"  , "O", z=8., a= 16.00*g/mole);
  
  //***Materials
  //Liquid Xenon
  LXe = new G4Material("LXe",z=54.,a=131.29*g/mole,density=3.020*g/cm3);
  //Aluminum
  Al = new G4Material("Al",z=13.,a=26.98*g/mole,density=2.7*g/cm3);
  //Vacuum
  Vacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole,
			  density=universe_mean_density,kStateGas,0.1*kelvin,
			  1.e-19*pascal); 
  //Air
  Air = new G4Material("Air", density= 1.29*mg/cm3, 2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);
  //Glass
  Glass = new G4Material("Glass", density=1.032*g/cm3,2);
  Glass->AddElement(C,91.533*perCent);
  Glass->AddElement(H,8.467*perCent);
  //Polystyrene
  Pstyrene = new G4Material("Polystyrene", density= 1.03*g/cm3, 2);
  Pstyrene->AddElement(C, 8);
  Pstyrene->AddElement(H, 8);
  //Fiber(PMMA)
  PMMA = new G4Material("PMMA", density=1190*kg/m3,3);
  PMMA->AddElement(H,nH_PMMA);
  PMMA->AddElement(C,nC_PMMA);
  PMMA->AddElement(O,2);
  //Cladding(polyethylene)
  Pethylene = new G4Material("Pethylene", density=1200*kg/m3,2);
  Pethylene->AddElement(H,nH_eth);
  Pethylene->AddElement(C,nC_eth);
  //Double cladding(flourinated polyethylene)
  fPethylene = new G4Material("fPethylene", density=1400*kg/m3,2);
  fPethylene->AddElement(H,nH_eth);
  fPethylene->AddElement(C,nC_eth);
  
  //***Material properties tables

  const G4int LXe_NUMENTRIES = 3;
  G4double LXe_Energy[LXe_NUMENTRIES]    = { 7.0*eV , 7.07*eV, 7.14*eV };

  G4double LXe_SCINT[LXe_NUMENTRIES] = { 0.1, 1.0, 0.1 };
  G4double LXe_RIND[LXe_NUMENTRIES]  = { 1.59 , 1.57, 1.54 };
  G4double LXe_ABSL[LXe_NUMENTRIES]  = { 35.*cm, 35.*cm, 35.*cm}; 
  LXe_mt = new G4MaterialPropertiesTable();
  LXe_mt->AddProperty("FASTCOMPONENT", LXe_Energy, LXe_SCINT, LXe_NUMENTRIES);
  LXe_mt->AddProperty("SLOWCOMPONENT", LXe_Energy, LXe_SCINT, LXe_NUMENTRIES);
  LXe_mt->AddProperty("RINDEX",        LXe_Energy, LXe_RIND,  LXe_NUMENTRIES);
  LXe_mt->AddProperty("ABSLENGTH",     LXe_Energy, LXe_ABSL,  LXe_NUMENTRIES);
  LXe_mt->AddConstProperty("SCINTILLATIONYIELD",12000./MeV); 
  LXe_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
  LXe_mt->AddConstProperty("FASTTIMECONSTANT",20.*ns);
  LXe_mt->AddConstProperty("SLOWTIMECONSTANT",45.*ns);
  LXe_mt->AddConstProperty("YIELDRATIO",1.0);
  LXe->SetMaterialPropertiesTable(LXe_mt);

  // Set the Birks Constant for the LXe scintillator

  LXe->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
  
  G4double Glass_RIND[LXe_NUMENTRIES]={1.49,1.49,1.49};
  G4double Glass_AbsLength[LXe_NUMENTRIES]={420.*cm,420.*cm,420.*cm};
  G4MaterialPropertiesTable *Glass_mt = new G4MaterialPropertiesTable();
  Glass_mt->AddProperty("ABSLENGTH",LXe_Energy,Glass_AbsLength,LXe_NUMENTRIES);
  Glass_mt->AddProperty("RINDEX",LXe_Energy,Glass_RIND,LXe_NUMENTRIES);
  Glass->SetMaterialPropertiesTable(Glass_mt);

  G4double Vacuum_Energy[LXe_NUMENTRIES]={2.0*eV,7.0*eV,7.14*eV};
  G4double Vacuum_RIND[LXe_NUMENTRIES]={1.,1.,1.};  
  G4MaterialPropertiesTable *Vacuum_mt = new G4MaterialPropertiesTable();
  Vacuum_mt->AddProperty("RINDEX", Vacuum_Energy, Vacuum_RIND,LXe_NUMENTRIES);
  Vacuum->SetMaterialPropertiesTable(Vacuum_mt);
  Air->SetMaterialPropertiesTable(Vacuum_mt);//Give air the same rindex

  const G4int WLS_NUMENTRIES = 4;
  G4double WLS_Energy[] = {2.00*eV,2.87*eV,2.90*eV,3.47*eV};
    
  G4double RIndexPstyrene[WLS_NUMENTRIES]={ 1.5, 1.5, 1.5, 1.5};
  G4double Absorption1[WLS_NUMENTRIES]={2.*cm, 2.*cm, 2.*cm, 2.*cm};
  G4double ScintilFast[WLS_NUMENTRIES]={0.00, 0.00, 1.00, 1.00};
  MPTPStyrene = new G4MaterialPropertiesTable();
  MPTPStyrene->AddProperty("RINDEX",WLS_Energy,RIndexPstyrene,WLS_NUMENTRIES);
  MPTPStyrene->AddProperty("ABSLENGTH",WLS_Energy,Absorption1,WLS_NUMENTRIES);
  MPTPStyrene->AddProperty("FASTCOMPONENT",WLS_Energy, ScintilFast,
			   WLS_NUMENTRIES);
  MPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);
  MPTPStyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPTPStyrene->AddConstProperty("FASTTIMECONSTANT", 10.*ns);
  Pstyrene->SetMaterialPropertiesTable(MPTPStyrene);

  // Set the Birks Constant for the Polystyrene scintillator

  Pstyrene->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  G4double RefractiveIndexFiber[WLS_NUMENTRIES]={ 1.60, 1.60, 1.60, 1.60};
  G4double AbsFiber[WLS_NUMENTRIES]={9.00*m,9.00*m,0.1*mm,0.1*mm};
  G4double EmissionFib[WLS_NUMENTRIES]={1.0, 1.0, 0.0, 0.0};
  G4MaterialPropertiesTable* MPTFiber = new G4MaterialPropertiesTable();
  MPTFiber->AddProperty("RINDEX",WLS_Energy,RefractiveIndexFiber,
			WLS_NUMENTRIES);
  MPTFiber->AddProperty("WLSABSLENGTH",WLS_Energy,AbsFiber,WLS_NUMENTRIES);
  MPTFiber->AddProperty("WLSCOMPONENT",WLS_Energy,EmissionFib,WLS_NUMENTRIES);
  MPTFiber->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
  PMMA->SetMaterialPropertiesTable(MPTFiber);

  G4double RefractiveIndexClad1[WLS_NUMENTRIES]={ 1.49, 1.49, 1.49, 1.49};
  G4MaterialPropertiesTable* MPTClad1 = new G4MaterialPropertiesTable();
  MPTClad1->AddProperty("RINDEX",WLS_Energy,RefractiveIndexClad1,
			WLS_NUMENTRIES);
  MPTClad1->AddProperty("ABSLENGTH",WLS_Energy,AbsFiber,WLS_NUMENTRIES);
  Pethylene->SetMaterialPropertiesTable(MPTClad1);

  G4double RefractiveIndexClad2[WLS_NUMENTRIES]={ 1.42, 1.42, 1.42, 1.42};
  G4MaterialPropertiesTable* MPTClad2 = new G4MaterialPropertiesTable();
  MPTClad2->AddProperty("RINDEX",WLS_Energy,RefractiveIndexClad2,
			WLS_NUMENTRIES);
  MPTClad2->AddProperty("ABSLENGTH",WLS_Energy,AbsFiber,WLS_NUMENTRIES);
  fPethylene->SetMaterialPropertiesTable(MPTClad2);
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4VPhysicalVolume* LXeDetectorConstruction::Construct(){
  DefineMaterials();
  return ConstructDetector();
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4VPhysicalVolume* LXeDetectorConstruction::ConstructDetector()
{
  //The experimental hall walls are all 1m away from housing walls
  G4double expHall_x = scint_x+d_mtl+1.*m;
  G4double expHall_y = scint_y+d_mtl+1.*m;
  G4double expHall_z = scint_z+d_mtl+1.*m;

  //Create experimental hall
  experimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,
                                             Vacuum,"expHall_log",0,0,0);
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
			      experimentalHall_log,"expHall",0,false,0);

  experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  
  //Place the main volume
  if(mainVolume){
    new LXeMainVolume(0,G4ThreeVector(),experimentalHall_log,false,0,this);
  }

  //Place the WLS slab
  if(WLSslab){
    G4VPhysicalVolume* slab = new LXeWLSSlab(0,G4ThreeVector(0.,0.,
					     -scint_z/2.-slab_z-1.*cm),
					     experimentalHall_log,false,0,
					     this);

    //Surface properties for the WLS slab
    G4OpticalSurface* ScintWrap = new G4OpticalSurface("ScintWrap");
    
    new G4LogicalBorderSurface("ScintWrap", slab,
			       experimentalHall_phys,
			       ScintWrap);
    
    ScintWrap->SetType(dielectric_metal);
    ScintWrap->SetFinish(polished);
    ScintWrap->SetModel(glisur);

    const G4int NUM = 2;
    
    G4double pp[NUM] = {2.0*eV, 3.5*eV};
    G4double reflectivity[NUM] = {1., 1.};
    G4double efficiency[NUM] = {0.0, 0.0};
    
    G4MaterialPropertiesTable* ScintWrapProperty 
      = new G4MaterialPropertiesTable();

    ScintWrapProperty->AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
    ScintWrapProperty->AddProperty("EFFICIENCY",pp,efficiency,NUM);
    ScintWrap->SetMaterialPropertiesTable(ScintWrapProperty);
  }

  return experimentalHall_phys;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeDetectorConstruction::SetDimensions(G4ThreeVector dims){
  this->scint_x=dims[0];
  this->scint_y=dims[1];
  this->scint_z=dims[2];
  updated=true;
}
 
 //_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeDetectorConstruction::SetHousingThickness(G4double d_mtl){
  this->d_mtl=d_mtl;
  updated=true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeDetectorConstruction::SetNX(G4int nx){
  this->nx=nx;
  updated=true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeDetectorConstruction::SetNY(G4int ny){
  this->ny=ny;
  updated=true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeDetectorConstruction::SetNZ(G4int nz){
  this->nz=nz;
  updated=true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeDetectorConstruction::SetPMTRadius(G4double outerRadius_pmt){
  this->outerRadius_pmt=outerRadius_pmt;
  updated=true;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeDetectorConstruction::SetDefaults(){
  //Resets to default values
  d_mtl=0.0635*cm;
  
  scint_x = 17.8*cm;
  scint_y = 17.8*cm;
  scint_z = 22.6*cm;

  nx = 2;
  ny = 2;
  nz = 3;

  outerRadius_pmt = 2.3*cm;

  sphereOn = true;
  refl=1.0;
  
  nfibers=15;
  WLSslab=false;
  mainVolume=true;
  slab_z=2.5*mm;

  G4UImanager::GetUIpointer()
    ->ApplyCommand("/LXe/detector/scintYieldFactor 1.");
  
  if(LXe_mt)LXe_mt->AddConstProperty("SCINTILLATIONYIELD",12000./MeV);
  if(MPTPStyrene)MPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);

  updated=true;
}


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeDetectorConstruction::UpdateGeometry(){

  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();
  G4SurfaceProperty::CleanSurfacePropertyTable();

  //define new one
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  updated=false;
}

void LXeDetectorConstruction::SetMainScintYield(G4double y){
  LXe_mt->AddConstProperty("SCINTILLATIONYIELD",y/MeV);  
}
 
void LXeDetectorConstruction::SetWLSScintYield(G4double y){
  MPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",y/MeV); 
}





