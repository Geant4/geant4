////////////////////////////////////////////////////////////////////////////////
//
#include "MLGeometryConstruction.hh"
#include "MLGeometryMessenger.hh"
#include "MLRunManager.hh"
#include "MLAnalysisManager.hh"
#include "MLFluenceAnalyser.hh"
#include "MLDoseAnalyser.hh"
#include "MLNIELAnalyser.hh"

#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include <strstream>
////////////////////////////////////////////////////////////////////////////////
//
MLGeometryConstruction::MLGeometryConstruction ()
  :shape(SLAB),defaultMaterial(NULL),solidWorld(NULL), solidSWorld(NULL),
   logicWorld(NULL),physiWorld(NULL),materialsManager(NULL),detectorSD(NULL)
{
  SetToDefault();
  ComputeParameters();
  materialsManager = new MLMaterial () ;
  colourManager = new MLColour ();

  // create commands for interactive definition of the calorimeter  
  geometryMessenger = new MLGeometryMessenger(this);
}
////////////////////////////////////////////////////////////////////////////////
//
MLGeometryConstruction::~MLGeometryConstruction ()
{
  delete geometryMessenger;
  delete materialsManager;
  delete colourManager;
  
//  delete detectorSD;
}
////////////////////////////////////////////////////////////////////////////////
//
G4VPhysicalVolume* MLGeometryConstruction::Construct ()
{
  return ConstructGeometry();
}
////////////////////////////////////////////////////////////////////////////////
//
void MLGeometryConstruction::DefineMaterials ()
{ 
//  materialsManager = new MLMaterial();
}
////////////////////////////////////////////////////////////////////////////////
//
G4VPhysicalVolume* MLGeometryConstruction::ConstructGeometry ()
{
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  ComputeParameters();

//
//
// Define volumes defining the World, either in the form of a box or sphere.
//
  defaultMaterial = materialsManager->GetMaterial(0); // Vacuum

  if (shape == SLAB) {
    solidWorld = new G4Box("World",
      WorldSizeXY/2., WorldSizeXY/2., WorldSizeZ/2.*1.01);
    logicWorld = new G4LogicalVolume(solidWorld, defaultMaterial,
      "World");
  } else {
    solidSWorld = new G4Sphere("World", 0, WorldSizeZ*1.01,
      0.*deg, 360.*deg, 0.*deg, 180.*deg);
    logicWorld = new G4LogicalVolume(solidSWorld, defaultMaterial,
      "World");
  }

  physiWorld = new G4PVPlacement(0,                     //no rotation
                                 G4ThreeVector(),       //at (0,0,0)
                                 "World",               //its name
                                 logicWorld,            //its logical volume
                                 NULL,                  //its mother  volume
                                 false,                 //no boolean operation
                                 0);                    //copy number

  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
//
//
// Now define the shielding stack - either slab or spherical.
//
  G4String tname = "Layer-" ;
  G4String lname, tnum;
  char x[10];

  logicLayer.clear();
  physiLayer.clear();
  G4VisAttributes* layer_vat;

  if (shape == SLAB) {
//
//
// Slab shield case.
//
    solidLayer.clear();
    for (G4int k=0; k<NbOfLayers; k++) {
      std::ostrstream os(x,10);
      os <<k+1 <<'\0';
      tnum  = x;
      lname = tname + tnum;

      solidLayer.push_back(
        new G4Box(lname, WorldSizeXY/2., WorldSizeXY/2., Layers[k].x()/2.)) ;
      logicLayer.push_back(new G4LogicalVolume(solidLayer[k],
        materialsManager->GetMaterial(G4int(Layers[k].y())), lname));

      G4double zcentre = 0.5*Layers[k].x() - WorldSizeZ/2.;
      for (G4int l = 0; l < k; l++) zcentre += Layers[l].x();

      physiLayer.push_back(new G4PVPlacement(0,		//no rotation
                           G4ThreeVector(0.,0.,zcentre),//its position
                           lname,                       //its name
                           logicLayer[k],               //its logical volume
                           physiWorld,                  //its mother
                           false,                       //no boulean operat
                           k) );                        //copy number
      
      layer_vat = new G4VisAttributes(
        colourManager->GetColour(G4int(Layers[k].z())));
      layer_vat->SetVisibility(true);
      logicLayer[k]->SetVisAttributes(layer_vat);
    }
  } else {
//
//
// Spherical shield case.
//
    solidSphere.clear();
    G4double r1 = 0.;
    G4double r2 = WorldSizeZ;
    for (G4int k=0; k<NbOfLayers; k++) {
      std::ostrstream os(x,10);
      os << k+1 << '\0';
      tnum  = x;
      lname = tname + tnum;
      r1    = r2 - Layers[k].x();
      if (r1 < 1e-10) r1 = 0.;

      solidSphere.push_back (new G4Sphere(lname, r1, r2, 0.*deg, 360.*deg,
        0.*deg, 180.*deg));
      logicLayer.push_back (new G4LogicalVolume(solidSphere[k],
        materialsManager->GetMaterial(G4int(Layers[k].y())), lname) );
      physiLayer.push_back (new G4PVPlacement(0,             //no rotation
                           G4ThreeVector(0.,0.,0.),          //its position
                           lname,                            //its name
                           logicLayer[k],                    //its logical volume
                           physiWorld,                       //its mother
                           false,                            //no boulean operat
                           k));                              //copy number
      
      layer_vat = new G4VisAttributes(
        colourManager->GetColour(G4int(Layers[k].z())));
      layer_vat->SetVisibility(true);
      logicLayer[k]->SetVisAttributes(layer_vat);
      r2 = r1;
    }
  }

  //
  // Energy Sensitive Detectors: defaul tall Layers
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if (!detectorSD) {
    detectorSD = new MLSD("detectorSD",this);
    SDman->AddNewDetector( detectorSD );
  }
//
//
// Set sensitive detector for layers to be analysed for PHS.
//
  G4int k = 0;
  for (size_t l=0; l<EdepLayers.size(); l++) {
    k = EdepLayers[l];
    logicLayer[k]->SetSensitiveDetector(detectorSD);
  }
//
//
// Now do the same for layers to be analysed for dose.
//
  MLAnalysisManager* analysisManager = MLAnalysisManager::getInstance();
  MLDoseAnalyser* doseAnalyser       = analysisManager->GetDoseAnalyser();
  for (G4int i = 0; i < doseAnalyser->GetNbOfDLayers(); i++) {
    k = doseAnalyser->GetDLayerIdx(i);
    if (logicLayer[k]->GetSensitiveDetector() == NULL)
      logicLayer[k]->SetSensitiveDetector(detectorSD);
  }
  return physiWorld;
}
////////////////////////////////////////////////////////////////////////////////
//
const G4LogicalVolume* MLGeometryConstruction::GetLogicalLayer (G4int l)
{
  return logicLayer[l];
}
////////////////////////////////////////////////////////////////////////////////
//
void MLGeometryConstruction::AddLayer (G4int ival, G4String mat, G4int kval,
  G4double tk)
{
  G4int j = materialsManager->GetMaterialIndex(mat);
  G4int k = kval-1;
  if (ival < 0 || ival > NbOfLayers) {
    G4cerr <<G4endl;
    G4cerr <<"AddLayer : position " <<ival
           <<" must be >= 0 and <= " <<NbOfLayers <<"." <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
    return;
  }

  if (j < 0 || j > materialsManager->GetNbOfMaterial()) {
    G4cerr <<G4endl;
    G4cerr <<"AddLayer : material : " <<mat <<" is not defined." <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
    return;
  }

  if (k < 0 || k >= colourManager->GetNbOfColours()) {
    G4cerr <<G4endl;
    G4cerr <<"AddLayer : colour index " <<k+1 <<" must be >0 and <"
           <<colourManager->GetNbOfColours() <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
    return;
  }

  G4ThreeVector layer(tk, G4double(j), G4double(k) );
  if (NbOfLayers == 0) {
    Layers.push_back(layer);
    NbOfLayers++;
    return;
  }
  else {
    Layers.insert(Layers.begin()+ival,layer);
    NbOfLayers++;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
//
void MLGeometryConstruction::DeleteLayer(G4int ival)
{
  // search the material by its name
  //
  if (ival > NbOfLayers || ival < 0) {
    G4cerr <<G4endl;
    G4cerr <<"DeleteLayer: layer number " <<ival <<" out of range." <<G4endl;
    G4cerr <<"-->  Command rejected." <<G4endl;
  } else if (ival == 0 ) {
    Layers.clear();
    NbOfLayers = 0;
    EdepLayers.clear();
    //    EdepLayers.push_back(0);
    FluxLayers.clear();
    //    FluxLayers.push_back(0);
    // need to take care of the particles seletcted as well 
    //    MLAnalysisManager::getInstance()->GetFluenceAnalyser()->DeleteSPart ("proton");
    // MLAnalysisManager::getInstance()->GetFluenceAnalyser()->DeleteSPart ("neutron");
  } else {
    G4int i = ival-1;
    Layers.erase (Layers.begin()+i);
    NbOfLayers--;
//
//
// Check whether the layer corresponds to an energy-deposition or flux layer.
//
    size_t j = 0;
    while (j < EdepLayers.size()) {
      if (EdepLayers[j] == i) {
        EdepLayers.erase(EdepLayers.begin()+j);
      } else if (EdepLayers[j] > i) {
        EdepLayers[j]--;
        j++;
      } else {
        j++;
      }
    }
    j = 0;
    while (j < FluxLayers.size()) {
      if (FluxLayers[j] == i) {
        FluxLayers.erase(FluxLayers.begin()+j);
      } else if (FluxLayers[j] > i) {
        FluxLayers[j]--;
        j++;
      } else {
        j++;
      }
    }

    MLDoseAnalyser* doseAnalyser = MLAnalysisManager::getInstance()
      ->GetDoseAnalyser();
    doseAnalyser->DeleteLayer(i);
    MLNIELAnalyser* nielAnalyser = MLAnalysisManager::getInstance()
      ->GetNIELAnalyser();
    nielAnalyser->DeleteLayer(i);
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLGeometryConstruction::ListLayer (G4int l)
{
  G4cout <<G4endl;
  G4cout <<"-------------------------------------------------------------"
         <<G4endl;
  G4cout <<"---> The geometry setup has " <<NbOfLayers <<" layers:"
         <<G4endl;

  if ( l == 0 ) {
    for (G4int i=0; i<NbOfLayers; i++) {
      G4cout <<"  Layer No. " <<std::setw(3) <<i+1 <<" " <<std::setw(20)
             << materialsManager->GetMaterial(G4int(Layers[i].y()))->GetName()
             <<": " <<std::setw(6) <<G4BestUnit(Layers[i].x(),"Length")
             <<G4endl;
    }
  } else {
    G4cout <<"  Layer No. " <<std::setw(3) <<l <<" " <<std::setw(12)
           <<materialsManager->GetMaterial(G4int(Layers[l-1].y()))->GetName()
           <<": " <<std::setw(6) <<G4BestUnit(Layers[l-1].x(),"Length")
           <<G4endl;
  }
  G4cout <<G4endl;
  G4cout <<"-------------------------------------------------------------"
         <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
//
void MLGeometryConstruction::SetToDefault ()
{
  static const G4double thick[26] = {
    0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5,
    2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
  Layers.clear();
  shape      = SLAB;
  NbOfLayers = 26;
  G4ThreeVector layer(0., 0., 0.);
  G4double kc = 0.;
  layer.setX(thick[0]*mm);
  layer.setY(2.);
  layer.setZ(kc);
  Layers.push_back (layer); 
  for (G4int i=1; i<NbOfLayers; i++) {
    layer.setX((thick[i] - thick[i-1])*mm);
    kc++;
    if (kc > 11.) kc =0.;
    layer.setZ(kc); 
    Layers.push_back (layer); 
  }
  EdepLayers.clear();
  EdepLayers.push_back(0);
  FluxLayers.clear();
  FluxLayers.push_back(0);
  // why these below don't work?
  //  MLAnalysisManager* analysisManager = MLAnalysisManager::getInstance();
  // MLFluenceAnalyser* fluenceAnalyser       = analysisManager->GetFluenceAnalyser();
  // fluenceAnalyser->AddSPart ("proton");
  // fluenceAnalyser->AddSPart ("neutron");
}
////////////////////////////////////////////////////////////////////////////////
//
void MLGeometryConstruction::UpdateGeometry ()
{
  //  static G4bool geom = false ;
  MLRunManager* runManager =  static_cast <MLRunManager*>
                              (MLRunManager::GetRunManager());
  // if (!geom) {
    //    runManager->Initialize();
  //   geom = true ;
  // } else {
    runManager->DefineWorldVolume(ConstructGeometry());
    //    runManager->CutOffHasBeenModified();
    //    runManager->Initialize();
    //  }
  // Get the pointer to the UI manager
  G4UImanager* UI = G4UImanager::GetUIpointer();
  // UI->ApplyCommand("/control/execute display.mac");
  G4double xpos = -WorldSizeZ/2.*1.001;
  if (shape == SPHERE) xpos *= 2.;
  G4String command = "/gps/centre ";
  char x[50];
  std::ostrstream os(x,50);
  os <<xpos <<'\0';
  G4String xs = x;
  if (shape == SPHERE) {
    UI->ApplyCommand(command+"0. "+xs+" 0. mm");
    UI->ApplyCommand("/gps/angrot1 1. 0. 0.");
    UI->ApplyCommand("/gps/angrot2 0. 0. 1.");
    UI->ApplyCommand("/gps/direction 0 1 0");
  } else {
    UI->ApplyCommand(command+"0. 0. "+xs+" mm");
    UI->ApplyCommand("/gps/direction 0 0 1");
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLGeometryConstruction::AddEdepLayer (G4int ival)
{
  // set the number of Absorbers
  //
  if (ival < 0 || ival > NbOfLayers) {
    G4cerr <<G4endl;
    G4cerr <<"AddEdepLayer : Layer number " <<ival
           <<" must be >= 0 and <= " <<NbOfLayers <<"." <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else if (ival == 0) {
    EdepLayers.clear();
    for ( G4int i = 0; i < NbOfLayers; i++) EdepLayers.push_back(i);
  } else if (std::find(EdepLayers.begin(),EdepLayers.end(),ival-1)
    != EdepLayers.end()) {
    G4cerr <<G4endl;
    G4cerr <<"AddEdepLayer : Layer number " <<ival << " already selected."
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else {
    EdepLayers.push_back(ival-1);

  }
}
////////////////////////////////////////////////////////////////////////////////
//
//
void MLGeometryConstruction::DeleteEdepLayer (G4int ival)
{
  // search the material by its name
  //
  if (ival < 0 || ival > NbOfLayers) {
    G4cerr <<G4endl;
    G4cerr <<"DeleteEdepLayer: layer number " <<ival
           <<" out of range." <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else  if ( ival == 0) {
    EdepLayers.clear();

  } else if (std::find(EdepLayers.begin(),EdepLayers.end(),ival-1)
    == EdepLayers.end()) {
    G4cerr <<G4endl;
    G4cerr <<"DeleteEdepLayer : " <<ival <<" layer has not been selected!"
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else {
    EdepLayers.erase(std::find(EdepLayers.begin(),EdepLayers.end(),ival-1));

  }
}
////////////////////////////////////////////////////////////////////////////////
//
//
void MLGeometryConstruction::ListEdepLayer ()
{
  G4cout <<G4endl;
  G4cout <<"-------------------------------------------------------------"
         <<G4endl;
  G4cout <<"PHS for energy deposition are recorded for the following layers: "
         <<G4endl;
  if (EdepLayers.size() == 0) {
    G4cout <<G4endl;
    G4cout <<"No layer has been selected for PHS energy deposition analysis."
           <<G4endl;

  } else {
    for (size_t k=0; k<EdepLayers.size(); k++) {
      G4int i = EdepLayers[k];
      G4cout <<"  Layer No. " <<std::setw(3) <<i+1 <<" " <<std::setw(20)
             <<materialsManager->GetMaterial(G4int(Layers[i].y()))->GetName()
             <<": " <<std::setw(6) <<G4BestUnit(Layers[i].x(),"Length")
             <<G4endl;
    }
  }
  G4cout <<G4endl;
  G4cout <<"-------------------------------------------------------------"
         <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
//
void MLGeometryConstruction::AddFluxLayer (G4int ival)
{
  // set the number of Absorbers
  //
  if (ival < 0 || (shape == SLAB && ival > NbOfLayers+1) ||
                  (shape == SPHERE && ival > NbOfLayers)) {
    G4cerr <<G4endl;
    G4cerr <<"AddFluxLayer : Layer index " <<ival;
    if (shape == SLAB) {
      G4cerr <<" must be >0 and and <= " <<NbOfLayers+1 <<"." <<G4endl;
    } else {
      G4cerr <<" must be >0 and and <= " <<NbOfLayers <<"." <<G4endl;
    }
    G4cerr <<"--> Command rejected." <<G4endl;

  } else if (ival == 0) {
    FluxLayers.clear();
    for ( G4int i = 0; i <= NbOfLayers; i++) FluxLayers.push_back(i);
  } else if (std::find(FluxLayers.begin(),FluxLayers.end(),ival-1)
    != FluxLayers.end()) {
    G4cerr <<G4endl;
    G4cerr <<"AddFluxLayer : Layer index " <<ival <<" already selected."
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  } else {
    FluxLayers.push_back(ival-1);

  }
}
////////////////////////////////////////////////////////////////////////////////
//
//
void MLGeometryConstruction::DeleteFluxLayer (G4int ival)
{
  // search the material by its name
  //
  if (ival > NbOfLayers+1) {
    G4cerr <<G4endl;
    G4cerr <<"DeleteFluxLayer: layer index " <<ival <<" out of range." <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else  if (ival == 0) {
    FluxLayers.clear();

  } else if (std::find(FluxLayers.begin(),FluxLayers.end(),ival-1)
      == FluxLayers.end()) {
    G4cerr <<G4endl;
    G4cerr <<"DeleteFluxLayer : layer index " <<ival
           <<" has not been selected." <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else {
    FluxLayers.erase(std::
      find(FluxLayers.begin(),FluxLayers.end(),ival-1));
    MLNIELAnalyser* nielAnalyser = MLAnalysisManager::getInstance()
      ->GetNIELAnalyser();
    nielAnalyser->DeleteLayer(ival-1);

  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLGeometryConstruction::ListFluxLayer ()
{
  G4cout <<G4endl;
  G4cout <<"-------------------------------------------------------------"
         <<G4endl;
  G4cout <<"Particle fluences to be recorded for the following layers: "
         <<G4endl;

  if (FluxLayers.size() == 0) {
    G4cout <<"No layer has been selected for fluence analysis." <<G4endl;
  } else {
    for (size_t k=0; k<FluxLayers.size(); k++) {
      G4int i = FluxLayers[k];
      if (i == NbOfLayers) {
        G4cout <<"  Boundary before layer No. " <<std::setw(3) <<i+1
               <<G4endl;
      } else {
        G4cout <<"  Boundary before layer No. " <<std::setw(3) <<i+1 <<"  "
               <<std::setw(20) <<materialsManager
                 ->GetMaterial(G4int(Layers[i].y()))->GetName() <<": "
               <<std::setw(6) <<G4BestUnit(Layers[i].x(),"Length")
               <<G4endl;
      }
    }
  }
  G4cout <<"-------------------------------------------------------------"
         <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
G4bool MLGeometryConstruction::IsAFluxDetector (G4int ival)
{
  return (std::find(FluxLayers.begin(),FluxLayers.end(),ival)
    != FluxLayers.end());
}
////////////////////////////////////////////////////////////////////////////////
//
const G4int MLGeometryConstruction::GetELayerPosition (G4int ival)
{
  for (size_t i = 0; i < EdepLayers.size(); i++) {
    if ( EdepLayers[i] == ival ) return G4int(i);
  }
  G4cout <<"GetELayerPosition: " <<ival <<" not in the list." <<G4endl;
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
//
const G4int MLGeometryConstruction::GetFLayerPosition(G4int ival)
{
  for (size_t i = 0; i < FluxLayers.size(); i++) {
    if ( FluxLayers[i] == ival ) return G4int(i);
  }
  G4cout <<"GetFLayerPosition: " <<ival << " not in the list." <<G4endl;
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
RPTofstream & operator << (RPTofstream &RPTFile, MLGeometryConstruction &)
{
//
//
// Initialise some variables.
//
  G4RunManager *runManager = G4RunManager::GetRunManager();
  MLGeometryConstruction *geometry =
    (MLGeometryConstruction*)(runManager->GetUserDetectorConstruction());

  G4VPhysicalVolume *tempPV      = NULL;
  G4PhysicalVolumeStore *PVStore = NULL;
  PVStore                        = G4PhysicalVolumeStore::GetInstance();
  G4double halfz                 = 0.;

  std::vector <G4Material*> mList ;
  G4Material* aMat = NULL;
  G4int k = 0;
  G4int i;
  for (i=0; i<G4int(PVStore->size());i++) {
    tempPV = (*PVStore)[i];
    aMat   = tempPV->GetLogicalVolume()->GetMaterial();
    k      = count(mList.begin(),mList.end(),aMat);
    if (k==0) mList.push_back(aMat);
  }

  RPTFile <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
  RPTFile <<"Material Definition:" <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
  RPTFile <<G4endl;

  RPTFile <<"There are " <<mList.size() <<" materials used:" <<G4endl;

  for (i=0; i<G4int(mList.size()); i++) {
    RPTFile <<mList[i] <<G4endl;
  }

  RPTFile <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
  RPTFile <<"Geometry Definition:" <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
  RPTFile <<G4endl;
  RPTFile <<"There are " <<PVStore->size() <<" physical volumes used"
          <<" (including the world volume which is PhysVol #1)."
          <<G4endl;

  RPTFile <<G4endl;
  if (geometry->GetShape() == SLAB) {
    RPTFile <<"PhyVol#  PhyVol Name     Start      Thickness  Material"
            <<G4endl;
    for (i=1; i<G4int(PVStore->size()); i++) {
      RPTFile <<std::setw(7) <<i+1 <<"  ";
      RPTFile.outG4String((*PVStore)[i]->GetName(),15);

      RPTFile <<" ";
      halfz = ((G4Box*) ((*PVStore)[i]->GetLogicalVolume()->GetSolid()))
        ->GetZHalfLength();
      RPTFile.outG4BestUnit(G4BestUnit((*PVStore)[i]
        ->GetObjectTranslation().z()-halfz,"Length"),10);
      RPTFile <<" ";
      RPTFile.outG4BestUnit(G4BestUnit(halfz*2., "Length"),10);
      RPTFile <<" " <<std::setw(15) <<(*PVStore)[i]->GetLogicalVolume()
                ->GetMaterial()->GetName()<<G4endl;
    }
  } else {
    RPTFile <<"PhyVol#  PhyVol Name     Inner Radius  Outer Radius  Material"
            <<G4endl;
    for (i=1; i<G4int(PVStore->size()); i++) {
      RPTFile <<std::setw(7) <<i+1 <<"  ";
      RPTFile.outG4String((*PVStore)[i]->GetName(),15);

      RPTFile <<" ";
      RPTFile.outG4BestUnit(G4BestUnit(((G4Sphere*) ((*PVStore)[i]
        ->GetLogicalVolume()->GetSolid()))->GetInsideRadius(),"Length"),12);
      RPTFile <<" ";
      RPTFile.outG4BestUnit(G4BestUnit(((G4Sphere*) ((*PVStore)[i]
        ->GetLogicalVolume()->GetSolid()))->GetOuterRadius(),"Length"),12);
      RPTFile <<" " <<std::setw(15) <<(*PVStore)[i]->GetLogicalVolume()
                ->GetMaterial()->GetName()<<G4endl;
    }
  }
  return RPTFile;
}
#endif
////////////////////////////////////////////////////////////////////////////////
