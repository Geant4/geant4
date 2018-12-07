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
/// \file biasing/B02/src/B02ImportanceDetectorConstruction.cc
/// \brief Implementation of the B02ImportanceDetectorConstruction class
//
//
//

#include "globals.hh"
#include <sstream>

#include "B02ImportanceDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// For Primitive Scorers
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSNofCollision.hh"
#include "G4PSPopulation.hh"
#include "G4PSTrackCounter.hh"
#include "G4PSTrackLength.hh"

// for importance biasing
#include "G4IStore.hh"

// for weight window technique
#include "G4WeightWindowStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B02ImportanceDetectorConstruction::
B02ImportanceDetectorConstruction(G4String worldName) 
:G4VUserParallelWorld(worldName),fLogicalVolumeVector()
{
  //  Construct();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B02ImportanceDetectorConstruction::~B02ImportanceDetectorConstruction()
{
  fLogicalVolumeVector.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B02ImportanceDetectorConstruction::Construct()
{  
  G4cout << " constructing parallel world " << G4endl;

  G4Material* dummyMat  = 0;

  //GetWorld methods create a clone of the mass world to the parallel world (!)
  // via the transportation manager
  fGhostWorld = GetWorld();
  G4cout << " B02ImportanceDetectorConstruction:: ghostWorldName = " 
         << fGhostWorld->GetName() << G4endl;
  G4LogicalVolume* worldLogical = fGhostWorld->GetLogicalVolume();
  fLogicalVolumeVector.push_back(worldLogical);

  //  fPVolumeStore.AddPVolume(G4GeometryCell(*pWorldVolume, 0));
  fPVolumeStore.AddPVolume(G4GeometryCell(*fGhostWorld, 0));

  // creating 18 slobs of 10 cm thicknes

  G4double innerRadiusShield = 0*cm;
  G4double outerRadiusShield = 100*cm;
  G4double heightShield       = 5*cm;
  G4double startAngleShield  = 0*deg;
  G4double spanningAngleShield    = 360*deg;

  G4Tubs *aShield = new G4Tubs("aShield",
                               innerRadiusShield,
                               outerRadiusShield,
                               heightShield,
                               startAngleShield,
                               spanningAngleShield);
  
  // logical parallel cells

  G4LogicalVolume *aShield_log_imp = 
    new G4LogicalVolume(aShield, dummyMat, "aShield_log_imp");
  fLogicalVolumeVector.push_back(aShield_log_imp);

  // physical parallel cells
  G4String name = "none";
  G4int i = 1;
  G4double startz = -85*cm; 
  //  for (i=1; i<=18; ++i) {
  for (i=1; i<=18; i++) {
   
    name = GetCellName(i);
    
    G4double pos_x = 0*cm;
    G4double pos_y = 0*cm;
    G4double pos_z = startz + (i-1) * (2*heightShield);
    G4VPhysicalVolume *pvol = 
      new G4PVPlacement(0, 
                        G4ThreeVector(pos_x, pos_y, pos_z),
                        aShield_log_imp, 
                        name, 
                        worldLogical, 
                        false, 
                        i);
    //                        0);
    G4GeometryCell cell(*pvol, i);
    //    G4GeometryCell cell(*pvol, 0);
    fPVolumeStore.AddPVolume(cell);
  }

  // filling the rest of the world volumr behind the concrete with
  // another slob which should get the same importance value as the 
  // last slob
  innerRadiusShield = 0*cm;
  //  outerRadiusShield = 110*cm; exceeds world volume!!!!
  outerRadiusShield = 100*cm;
  //  heightShield       = 10*cm;
  heightShield       = 5*cm;
  startAngleShield  = 0*deg;
  spanningAngleShield    = 360*deg;

  G4Tubs *aRest = new G4Tubs("Rest",
                             innerRadiusShield,
                             outerRadiusShield,
                             heightShield,
                             startAngleShield,
                             spanningAngleShield);
  
  G4LogicalVolume *aRest_log = 
    new G4LogicalVolume(aRest, dummyMat, "aRest_log");

  fLogicalVolumeVector.push_back(aRest_log);

  name = GetCellName(19);
    
  G4double pos_x = 0*cm;
  G4double pos_y = 0*cm;
  //  G4double pos_z = 100*cm;
  G4double pos_z = 95*cm;
  G4VPhysicalVolume *pvol = 
    new G4PVPlacement(0, 
                      G4ThreeVector(pos_x, pos_y, pos_z),
                      aRest_log, 
                      name, 
                      worldLogical, 
                      false, 
                      19);
  //                      0);
  G4GeometryCell cell(*pvol, 19);
  //  G4GeometryCell cell(*pvol, 0);
  fPVolumeStore.AddPVolume(cell);

  SetSensitive();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4VPhysicalVolume &B02ImportanceDetectorConstruction::
GetPhysicalVolumeByName(const G4String& name) const {
  return *fPVolumeStore.GetPVolume(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String B02ImportanceDetectorConstruction::ListPhysNamesAsG4String(){
  G4String names(fPVolumeStore.GetPNames());
  return names;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String B02ImportanceDetectorConstruction::GetCellName(G4int i) {
  std::ostringstream os;
  os << "cell_";
  if (i<10) {
    os << "0";
  }
  os << i;
  G4String name = os.str();
  return name;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4GeometryCell B02ImportanceDetectorConstruction::GetGeometryCell(G4int i){
  G4String name(GetCellName(i));
  const G4VPhysicalVolume *p=0;
  p = fPVolumeStore.GetPVolume(name);
  if (p) {
    return G4GeometryCell(*p,0);
  }
  else {
    G4cout << "B02ImportanceDetectorConstruction::GetGeometryCell: " << G4endl
           << " couldn't get G4GeometryCell" << G4endl;
    return G4GeometryCell(*fGhostWorld,-2);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume &B02ImportanceDetectorConstruction::
GetWorldVolumeAddress() const{
   return *fGhostWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *B02ImportanceDetectorConstruction::GetWorldVolume() {
  return fGhostWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B02ImportanceDetectorConstruction::SetSensitive(){

  //  -------------------------------------------------
  //   The collection names of defined Primitives are
  //   0       ConcreteSD/Collisions
  //   1       ConcreteSD/CollWeight
  //   2       ConcreteSD/Population
  //   3       ConcreteSD/TrackEnter
  //   4       ConcreteSD/SL
  //   5       ConcreteSD/SLW
  //   6       ConcreteSD/SLWE
  //   7       ConcreteSD/SLW_V
  //   8       ConcreteSD/SLWE_V
  //  -------------------------------------------------

  // moved to ConstructSD() for MT compliance

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void B02ImportanceDetectorConstruction::ConstructSD()
{

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  //
  // Sensitive Detector Name
  G4String concreteSDname = "ConcreteSD";

  //------------------------
  // MultiFunctionalDetector
  //------------------------
  //
  // Define MultiFunctionalDetector with name.
  G4MultiFunctionalDetector* MFDet = 
                         new G4MultiFunctionalDetector(concreteSDname);
  SDman->AddNewDetector( MFDet );                 // Register SD to SDManager

  G4String fltName,particleName;
  G4SDParticleFilter* neutronFilter = 
      new G4SDParticleFilter(fltName="neutronFilter", particleName="neutron");

  MFDet->SetFilter(neutronFilter);

  for (std::vector<G4LogicalVolume *>::iterator it =  
                                                fLogicalVolumeVector.begin();
       it != fLogicalVolumeVector.end(); it++){
    //      (*it)->SetSensitiveDetector(MFDet);
    SetSensitiveDetector((*it)->GetName(), MFDet);
  }

  G4String psName;
  G4PSNofCollision*   scorer0 = new G4PSNofCollision(psName="Collisions");  
  MFDet->RegisterPrimitive(scorer0);

  G4PSNofCollision*   scorer1 = new G4PSNofCollision(psName="CollWeight");  
  scorer1->Weighted(true);
  MFDet->RegisterPrimitive(scorer1);

  G4PSPopulation*   scorer2 = new G4PSPopulation(psName="Population");  
  MFDet->RegisterPrimitive(scorer2);

  G4PSTrackCounter* scorer3 = 
                new G4PSTrackCounter(psName="TrackEnter",fCurrent_In);  
  MFDet->RegisterPrimitive(scorer3);

  G4PSTrackLength* scorer4 = new G4PSTrackLength(psName="SL");  
  MFDet->RegisterPrimitive(scorer4);

  G4PSTrackLength* scorer5 = new G4PSTrackLength(psName="SLW");  
  scorer5->Weighted(true);
  MFDet->RegisterPrimitive(scorer5);

  G4PSTrackLength* scorer6 = new G4PSTrackLength(psName="SLWE");  
  scorer6->Weighted(true);
  scorer6->MultiplyKineticEnergy(true);
  MFDet->RegisterPrimitive(scorer6);

  G4PSTrackLength* scorer7 = new G4PSTrackLength(psName="SLW_V");  
  scorer7->Weighted(true);
  scorer7->DivideByVelocity(true);
  MFDet->RegisterPrimitive(scorer7);

  G4PSTrackLength* scorer8 = new G4PSTrackLength(psName="SLWE_V");  
  scorer8->Weighted(true);
  scorer8->MultiplyKineticEnergy(true);
  scorer8->DivideByVelocity(true);
  MFDet->RegisterPrimitive(scorer8);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VIStore* B02ImportanceDetectorConstruction::CreateImportanceStore()
{

  G4cout << " B02ImportanceDetectorConstruction:: Creating Importance Store " 
         << G4endl;
  if (!fPVolumeStore.Size())
  {
    G4Exception("B02ImportanceDetectorConstruction::CreateImportanceStore"
               ,"exampleB02_0001",RunMustBeAborted
               ,"no physical volumes created yet!");
  }

  // creating and filling the importance store
  
  //  G4IStore *istore = new G4IStore(*fWorldVolume);
  
  G4IStore *istore = G4IStore::GetInstance(GetName());

  G4GeometryCell gWorldVolumeCell(GetWorldVolumeAddress(), 0);

  G4double imp =1;

  istore->AddImportanceGeometryCell(1, gWorldVolumeCell);

  // set importance values and create scorers 
  G4int cell(1);
  for (cell=1; cell<=18; cell++) {
    G4GeometryCell gCell = GetGeometryCell(cell);
    G4cout << " adding cell: " << cell 
           << " replica: " << gCell.GetReplicaNumber() 
           << " name: " << gCell.GetPhysicalVolume().GetName() << G4endl;
    imp = std::pow(2.0,cell-1);

    G4cout << "Going to assign importance: " << imp << ", to volume: " 
           << gCell.GetPhysicalVolume().GetName() << G4endl;
    //x    aIstore.AddImportanceGeometryCell(imp, gCell);
    istore->AddImportanceGeometryCell(imp, gCell.GetPhysicalVolume(), cell);
  }

  // creating the geometry cell and add both to the store
  //  G4GeometryCell gCell = GetGeometryCell(18);

  // create importance geometry cell pair for the "rest"cell
  // with the same importance as the last concrete cell 
  G4GeometryCell gCell = GetGeometryCell(19);
  //  G4double imp = std::pow(2.0,18); 
  imp = std::pow(2.0,17); 
  istore->AddImportanceGeometryCell(imp, gCell.GetPhysicalVolume(), 19);

  return istore;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VWeightWindowStore *B02ImportanceDetectorConstruction::
CreateWeightWindowStore()
{
  G4cout << " B02ImportanceDetectorConstruction:: Creating Importance Store " 
         << G4endl;
  if (!fPVolumeStore.Size())
  {
    G4Exception("B02ImportanceDetectorConstruction::CreateWeightWindowStore"
               ,"exampleB02_0002",RunMustBeAborted
               ,"no physical volumes created yet!");
  }

  // creating and filling the importance store
  
  //  G4IStore *istore = new G4IStore(*fWorldVolume);
  
  G4WeightWindowStore *wwstore = G4WeightWindowStore::GetInstance(GetName());

  // create one energy region covering the energies of the problem
  //
  std::set<G4double, std::less<G4double> > enBounds;
  enBounds.insert(1 * GeV);
  wwstore->SetGeneralUpperEnergyBounds(enBounds);
  
  G4int n = 0;
  G4double lowerWeight =1;
  std::vector<G4double> lowerWeights;

  lowerWeights.push_back(1);
  G4GeometryCell gWorldCell(GetWorldVolumeAddress(),0);
  wwstore->AddLowerWeights(gWorldCell, lowerWeights);

  G4int cell(1);
  for (cell=1; cell<=18; cell++) {
    G4GeometryCell gCell = GetGeometryCell(cell);
    G4cout << " adding cell: " << cell 
           << " replica: " << gCell.GetReplicaNumber() 
           << " name: " << gCell.GetPhysicalVolume().GetName() << G4endl;

      lowerWeight = 1./std::pow(2., n++);
      G4cout << "Going to assign lower weight: " << lowerWeight 
             << ", to volume: " 
             << gCell.GetPhysicalVolume().GetName() << G4endl;
      lowerWeights.clear();
      lowerWeights.push_back(lowerWeight);
      wwstore->AddLowerWeights(gCell, lowerWeights);
  }

  // the remaining part pf the geometry (rest) gets the same
  // lower weight bound  as the last conrete cell
  //

  // create importance geometry cell pair for the "rest"cell
  // with the same importance as the last concrete cell 
  G4GeometryCell gCell = GetGeometryCell(19);
  wwstore->AddLowerWeights(gCell,  lowerWeights);

  return wwstore;

}
