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
/// \file biasing/B01/src/B01DetectorConstruction.cc
/// \brief Implementation of the B01DetectorConstruction class
//
//
//

#include "G4Types.hh"
#include <sstream>
#include <set>
#include "globals.hh"

#include "B01DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
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

B01DetectorConstruction::B01DetectorConstruction() :
  G4VUserDetectorConstruction(),
  fLogicalVolumeVector(),fPhysicalVolumeVector()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B01DetectorConstruction::~B01DetectorConstruction()
{
  fLogicalVolumeVector.clear();
  fPhysicalVolumeVector.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B01DetectorConstruction::Construct()
{
  G4double pos_x;
  G4double pos_y;
  G4double pos_z; 

  G4double density, pressure, temperature;
  G4double A;
  G4int Z;

  G4String name, symbol;
  G4double z;
  G4double fractionmass;

  A = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , Z= 1, A);

  A = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , Z = 6, A);

  A = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , Z= 8, A);

  A = 22.99*g/mole; 
  G4Element* elNa  = new G4Element(name="Natrium"  ,symbol="Na" , Z=11 , A);

  A = 200.59*g/mole; 
  G4Element* elHg  = new G4Element(name="Hg"  ,symbol="Hg" , Z=80, A);

  A = 26.98*g/mole; 
  G4Element* elAl  = new G4Element(name="Aluminium"  ,symbol="Al" , Z=13, A);

  A = 28.09*g/mole;
  G4Element* elSi  = new G4Element(name="Silicon", symbol="Si", Z=14, A);

  A = 39.1*g/mole; 
  G4Element* elK  = new G4Element(name="K"  ,symbol="K" , Z=19 , A);

  A = 69.72*g/mole; 
  G4Element* elCa  = new G4Element(name="Calzium"  ,symbol="Ca" , Z=31 , A);

  A = 55.85*g/mole;
  G4Element* elFe = new G4Element(name="Iron"    ,symbol="Fe", Z=26, A);

  density     = universe_mean_density;            //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material *Galactic = 
    new G4Material(name="Galactic", z=1., A=1.01*g/mole, density,
                   kStateGas,temperature,pressure);

  density = 2.03*g/cm3;
  G4Material* Concrete = new G4Material("Concrete", density, 10);
  Concrete->AddElement(elH , fractionmass= 0.01);
  Concrete->AddElement(elO , fractionmass= 0.529);
  Concrete->AddElement(elNa , fractionmass= 0.016);
  Concrete->AddElement(elHg , fractionmass= 0.002);
  Concrete->AddElement(elAl , fractionmass= 0.034);
  Concrete->AddElement(elSi , fractionmass= 0.337);
  Concrete->AddElement(elK , fractionmass= 0.013);
  Concrete->AddElement(elCa , fractionmass= 0.044);
  Concrete->AddElement(elFe , fractionmass= 0.014);
  Concrete->AddElement(elC , fractionmass= 0.001);

  /////////////////////////////
  // world cylinder volume
  ////////////////////////////

  // world solid

  G4double innerRadiusCylinder = 0*cm;
  G4double outerRadiusCylinder = 100*cm; 
  G4double heightCylinder       = 100*cm;
  G4double startAngleCylinder  = 0*deg;
  G4double spanningAngleCylinder    = 360*deg;

  G4Tubs *worldCylinder = new G4Tubs("worldCylinder",
                                     innerRadiusCylinder,
                                     outerRadiusCylinder,
                                     heightCylinder,
                                     startAngleCylinder,
                                     spanningAngleCylinder);

  // logical world

  G4LogicalVolume *worldCylinder_log = 
    new G4LogicalVolume(worldCylinder, Galactic, "worldCylinder_log");
  fLogicalVolumeVector.push_back(worldCylinder_log); 

  name = "shieldWorld";
  fWorldVolume = new 
    G4PVPlacement(0, G4ThreeVector(0,0,0), worldCylinder_log,
                  name, 0, false, 0);

  fPhysicalVolumeVector.push_back(fWorldVolume);

  // creating 18 slabs of 10 cm thick concrete

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
  
  // logical shield

  G4LogicalVolume *aShield_log =
    new G4LogicalVolume(aShield, Concrete, "aShield_log");
  fLogicalVolumeVector.push_back(aShield_log);

  G4VisAttributes* pShieldVis = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  pShieldVis->SetForceSolid(true);
  aShield_log->SetVisAttributes(pShieldVis);

  // physical shields

  G4int i;
  G4double startz = -85*cm; 
  for (i=1; i<=18; i++)
  {
    name = GetCellName(i);
    pos_x = 0*cm;
    pos_y = 0*cm;
    pos_z = startz + (i-1) * (2*heightShield);
    G4VPhysicalVolume *pvol = 
      new G4PVPlacement(0, 
                        G4ThreeVector(pos_x, pos_y, pos_z),
                        aShield_log, 
                        name, 
                        worldCylinder_log, 
                        false, 
                        i); 
    fPhysicalVolumeVector.push_back(pvol);
  }

  // filling the rest of the world volume behind the concrete with
  // another slab which should get the same importance value 
  // or lower weight bound as the last slab
  //
  innerRadiusShield = 0*cm;
  outerRadiusShield = 100*cm;
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
    new G4LogicalVolume(aRest, Galactic, "aRest_log");
  fLogicalVolumeVector.push_back(aRest_log);
  name = "rest";
    
  pos_x = 0*cm;
  pos_y = 0*cm;
  pos_z = 95*cm;
  G4VPhysicalVolume *pvol_rest = 
    new G4PVPlacement(0, 
                      G4ThreeVector(pos_x, pos_y, pos_z),
                      aRest_log, 
                      name, 
                      worldCylinder_log, 
                      false, 
                      19); // i=19

  fPhysicalVolumeVector.push_back(pvol_rest);

  SetSensitive();
  return fWorldVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VIStore* B01DetectorConstruction::CreateImportanceStore()
{
  G4cout << " B01DetectorConstruction:: Creating Importance Store " << G4endl;
  if (!fPhysicalVolumeVector.size())
  {
    G4Exception("B01DetectorConstruction::CreateImportanceStore"
               ,"exampleB01_0001",RunMustBeAborted
               ,"no physical volumes created yet!");
  }

  fWorldVolume = fPhysicalVolumeVector[0];

  // creating and filling the importance store
  
  G4IStore *istore = G4IStore::GetInstance();

  G4int n = 0;
  G4double imp =1;
  istore->AddImportanceGeometryCell(1,  *fWorldVolume);
  for (std::vector<G4VPhysicalVolume *>::iterator
       it =  fPhysicalVolumeVector.begin();
       it != fPhysicalVolumeVector.end() - 1; it++)
  {
    if (*it != fWorldVolume)
    {
      imp = std::pow(2., n++);
      G4cout << "Going to assign importance: " << imp << ", to volume: " 
             << (*it)->GetName() << G4endl;
      istore->AddImportanceGeometryCell(imp, *(*it),n);
    }
  }

  // the remaining part pf the geometry (rest) gets the same
  // importance as the last conrete cell
  //
  istore->AddImportanceGeometryCell(imp, 
          *(fPhysicalVolumeVector[fPhysicalVolumeVector.size()-1]),++n);
  
  return istore;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VWeightWindowStore *B01DetectorConstruction::CreateWeightWindowStore()
{
  if (!fPhysicalVolumeVector.size())
  {
    G4Exception("B01DetectorConstruction::CreateWeightWindowStore"
               ,"exampleB01_0002",RunMustBeAborted
               ,"no physical volumes created yet!");
  }

  fWorldVolume = fPhysicalVolumeVector[0];

  // creating and filling the weight window store
  
  G4WeightWindowStore *wwstore = G4WeightWindowStore::GetInstance();
  
  // create one energy region covering the energies of the problem
  //
  std::set<G4double, std::less<G4double> > enBounds;
  enBounds.insert(1 * GeV);
  wwstore->SetGeneralUpperEnergyBounds(enBounds);

  G4int n = 0;
  G4double lowerWeight =1;
  std::vector<G4double> lowerWeights;

  lowerWeights.push_back(1);
  G4GeometryCell gWorldCell(*fWorldVolume,0);
  wwstore->AddLowerWeights(gWorldCell, lowerWeights);

  for (std::vector<G4VPhysicalVolume *>::iterator
       it =  fPhysicalVolumeVector.begin();
       it != fPhysicalVolumeVector.end() - 1; it++)
  {
    if (*it != fWorldVolume)
    {
      lowerWeight = 1./std::pow(2., n++);
      G4cout << "Going to assign lower weight: " << lowerWeight 
             << ", to volume: " 
             << (*it)->GetName() << G4endl;
      G4GeometryCell gCell(*(*it),n);
      lowerWeights.clear();
      lowerWeights.push_back(lowerWeight);
      wwstore->AddLowerWeights(gCell, lowerWeights);
    }
  }

  // the remaining part pf the geometry (rest) gets the same
  // lower weight bound  as the last conrete cell
  //
  G4GeometryCell
    gRestCell(*(fPhysicalVolumeVector[fPhysicalVolumeVector.size()-1]), ++n);
  wwstore->AddLowerWeights(gRestCell,  lowerWeights);

  return wwstore;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String B01DetectorConstruction::GetCellName(G4int i)
{
  std::ostringstream os;
  os << "cell_";
  if (i<10)
  {
    os << "0";
  }
  os << i ;
  G4String name = os.str();
  return name;
}

G4VPhysicalVolume *B01DetectorConstruction::GetWorldVolume() {
   return fWorldVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B01DetectorConstruction::SetSensitive(){

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

  // moved to ConstructSDandField() for MT compliance

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B01DetectorConstruction::ConstructSDandField()
{

  //  Sensitive Detector Manager.
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
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

  G4PSTrackCounter* scorer3 = new G4PSTrackCounter(psName="TrackEnter"
                                                  ,fCurrent_In);  
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
