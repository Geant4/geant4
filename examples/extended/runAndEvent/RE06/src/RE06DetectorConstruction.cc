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
/// \file RE06/src/RE06DetectorConstruction.cc
/// \brief Implementation of the RE06DetectorConstruction class
//
// 

#include "RE06DetectorConstruction.hh"

#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSNofSecondary.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofStep.hh"
#include "G4PSMinKinEAtGeneration.hh"
#include "G4VSDFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4ios.hh"

#include "RE06DetectorMessenger.hh"
#include "RE06PrimaryGeneratorAction.hh"
#include "RE06ParallelWorld.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4bool RE06DetectorConstruction::fConstructedSDandField = false;

RE06DetectorConstruction::RE06DetectorConstruction()
: G4VUserDetectorConstruction(),
  fNumberOfLayers(40),
  fTotalThickness (2.0*m),
  fLayerThickness(0.),
  fConstructed(false),
  fWorldMaterial(0),
  fAbsorberMaterial(0),
  fGapMaterial(0),
  fLayerSolid(0),
  fGapSolid(0),
  fWorldLogical(0),
  fWorldPhysical(0),
  fSerial(false),
  fDetectorMessenger(0),
  fVerboseLevel(1)
{
  fLayerThickness = fTotalThickness / fNumberOfLayers;

  for(size_t i=0;i<3;i++)
  {
    fCalorLogical[i] = 0;
    fLayerLogical[i] = 0;
    fGapLogical[i] = 0;
    fCalorPhysical[i] = 0;
    fLayerPhysical[i] = 0;
    fGapPhysical[i] = 0;
  }

  fCalName[0] = "Calor-A";
  fCalName[1] = "Calor-B";
  fCalName[2] = "Calor-C";

  fDetectorMessenger = new RE06DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE06DetectorConstruction::~RE06DetectorConstruction()
{ delete fDetectorMessenger;}

G4VPhysicalVolume* RE06DetectorConstruction::Construct()
{
  if(!fConstructed)
  {
    fConstructed = true;
    DefineMaterials();
    SetupGeometry();
  }
  if (GetVerboseLevel()>0) {
    PrintCalorParameters();
  }
  return fWorldPhysical;
}

void RE06DetectorConstruction::ConstructSDandField() 
{
  if(!fConstructedSDandField)
  {
    fConstructedSDandField = true;
    SetupDetectors();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06DetectorConstruction::DefineMaterials()
{ 
  G4String name, symbol;             //a=mass of a mole;
  G4double a, z, density;            //z=mean number of protons;  
  G4int iz;                          //iz=number of protons  in an isotope; 
  G4int n;                           // n=number of nucleons in an isotope;

  G4int ncomponents, natoms;
  G4double abundance, fractionmass;
  G4double temperature, pressure;

  //
  // define Elements
  //

  a = 1.01*g/mole;
  G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

  a = 12.01*g/mole;
  G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  a = 14.01*g/mole;
  G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  //
  // define an Element from isotopes, by relative abundance 
  //

  G4Isotope* U5 = new G4Isotope(name="U235", iz=92, n=235, a=235.01*g/mole);
  G4Isotope* U8 = new G4Isotope(name="U238", iz=92, n=238, a=238.03*g/mole);

  G4Element* U  = new G4Element(name="enriched Uranium",symbol="U",ncomponents=2);
  U->AddIsotope(U5, abundance= 90.*perCent);
  U->AddIsotope(U8, abundance= 10.*perCent);

  //
  // define simple materials
  //

  new G4Material(name="Aluminium", z=13., a=26.98*g/mole, density=2.700*g/cm3);
  new G4Material(name="Silicon", z=14., a= 28.09*g/mole, density= 2.33*g/cm3);
  new G4Material(name="Iron", z=26., a=55.85*g/mole, density=7.87*g/cm3);
  new G4Material(name="ArgonGas",z=18., a= 39.95*g/mole, density=1.782*mg/cm3);
  new G4Material(name="He", z=2., a=4.0*g/mole, density=0.1786e-03*g/cm3);

  density = 1.390*g/cm3;
  a = 39.95*g/mole;
  G4Material* lAr = new G4Material(name="liquidArgon", z=18., a, density);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  G4Material* Pb = new G4Material(name="Lead"     , z=82., a, density);

  //
  // define a material from elements.   case 1: chemical molecule
  //
 
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);

  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(C, natoms=9);
  Sci->AddElement(H, natoms=10);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  //
  // examples of vacuum
  //

  density     = universe_mean_density;
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material* Vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole,
                                    density,kStateGas,temperature,pressure);

  if (GetVerboseLevel()>1) {
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  }

  //default materials of the calorimeter
  fWorldMaterial    = Vacuum;
  fAbsorberMaterial = Pb;
  fGapMaterial      = lAr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06DetectorConstruction::SetupGeometry()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",2.*m,2.*m,fTotalThickness*2.);
  fWorldLogical = new G4LogicalVolume(worldSolid,fWorldMaterial,"World");
  fWorldPhysical = new G4PVPlacement(0,G4ThreeVector(),fWorldLogical,"World",
                        0,false,0);
  
  //                               
  // Calorimeter
  //  
  G4VSolid* calorSolid = new G4Box("Calor",0.5*m,0.5*m,fTotalThickness/2.);
  G4int i;
  for(i=0;i<3;i++)
  {
    fCalorLogical[i] 
      = new G4LogicalVolume(calorSolid,fAbsorberMaterial,fCalName[i]);
    if(fSerial)
    {
      fCalorPhysical[i] = new G4PVPlacement(0,
                 G4ThreeVector(0.,0.,G4double(i-1)*fTotalThickness),
                 fCalorLogical[i],fCalName[i],fWorldLogical,false,i);
    }
    else
    {
      fCalorPhysical[i] = new G4PVPlacement(0,
                 G4ThreeVector(0.,G4double(i-1)*m,0.),
                 fCalorLogical[i],fCalName[i],fWorldLogical,false,i);
    }
  }
 
  //                                 
  // Layers --- as absorbers
  //
  fLayerSolid = new G4Box("Layer",0.5*m,0.5*m,fLayerThickness/2.);
  for(i=0;i<3;i++)
  {
    fLayerLogical[i] 
      = new G4LogicalVolume(fLayerSolid,fAbsorberMaterial,fCalName[i]+"_LayerLog");
    fLayerPhysical[i] 
      = new G4PVReplica(fCalName[i]+"_Layer",fLayerLogical[i],fCalorLogical[i],
                        kZAxis,fNumberOfLayers,fLayerThickness);
  }
   
  //
  // Gap
  //
  fGapSolid = new G4Box("Gap",0.5*m,0.5*m,fLayerThickness/4.);
  for(i=0;i<3;i++)
  {
    fGapLogical[i] = new G4LogicalVolume(fGapSolid,fGapMaterial,fCalName[i]+"_Gap");
    fGapPhysical[i] = new G4PVPlacement(0,G4ThreeVector(0.,0.,fLayerThickness/4.),
                fGapLogical[i],fCalName[i]+"_gap",fLayerLogical[i],false,0);
  }

  //
  // Regions
  //
  for(i=0;i<3;i++)
  {
    G4Region* aRegion = new G4Region(fCalName[i]);
    fCalorLogical[i]->SetRegion(aRegion);
    aRegion->AddRootLogicalVolume(fCalorLogical[i]);
  }

  //                                        
  // Visualization attributes
  //
  fWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  for(i=0;i<3;i++)
  { 
    fCalorLogical[i]->SetVisAttributes(simpleBoxVisAtt);
    fLayerLogical[i]->SetVisAttributes(simpleBoxVisAtt);
    fGapLogical[i]->SetVisAttributes(simpleBoxVisAtt);
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06DetectorConstruction::SetupDetectors()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  G4String filterName, particleName;

  G4SDParticleFilter* gammaFilter 
    = new G4SDParticleFilter(filterName="gammaFilter",particleName="gamma");
  G4SDParticleFilter* electronFilter 
    = new G4SDParticleFilter(filterName="electronFilter",particleName="e-");
  G4SDParticleFilter* positronFilter 
    = new G4SDParticleFilter(filterName="positronFilter",particleName="e+");
  G4SDParticleFilter* epFilter 
    = new G4SDParticleFilter(filterName="epFilter");
  epFilter->add(particleName="e-");
  epFilter->add(particleName="e+");


  for(G4int i=0;i<3;i++)
  {
   for(G4int j=0;j<2;j++)
   {
    // Loop counter j = 0 : absorber
    //                = 1 : gap
    G4String detName = fCalName[i];
    if(j==0)
    { detName += "_abs"; }
    else
    { detName += "_gap"; }
    G4MultiFunctionalDetector* det = new G4MultiFunctionalDetector(detName);
    G4SDManager::GetSDMpointer()->AddNewDetector(det);

    // The second argument in each primitive means the "level" of geometrical 
    // hierarchy, the copy number of that level is used as the key of the 
    // G4THitsMap.
    // For absorber (j = 0), the copy number of its own physical volume is used.
    // For gap (j = 1), the copy number of its mother physical volume is used, 
    // since there is only one physical volume of gap is placed with respect 
    // to its mother.
    G4VPrimitiveScorer* primitive;
    primitive = new G4PSEnergyDeposit("eDep",j);
    det->RegisterPrimitive(primitive);
    primitive = new G4PSNofSecondary("nGamma",j);
    primitive->SetFilter(gammaFilter);
    det->RegisterPrimitive(primitive);
    primitive = new G4PSNofSecondary("nElectron",j);
    primitive->SetFilter(electronFilter);
    det->RegisterPrimitive(primitive);
    primitive = new G4PSNofSecondary("nPositron",j);
    primitive->SetFilter(positronFilter);
    det->RegisterPrimitive(primitive);
    primitive = new G4PSMinKinEAtGeneration("minEkinGamma",j);
    primitive->SetFilter(gammaFilter);
    det->RegisterPrimitive(primitive);
    primitive = new G4PSMinKinEAtGeneration("minEkinElectron",j);
    primitive->SetFilter(electronFilter);
    det->RegisterPrimitive(primitive);
    primitive = new G4PSMinKinEAtGeneration("minEkinPositron",j);
    primitive->SetFilter(positronFilter);
    det->RegisterPrimitive(primitive);
    primitive = new G4PSTrackLength("trackLength",j);
    primitive->SetFilter(epFilter);
    det->RegisterPrimitive(primitive);
    primitive = new G4PSNofStep("nStep",j);
    primitive->SetFilter(epFilter);
    det->RegisterPrimitive(primitive);

    if(j==0)
     { SetSensitiveDetector(fLayerLogical[i], det); }
    else
    { SetSensitiveDetector(fGapLogical[i], det);}
   }
  }
  G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
}

void RE06DetectorConstruction::PrintCalorParameters() const
{
  G4cout 
    << "--------------------------------------------------------" << G4endl;
  if(fSerial)
  { G4cout << " Calorimeters are placed in serial." << G4endl; }
  else
  { G4cout << " Calorimeters are placed in parallel." << G4endl; }
  G4cout 
    << " Absorber is made of " << fAbsorberMaterial->GetName() << G4endl
    << " Gap is made of " << fGapMaterial->GetName() << G4endl
    << "--------------------------------------------------------" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if(pttoMaterial)
  {
    fAbsorberMaterial = pttoMaterial;
    if(fConstructed) for(size_t i=0;i<3;i++)
    {
      fCalorLogical[i]->SetMaterial(fAbsorberMaterial);
      fLayerLogical[i]->SetMaterial(fAbsorberMaterial);
    }
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    if (GetVerboseLevel()>1) {
      PrintCalorParameters();
    }
  }
  else
  { 
    G4cerr 
      << materialChoice << " is not defined. - Command is ignored." << G4endl; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String RE06DetectorConstruction::GetAbsorberMaterial() const
{ return fAbsorberMaterial->GetName(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06DetectorConstruction::SetGapMaterial(G4String materialChoice)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);  
  if(pttoMaterial)
  {
    fGapMaterial = pttoMaterial;
    if(fConstructed) for(size_t i=0;i<3;i++)
    { fGapLogical[i]->SetMaterial(fGapMaterial); }
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    if (GetVerboseLevel()>1) {
      PrintCalorParameters();
    }
   }
  else
  { 
    G4cerr 
      << materialChoice << " is not defined. - Command is ignored." << G4endl; 
  }
}

G4String RE06DetectorConstruction::GetGapMaterial() const
{ return fGapMaterial->GetName(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06DetectorConstruction::SetSerialGeometry(G4bool serial)
{
  if(fSerial==serial) return;
  fSerial=serial;
  RE06PrimaryGeneratorAction* gen = (RE06PrimaryGeneratorAction*)
          (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  if(gen) gen->SetSerial(fSerial);
  if(!fConstructed) return;
  for(G4int i=0;i<3;i++)
  {
    if(fSerial)
    { 
      fCalorPhysical[i]
        ->SetTranslation(G4ThreeVector(0.,0.,G4double(i-1)*2.*m));
    }
    else
    { 
      fCalorPhysical[i]
        ->SetTranslation(G4ThreeVector(0.,G4double(i-1)*m,0.)); 
    }
  }
  ((RE06ParallelWorld*)GetParallelWorld(0))->SetSerialGeometry(serial);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06DetectorConstruction::SetNumberOfLayers(G4int nl)
{
  fNumberOfLayers = nl;
  fLayerThickness = fTotalThickness/fNumberOfLayers;
  if(!fConstructed) return;

  fLayerSolid->SetZHalfLength(fLayerThickness/2.);
  fGapSolid->SetZHalfLength(fLayerThickness/4.);
  for(size_t i=0;i<3;i++)
  { 
    fCalorLogical[i]->RemoveDaughter(fLayerPhysical[i]);
    delete fLayerPhysical[i];
    fLayerPhysical[i] 
      = new G4PVReplica(fCalName[i]+"_Layer",fLayerLogical[i],fCalorLogical[i],
                        kZAxis,fNumberOfLayers,fLayerThickness);
    fGapPhysical[i]->SetTranslation(G4ThreeVector(0.,0.,fLayerThickness/4.));
  }
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void   RE06DetectorConstruction::AddMaterial()
{
  static G4bool isAdded = false;   

  if( isAdded ) return;

  G4String name, symbol;             //a=mass of a mole;
  G4double a, z, density;            //z=mean number of protons;  

  G4int ncomponents, natoms;

  //
  // define simple materials
  //

  new G4Material(name="Copper", z=29., a=63.546*g/mole, density=8.96*g/cm3);
  new G4Material(name="Tungsten", z=74., a=183.84*g/mole, density=19.3*g/cm3);

  G4Element* C = G4Element::GetElement("Carbon");
  G4Element* O = G4Element::GetElement("Oxygen");
  

  G4Material* CO2 = 
    new G4Material("CarbonicGas", density= 27.*mg/cm3, ncomponents=2,
                   kStateGas, 325.*kelvin, 50.*atmosphere);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);

  isAdded = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
