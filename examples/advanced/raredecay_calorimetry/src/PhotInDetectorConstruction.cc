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
// $Id: PhotInDetectorConstruction.cc,v 1.5 2006/06/29 16:25:09 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// ***************************************************************************

//#define debug

#include "PhotInDetectorConstruction.hh"

// all BODY variables must be initialized in this constructor
PhotInDetectorConstruction::PhotInDetectorConstruction(G4double x,G4double y,G4double z):
 numberOfLayers(PhotInNOfLayers),numberOfSlabs(PhotInNOfSlabs),
 samplingFraction(PhotInSampFract),xHD(x),yHD(y),zHD(z),serial(false),
 worldMaterial(0),absorberMaterial(0),gapMaterial(0),layerSolid(0),slabSolid(0)
{
#ifdef debug
  G4cout<<"PhotInDetectorConstruction::Constructor is called"<<G4endl;
#endif
  layerParam = new PhotInLayerParameterisation;    // should it be deleted in Distructor?MK
  layerParam->SetNumberOfLayers(numberOfLayers);   // 9 by default (can be changed)
  layerParam->SetHalfTotalThickness(zHD);          // transfer totalThickness of theSection
  G4double hlth=zHD/numberOfLayers;                // Calculate thickness of one layer
  gapParam = new PhotInGapParameterisation;        // should it be deleted in Distructor?MK
  gapParam->SetNumberOfSlabs(numberOfSlabs);       // 9 by default (can be changed)
  gapParam->SetHalfTotalWidth(yHD);                // Transfer
  gapParam->SetHalfTotalLayerThickness(hlth);      // Transfer
  gapParam->SetSamplingFraction(samplingFraction); // Transfer
  PhotInCalorimeterSD::SetNumberOfLayers(numberOfLayers);// Construction of sensitiveLayers
  PhotInCalorimeterSD::SetNumberOfSlabs(numberOfSlabs);  // Construction of sensitiveSlabs
  DefineMaterials();                               // Material factory

  for(G4int i=0; i<PhotInNumSections ;i++)
  {
    calorLogical[i]=0;
    layerLogical[i]=0;
    slabLogical[i]=0;
    calorPhysical[i]=0;
    layerPhysical[i]=0;
    slabPhysical[i]=0;
  }
}

PhotInDetectorConstruction::~PhotInDetectorConstruction()
{
  //delete gapParam;           // Try to open later to be sure that it is not delited by G4
  //delete layerParam;         // Try to open later to be sure that it is not delited by G4
}

void PhotInDetectorConstruction::DefineMaterials() // Material factory (@@ can be in Const)
{ 
  G4String name, symbol;             // a = mean_atomic_mass, z = mean_number_of_protons
  G4double a, z, density;            //iz=number of protons  in an isotope; 
  G4int iz, in;                      //in=number of nucleons in an isotope;
                                    
  G4int ncomponents, natoms;
  G4double abundance, fractionmass;
  G4double temperature, pressure;
  //
  // define Elements
  //
  G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a=1.01*g/mole);
  G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a=12.01*g/mole);
  G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a=14.01*g/mole);
  G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a=16.00*g/mole);
  //
  // define an Element from isotopes, by relative abundance 
  //
  G4Isotope* U5 = new G4Isotope(name="U235", iz=92, in=235, a=235.01*g/mole); // @@ Corr.
  G4Isotope* U8 = new G4Isotope(name="U238", iz=92, in=238, a=238.03*g/mole); // @@ Corr.

  G4Element* U  = new G4Element(name="enriched Uranium", symbol="U", ncomponents=2);
  U->AddIsotope(U5, abundance= 90.*perCent);
  U->AddIsotope(U8, abundance= 10.*perCent);
  //
  // define simple materials
  //
  G4Material* Al = new G4Material(name="Aluminium",z=13.,a=26.98*g/mole,density=2.7*g/cm3);

  //density = 1.390*g/cm3;
  //a = 39.95*g/mole;
  //G4Material* lAr = new G4Material(name="liquidArgon", z=18., a, density);

  //density = 11.35*g/cm3;
  //a = 207.19*g/mole;
  //G4Material* Pb = new G4Material(name="Lead"     , z=82., a, density);
  //
  // define complex materials. case 1: chemical molecule
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
  // define complex material.   case 2: mixture by fractional mass
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
  a=1.01*g/mole; // Kosmic
		z=1.;          // Hydrogrn
  G4Material* Vacuum=new G4Material(name="Vac",z,a,density,kStateGas,temperature,pressure);

#ifdef debug
  G4cout<<"PhotInDetectorConstruction::DefineMaterials:"<<*(G4Material::GetMaterialTable())
        <<G4endl;
#endif

  //default materials of the calorimeter
  worldMaterial    = Vacuum;
  absorberMaterial = Al;
  gapMaterial      = Sci;
  layerParam->SetAbsorberMaterial(absorberMaterial);
  gapParam->SetGapMaterial(gapMaterial);
}

G4VPhysicalVolume* PhotInDetectorConstruction::Construct()
{
  // Definition of Logical and Physical volumes
  G4LogicalVolume*           worldLogical;     // One logical volume for the World
  G4VPhysicalVolume*         worldPhysical;    // One Physical volume for the World

  G4double shiftY=yHD+yHD;        // shift perpendicular to slabs
  G4double shiftZ=zHD+zHD;        // shift perpendicular to layers
  G4double wxHD=xHD+0.5*m;        // along slabs
  G4double wyHD=shiftY+yHD+0.5*m; // perpendicular to slabs
  G4double wzHD=shiftZ+zHD+0.5*m; // perpendicular to layers
  //     
  // World volume: Solid=BOX, make Logical, make Physical
  // ====================================================
  G4VSolid* worldSolid = new G4Box("World",wxHD,wyHD,wzHD);
  worldLogical = new G4LogicalVolume(worldSolid,worldMaterial,"World");
  worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,false,0);
  //                               
  // Calorimeter Section: Solid=BOX, make 3 Logical and 3 serial/parallel Physical
  // =============================================================================
  G4VSolid* calorSolid = new G4Box("CalorSect",xHD,yHD,zHD);
  G4int i; // If i is used many times it is btter to define it as external var. for LOOPs
  for(i=0; i<PhotInNumSections; i++)
  {
    calorLogical[i] = new G4LogicalVolume(calorSolid,absorberMaterial,PhotInCalName[i]);
    if(serial)                              // if true = one section after another (serial)
      calorPhysical[i] = new G4PVPlacement(0, G4ThreeVector(0.,0.,(i-1)*shiftZ),
                                    calorLogical[i],PhotInCalName[i],worldLogical,false,i);
    else                                    // *DEFAULT* false = side by side (parallel)
      calorPhysical[i] = new G4PVPlacement(0, G4ThreeVector(0.,G4double(i-1)*shiftY,0.),
                                    calorLogical[i],PhotInCalName[i],worldLogical,false,i);
  }
  //                                 
  // Layers & Slabs: Solid=BOX. First Active (samplFract thick div in slabs), then Absorber
  //
  G4double layerZ=zHD/numberOfLayers;     // Half thickness of one layer
  G4double slabY=yHD/numberOfSlabs;       // Half width of one active slab
  layerSolid = new G4Box("Layer",xHD,yHD,layerZ);
  slabSolid = new G4Box("Gap",xHD,slabY,layerZ*samplingFraction);
  for(i=0; i<PhotInNumSections; i++)                   // Construct layers for all sections
  {
    layerLogical[i] = new G4LogicalVolume(layerSolid,absorberMaterial,"Layer");
    layerPhysical[i] = new G4PVParameterised("Layer", layerLogical[i], calorLogical[i],
                                             kZAxis, numberOfLayers, layerParam);
    slabLogical[i] = new G4LogicalVolume(slabSolid,gapMaterial,"Slab");
    slabPhysical[i] = new G4PVParameterised("Slab", slabLogical[i], layerLogical[i],
                                             kYAxis, numberOfSlabs, gapParam);
  }
  //
  // Regions
  //
  for(i=0; i<PhotInNumSections; i++)
  {
    G4Region* aRegion = new G4Region(PhotInRegName[i]);
    calorLogical[i]->SetRegion(aRegion);            // Mutual
    aRegion->AddRootLogicalVolume(calorLogical[i]); // definition
  }
  //                               
  // Sensitive Detectors: Absorber and Gap
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer(); // Manager of sensitive detectors
  for(i=0; i<PhotInNumSections ;i++)
  {
    G4VSensitiveDetector* calorSD = new PhotInCalorimeterSD(PhotInDetName[i]); // Sections
    SDman->AddNewDetector(calorSD); // all section is a sensitive detector
    layerLogical[i]->SetSensitiveDetector(calorSD);// Make layers to be sensitive detectors
    slabLogical[i]->SetSensitiveDetector(calorSD); // Make slabs to be sensitive detectors
  }
  //                                        
  // Visualization attributes
  //
  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.,1.,1.));
  simpleBoxVisAtt->SetVisibility(true);
  for(i=0;i<PhotInNumSections;i++)
  { 
    calorLogical[i]->SetVisAttributes(simpleBoxVisAtt);
    layerLogical[i]->SetVisAttributes(simpleBoxVisAtt);
    slabLogical[i]->SetVisAttributes(simpleBoxVisAtt);
  }
#ifdef debug
  G4cout<<"PhotInDetectorConstruction::Construct:"<<G4endl;
  PrintCalorParameters();
#endif
  return worldPhysical;
}

void PhotInDetectorConstruction::PrintCalorParameters() const
{
  G4cout << "----PhotInDetectorConstruction::PrintCalorParameters()-----" << G4endl;
  if(serial) G4cout << " Calorimeters are placed in serial." << G4endl;
  else       G4cout << " Calorimeters are placed in parallel." << G4endl;
  G4cout << " Absorber is made of " << absorberMaterial->GetName() << G4endl;
  G4cout << " Gap is made of " << gapMaterial->GetName() << G4endl;
  G4cout << "--------------------------------------------------------" << G4endl;
}

void PhotInDetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
#ifdef debug
  G4cout<<"PhotInDetectorConstruction::SetAbsorberMaterial: "<<materialChoice<<G4endl;
#endif
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if(pttoMaterial)
  {
    absorberMaterial = pttoMaterial;
    layerParam->SetAbsorberMaterial(pttoMaterial);
    for(G4int i=0; i<PhotInNumSections; i++)
    {
      calorLogical[i]->SetMaterial(absorberMaterial);
      layerLogical[i]->SetMaterial(absorberMaterial);
    }
  }
  else G4cerr<<"PhotInDetectorConst::SetAbsM:"<<materialChoice<<" is not defined."<<G4endl;
}

void PhotInDetectorConstruction::SetGapMaterial(G4String materialChoice)
{
#ifdef debug
  G4cout<<"PhotInDetectorConstruction::SetGapMaterial: "<<materialChoice<<G4endl;
#endif
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);  
  if(pttoMaterial)
  {
    gapMaterial = pttoMaterial;
    gapParam->SetGapMaterial(pttoMaterial);
  }
  else G4cerr<<"PhotInDetectorConst::SetGapM:"<<materialChoice<<" is not defined."<<G4endl;
}

void PhotInDetectorConstruction::SetSerialGeometry(G4bool ser)
{
#ifdef debug
  G4cout<<"PhotInDetectorConstruction::SetSerialGeometry: "<<ser<<G4endl;
#endif
  if(serial==ser) return;         // Do nothing if serialization is the same
  serial=ser;
  G4double shiftZ=zHD+zHD;        // shift perpendicular to layers
  G4double shiftY=yHD+yHD;        // shift perpendicular to slabs
  for(G4int i=0; i<PhotInNumSections; i++)
  {
    if(serial) calorPhysical[i]->SetTranslation(G4ThreeVector(0.,0.,(i-1)*shiftZ));
    else calorPhysical[i]->SetTranslation(G4ThreeVector(0.,(i-1)*shiftY,0.));
  }
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void PhotInDetectorConstruction::SetNumberOfLayers(G4int nl)
{
#ifdef debug
  G4cout<<"PhotInDetectorConstruction::SetNumberOfLayers: "<<nl<<G4endl;
#endif
  numberOfLayers = nl;
  layerSolid->SetZHalfLength(zHD/numberOfLayers);
  layerParam->SetNumberOfLayers(nl);
  for(G4int i=0; i<PhotInNumSections; i++)
  { 
    if(layerPhysical[i])
    {
      if(slabPhysical[i])
      {
        layerLogical[i]->RemoveDaughter(slabPhysical[i]);
        delete slabPhysical[i];
      }
      calorLogical[i]->RemoveDaughter(layerPhysical[i]);
      delete layerPhysical[i];
    }
    layerPhysical[i] =  new G4PVParameterised("Layer", layerLogical[i], calorLogical[i],
                                              kZAxis , numberOfLayers , layerParam);
    slabPhysical[i] =  new G4PVParameterised("Slab", slabLogical[i], layerLogical[i],
                                              kYAxis , numberOfSlabs , gapParam);
  }
  PhotInCalorimeterSD::SetNumberOfLayers(nl);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void PhotInDetectorConstruction::SetNumberOfSlabs(G4int sn)
{
#ifdef debug
  G4cout<<"PhotInDetectorConstruction::SetNumberOfSlabs: "<<sn<<G4endl;
#endif
  numberOfSlabs = sn;
  slabSolid->SetZHalfLength(yHD/sn);
  gapParam->SetNumberOfSlabs(sn);
  for(G4int i=0; i<PhotInNumSections; i++)
  { 
    if(slabPhysical[i])
    {
      layerLogical[i]->RemoveDaughter(slabPhysical[i]);
      delete slabPhysical[i];
    }
    slabPhysical[i] =  new G4PVParameterised("Slab", slabLogical[i], layerLogical[i],
                                              kYAxis , numberOfSlabs , gapParam);
  }
  PhotInCalorimeterSD::SetNumberOfSlabs(sn);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void PhotInDetectorConstruction::CreateMaterial(G4String materialChoice)
{
#ifdef debug
  G4cout<<"PhotInDetectorConstruction::CreateMaterial: "<<materialChoice<<G4endl;
#endif
  if(G4Material::GetMaterial(materialChoice) != 0) return;
   
  G4double a, z, density;  
  // List of possible new materials (Material Factory)
  if (materialChoice == "Silicon")
    new G4Material("Silicon", z=14., a= 28.09*g/mole, density= 2.33*g/cm3);
  else if  (materialChoice =="Iron")
    new G4Material("Iron", z=26., a=55.85*g/mole, density=7.87*g/cm3);
  else if  (materialChoice =="ArgonGas")
    new G4Material("ArgonGas",z=18., a= 39.95*g/mole, density=1.782*mg/cm3);
  else if  (materialChoice =="He")
    new G4Material("He", z=2., a=4.0*g/mole, density=0.1786e-03*g/cm3);
  else G4cerr<<"**PhotInDetectorConstruction::CreateMaterial: No "<<materialChoice<<G4endl;
}
