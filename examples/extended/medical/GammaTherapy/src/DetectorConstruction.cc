//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------
//      GEANT4 ibrem test
//
// Authors: V.Grichine, V.Ivanchenko
//
// Modified:
//
// 18-02-03 V.Ivanchenko create
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "PhantomSD.hh"
#include "TargetSD.hh"
#include "CheckVolumeSD.hh"
#include "Histo.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"
#include "PhantomSD.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"


#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "globals.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction()
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::Initialise()
{
  logicTarget1 = 0;
  logicTarget2 = 0;
  DefineMaterials();
  dMessenger = new DetectorMessenger(this);

  checkSD = new CheckVolumeSD("checkSD");
  (G4SDManager::GetSDMpointer())->AddNewDetector( checkSD );
  calorimeterSD = new PhantomSD("phantomSD");
  (G4SDManager::GetSDMpointer())->AddNewDetector( calorimeterSD );
  targetSD = new TargetSD("targetSD");
  (G4SDManager::GetSDMpointer())->AddNewDetector( targetSD );

  fDistanceVacuumTarget = 30.*mm,

  fDelta                = 0.001*mm;

  fTargetRadius         = 100.*mm;
  fTarget1Z             = 9.*mm;
  fTarget2Z             = 6.*mm;

  fGasVolumeRadius      = 210.*mm;
  fGasVolumeZ           = 690.*mm;
  fMylarVolumeZ         = 0.02*mm;

  fCheckVolumeZ         = 0.1*mm;
  fCheckShiftZ          = 200.*mm;

  fAbsorberRadius       = 200.*mm;
  fPhantomRadius        = 300.*mm;
  fPhantomZ             = 300.*mm;

  fAirZ                 = 210.*mm;
  fAbsorberShiftZ       = 70.*mm;
  fWindowZ              = 0.05*mm;
}


///////////////////////////////////////////////////////////////

void DetectorConstruction::InitialiseGeometryParameters()
{
  // Volumee sizes

  G4double factor = 1.2;

  fWorldXY       = factor*std::max(fPhantomRadius,fGasVolumeRadius);
  G4double nz    = (G4int)((Histo::GetPointer())->GetNumberDivZ());
  fAbsorberZ     = fPhantomZ/nz;
  fGasVolumeZ    = 1000.*mm - fAbsorberShiftZ - fAirZ - fTarget1Z - fTarget2Z;

  G4double ztot  = fGasVolumeZ + fAirZ + fPhantomZ + fDistanceVacuumTarget;
  fTargetVolumeZ = fDistanceVacuumTarget + fTarget2Z + fTarget1Z + fDelta;
  fWorldZ  = factor*ztot*0.5;

  if(fCheckShiftZ < fDelta) fCheckShiftZ = fDelta;
  if(fCheckShiftZ > fAirZ - fCheckVolumeZ -fDelta)
      fCheckShiftZ = fAirZ - fCheckVolumeZ -fDelta;

  // Z position of volumes from upstream to downstream

  fWindowPosZ      =  -(ztot + fWindowZ)*0.5;
  fGeneratorPosZ   =  fWindowPosZ - 0.5*fWindowZ - fDelta;

  fTargetVolumePosZ=  -0.5*(ztot - fTargetVolumeZ);
  fTarget1PosZ     =  -0.5*(fTargetVolumeZ - fTarget1Z) + fDistanceVacuumTarget;
  fTarget2PosZ     =  fTarget1PosZ + 0.5*(fTarget2Z + fTarget1Z);

  fGasVolumePosZ   =  fTargetVolumePosZ + 0.5*(fTargetVolumeZ + fGasVolumeZ);
  fCheckVolumePosZ =  fGasVolumePosZ + 0.5*(fGasVolumeZ + fCheckVolumeZ)
                                +  fCheckShiftZ;
  fMylarPosZ       =  fGasVolumePosZ + 0.5*(fGasVolumeZ + fMylarVolumeZ) + fDelta;

  fPhantomPosZ     =  fGasVolumePosZ + 0.5*(fGasVolumeZ + fPhantomZ) + fAirZ;
  fAbsorberPosZ    =  fAbsorberShiftZ - 0.5*(fPhantomZ - fAbsorberZ);
  (Histo::GetPointer())->SetAbsorberZ(fPhantomZ);
  (Histo::GetPointer())->SetAbsorberR(fAbsorberRadius);
  (Histo::GetPointer())->SetScoreZ(fAbsorberShiftZ);
  G4double shiftZPh = fPhantomPosZ-0.5*fPhantomZ;
  calorimeterSD->setShiftZ(shiftZPh);
  G4cout << "===================================================" << G4endl;
  G4cout << "#                  IBREM geometry                 #" << G4endl;
  G4cout << "===================================================" << G4endl;
  G4cout << "  World   width= " << fWorldZ/mm << " mm " << G4endl;
  G4cout << "  Window  width= " << fWindowZ/mm << " mm       position = "
                                << fWindowPosZ/mm << " mm:" << G4endl;
  G4cout << "  TargetV width= " << fTargetVolumeZ/mm << " mm       position = "
                                << fTargetVolumePosZ/mm << " mm:" << G4endl;
  G4cout << "  Target1 width= " << fTarget1Z/mm << " mm       position = "
                                << fTarget1PosZ/mm << " mm:" << G4endl;
  G4cout << "  Target2 width= " << fTarget2Z/mm << " mm       position = "
                                << fTarget2PosZ/mm << " mm:" << G4endl;
  G4cout << "  Gas     width= " << fGasVolumeZ/mm << " mm     position = "
                                << fGasVolumePosZ/mm << " mm:" << G4endl;
  G4cout << "  Mylar   width= " << fMylarVolumeZ/mm << " mm    position = "
                                << fMylarPosZ/mm << " mm:" << G4endl;
  G4cout << "  Check   width= " << fCheckVolumeZ/mm << " mm     position = "
                                << fCheckVolumePosZ/mm << " mm:" << G4endl;
  G4cout << "  Air     width= " << fAirZ/mm << " mm " << G4endl;
  G4cout << "  Phantom width= " << fPhantomZ/mm << " mm     position = "
                                << fPhantomPosZ/mm << " mm:" << G4endl;
  G4cout << "  Absorb  width= " << fAbsorberZ/mm << " mm       position = "
                                << fAbsorberPosZ/mm << " mm:" << G4endl;
  G4cout << "  Absorb  shift= " << shiftZPh/mm << " mm " << G4endl;
  G4cout << "  Target1        " << fTarget1Material->GetName() << G4endl;
  G4cout << "  Target2        " << fTarget2Material->GetName() << G4endl;
  G4cout << "  Phantom        " << fAbsorberMaterial->GetName() << G4endl;
  G4cout << "===================================================" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::DefineMaterials()
{

  G4String name, symbol;             //a=mass of a mole;
  G4double a, z, density;            //z=mean number of protons;

  G4int    iz,n,ncomponents, natoms;
  G4double fractionmass;
  G4double temperature, pressure;

  //  std::vector<G4Material*> list;
  G4Material* ma = 0;

//
// define Elements
//

  //  a = 1.01*g/mole;
  // G4Element* elH  = new G4Element(name="Hydrogen",symbol="H", z= 1., a);

  a = 1.01*g/mole;
  G4Isotope* ih1 = new G4Isotope("Hydrogen",iz=1,n=1,a);

  a = 2.01*g/mole;
  G4Isotope* ih2 = new G4Isotope("Deuterium",iz=1,n=2,a);

  G4Element* elH = new G4Element(name="Hydrogen",symbol="H",2);
  elH->AddIsotope(ih1,.999);
  elH->AddIsotope(ih2,.001);


  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  a = 12.00*g/mole;
  G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  a = 39.948*g/mole;
  G4Element* elAr = new G4Element(name="Argon", symbol="Ar", z=18., a);

   a = 69.723*g/mole;
  G4Element* elGa  = new G4Element(name="Gallium"  ,symbol="Ga" , z= 31., a);

  a = 74.9216*g/mole;
  G4Element* elAs  = new G4Element(name="Arsenicum"  ,symbol="As" , z= 33., a);

  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);

  G4Element*   I  = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);


//
// define simple materials
//
  density = 1.848*g/cm3;
  a = 9.01*g/mole;
  ma = new G4Material(name="Be", z=4., a, density);
  fTarget1Material = ma;
  fWindowMaterial  = ma;

  density = 19.35*g/cm3;
  a = 183.85*g/mole;
  ma = new G4Material(name="W", z=74., a, density);
  fTarget2Material = ma;

  density = 8.960*g/cm3;
  a = 63.55*g/mole;
  ma = new G4Material(name="Cu"   , z=29., a, density);

  // Helium as light gas volume

  density = 0.178*mg/cm3 ;
  a = 4.0026*g/mole ;
  G4Material* He  = new G4Material(name="He",z=2., a, density );
  fLightMaterial = He;

  density = 2.699*g/cm3;
  a = 26.98*g/mole;
  ma = new G4Material(name="Aluminum", z=13., a, density);

  density = 2.265*g/cm3;
  a = 12.0107*g/mole;
  ma = new G4Material(name="Carbon", z=6., a, density);

  density = 2.330*g/cm3;
  a = 28.09*g/mole;
  ma = new G4Material(name="Silicon", z=14., a, density);

  density = 1.390*g/cm3;
  a = 39.95*g/mole;
  ma = new G4Material(name="LiquidArgon", z=18., a, density);

  density = 3.02*g/cm3;
  a = 131.29*g/mole;
  ma = new G4Material(name="LiquidXenon", z=54., a, density);

  density = 7.874*g/cm3;
  a = 55.85*g/mole;
  ma = new G4Material(name="Iron"   , z=26., a, density);

  density = 5.323*g/cm3;
  a = 72.61*g/mole;
  ma = new G4Material(name="Germanium", z=32., a, density);

  density = 19.32*g/cm3;
  a =196.97*g/mole;
  ma = new G4Material(name="Gold"   , z=79., a, density);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  ma = new G4Material(name="Lead"     , z=82., a, density);

//
// define a material from elements.   case 1: chemical molecule
//


  density = 1.000*g/cm3;
  ma = new G4Material("Water", density, 2);
  ma->SetChemicalFormula("H_2O");
  ma->AddElement(elH, natoms=2);
  ma->AddElement(elO, natoms=1);
  G4double exc = ma->GetIonisation()->FindMeanExcitationEnergy("H_2O");
  ma->GetIonisation()->SetMeanExcitationEnergy(exc);
  fAbsorberMaterial = ma;

  density = 0.0006672*g/cm3;
  ma = new G4Material("Methane", density, 2);
  ma->SetChemicalFormula("CH_4");
  ma->AddElement(elH, natoms=4);
  ma->AddElement(elC, natoms=1);

  ma = new G4Material("Graphite", 2.265*g/cm3, 1);
  ma->SetChemicalFormula("Graphite");
  ma->AddElement( elC, 1 );

  density = 5.3176*g/cm3;
  ma = new G4Material("GaAs", density, ncomponents=2);
  ma->SetChemicalFormula("GaAS");
  ma->AddElement(elGa, natoms=1);
  ma->AddElement(elAs, natoms=1);

  ma = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  ma->SetChemicalFormula("C_2H_6");
  ma->AddElement(elH,6);
  ma->AddElement(elC,2);

  ma = new G4Material ("CsI" , 4.51*g/cm3, 2);
  ma->SetChemicalFormula("CsI");
  ma->AddElement(Cs,1);
  ma->AddElement(I,1);

//
// define a material from elements.   case 2: mixture by fractional mass
//

// Dry air (average composition)

  density = 1.25053*mg/cm3 ;       // STP
  G4Material* Nitrogen = new G4Material(name="N2"  , density, ncomponents=1);
  Nitrogen->AddElement(elN, 2);

  density = 1.4289*mg/cm3 ;       // STP
  G4Material* Oxygen = new G4Material(name="O2"  , density, ncomponents=1);
  Oxygen->AddElement(elO, 2);

  density = 1.7836*mg/cm3 ;       // STP
  G4Material* Argon = new G4Material(name="Argon"  , density, ncomponents=1);
  Argon->AddElement(elAr, 1);

  density = 1.2928*mg/cm3 ;       // STP
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=3);
  Air->AddMaterial( Nitrogen, fractionmass = 0.7557 ) ;
  Air->AddMaterial( Oxygen,   fractionmass = 0.2315 ) ;
  Air->AddMaterial( Argon,    fractionmass = 0.0128 ) ;

  fWorldMaterial = Air;

  density = 1.39*g/cm3;
  ma = new G4Material("Mylar"  , density, ncomponents=3);
  ma->AddElement(elC, natoms=10);
  ma->AddElement(elH, natoms=18);
  ma->AddElement(elO, natoms=5);
  fMylar = ma;

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  a = 1.01*g/mole;
  z = 1.0;
  ma = new G4Material("Vacuum", z, a, density,
                                      kStateGas,temperature,pressure);
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
   // Cleanup old geometry
   G4PhysicalVolumeStore::GetInstance()->Clean();
   G4LogicalVolumeStore::GetInstance()->Clean();
   G4SolidStore::GetInstance()->Clean();

   //
   InitialiseGeometryParameters();

   //
   // World
   //

  G4VPhysicalVolume* pv;

  G4Box* solidWorld = new G4Box("World",fWorldXY,fWorldXY,fWorldZ);

  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,
                                                               fWorldMaterial,"World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,G4ThreeVector(),"World",
                                                               logicWorld,0,false,0);

  // Be Vacuum window
  G4Tubs* solidWin = new G4Tubs("Window",0.,fTargetRadius*0.25,0.5*fWindowZ,0.,twopi);
  G4LogicalVolume* logicWin = new G4LogicalVolume(solidWin,fWindowMaterial,"Window");
  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,fWindowPosZ),"Window",logicWin,
                                                               physWorld,false,0);

  // Target Volume
  G4Tubs* solidTGVolume = new G4Tubs("TargetVolume",0.,fTargetRadius,
                                                               0.5*fTargetVolumeZ,0.,twopi);
  G4LogicalVolume* logicTGVolume = new G4LogicalVolume(solidTGVolume,
                                                               fLightMaterial,"TargetVolume");
  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,fTargetVolumePosZ),logicTGVolume,"TargetVolume",
                                                               logicWorld,false,0);

  // Target 1
  G4Tubs* solidTarget1 = new G4Tubs("Target1",0.,fTargetRadius*0.5,0.5*fTarget1Z,0.,twopi);
  logicTarget1 = new G4LogicalVolume(solidTarget1,fTarget1Material,"Target1");
  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,fTarget1PosZ),logicTarget1,"Target1",
                                                               logicTGVolume,false,0);
  (Histo::GetPointer())->SetTarget1(pv);
  logicTarget1->SetSensitiveDetector(targetSD);


  // Target 2 (for combined targets)
  G4Tubs* solidTarget2 = new G4Tubs("Target2",0.,fTargetRadius*0.5,0.5*fTarget2Z,0.,twopi);
  logicTarget2 = new G4LogicalVolume(solidTarget2,fTarget2Material,"Target2");
  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,fTarget2PosZ),logicTarget2,"Target2",
                                                               logicTGVolume,false,0);

  (Histo::GetPointer())->SetTarget2(pv);
  logicTarget2->SetSensitiveDetector(targetSD);

  // Gas Volume

  G4Tubs* solidGasVolume = new G4Tubs("GasVolume",0.,fGasVolumeRadius,
                                                               0.5*fGasVolumeZ,0.,twopi);
  G4LogicalVolume* logicGasVolume = new G4LogicalVolume(solidGasVolume,
                                                               fLightMaterial,"GasVolume");
  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,fGasVolumePosZ),
                                                               "GasVolume",logicGasVolume,
                                                               physWorld,false,0);
  (Histo::GetPointer())->SetGasVolume(pv);

  // Mylar window

  G4Tubs* sMylarVolume = new G4Tubs("Mylar",0.,fGasVolumeRadius,
                                                               0.5*fMylarVolumeZ,0.,twopi);
  G4LogicalVolume* lMylarVolume = new G4LogicalVolume(sMylarVolume,
                                                               fMylar,"Mylar");
  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,fMylarPosZ),"Mylar",lMylarVolume,
                                                               physWorld,false,0);


  // Check Volume

  G4Tubs* solidCheckVolume = new G4Tubs("CheckVolume",0.,fGasVolumeRadius,
                                                               0.5*fCheckVolumeZ,0.,twopi);
  G4LogicalVolume* logicCheckVolume = new G4LogicalVolume(solidCheckVolume,
                                                              fWorldMaterial,"CheckVolume");
  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,fCheckVolumePosZ),
                                                              "CheckVolume",logicCheckVolume,
                                                              physWorld,false,0);
  (Histo::GetPointer())->SetCheckVolume(pv);
  logicCheckVolume->SetSensitiveDetector(checkSD);

  // Phantom

  G4Box* solidPhantom = new G4Box("Phantom",fPhantomRadius,fPhantomRadius,
                                                              0.5*fPhantomZ);
  G4LogicalVolume* logicPhantom = new G4LogicalVolume(solidPhantom,
                                                              fAbsorberMaterial,"Phantom");
  G4VPhysicalVolume* physPhantom = new G4PVPlacement(0,
                                                              G4ThreeVector(0.,0.,fPhantomPosZ),
                                                              "Phantom",logicPhantom,
                                                              physWorld,false,0);

  G4Tubs* solidPh = new G4Tubs("PhantomSD",0.,fAbsorberRadius,0.5*fPhantomZ,0.,twopi);
  G4LogicalVolume* logicPh = new G4LogicalVolume(solidPh,fAbsorberMaterial,"PhantomSD");
  G4VPhysicalVolume* physPh = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
                                                              "Phantom",logicPh,
                                                              physPhantom,false,0);
  G4cout << "Phantom R= " << fAbsorberRadius << " dz= " << 0.5*fPhantomZ << G4endl;

  // Sensitive Absorber

  G4double absWidth = 0.5*fAbsorberZ;
  G4Tubs* solidAbsorber = new G4Tubs("Absorber",0.,fAbsorberRadius,absWidth,0.,twopi);
  G4LogicalVolume* logicAbsorber = new G4LogicalVolume(solidAbsorber,
                                                              fAbsorberMaterial,"Absorber");
  G4cout << "Absorber R= " << fAbsorberRadius << " dz= " << absWidth << " posZ= " << fAbsorberPosZ<< G4endl;

  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,fAbsorberPosZ),"Absorber",logicAbsorber,
                                                              physPh,false,0);
  (Histo::GetPointer())->SetPhantom(physPh);
  G4int numR = (Histo::GetPointer())->GetNumberDivR();
  G4double stepR = fAbsorberRadius/(G4double)numR;

  G4double r1 = 0.0;
  G4double r2 = 0.0;
  G4Tubs* solidRing;
  G4LogicalVolume* logicRing;

  for(G4int k=0; k<numR; k++) {
    r2 = r1 + stepR;
    if(k == numR-1) r2 = fAbsorberRadius;
//    G4cout << "New ring r1= " << r1 << " r2= " << r2 << " dz= " << absWidth << G4endl;
    solidRing = new G4Tubs("Ring",r1,r2,absWidth,0.,twopi);
    logicRing = new G4LogicalVolume(solidRing,fAbsorberMaterial,"Ring");
    logicRing->SetSensitiveDetector(calorimeterSD);
    logicRing->SetVisAttributes(G4VisAttributes::Invisible);
    pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicRing,"Ring",
                                                              logicAbsorber,false,k);
    r1 = r2;
  }

  //
  // Sensitive Detectors: Absorber
  //

  logicPh->SetSensitiveDetector(calorimeterSD);
  logicAbsorber->SetSensitiveDetector(calorimeterSD);

  //
  // Visualization attributes
  //
  G4VisAttributes* VisAtt = 0;
  VisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  VisAtt->SetVisibility(true);
  logicAbsorber->SetVisAttributes(VisAtt);

  VisAtt= new G4VisAttributes(G4Colour(1.0,1.0,2.0));
  VisAtt->SetVisibility(true);
  logicPhantom->SetVisAttributes(VisAtt);

  VisAtt= new G4VisAttributes(G4Colour(1.0,0.0,2.0));
  VisAtt->SetVisibility(true);
  logicPh->SetVisAttributes(VisAtt);

  VisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  VisAtt->SetVisibility(true);
  logicAbsorber->SetVisAttributes(VisAtt);

  VisAtt= new G4VisAttributes(G4Colour(0.1,1.0,2.0));
  VisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(VisAtt);

  VisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  VisAtt->SetVisibility(true);
  logicGasVolume->SetVisAttributes(VisAtt);

  VisAtt= new G4VisAttributes(G4Colour(0.0,0.5,1.0));
  VisAtt->SetVisibility(true);
  logicTarget1->SetVisAttributes(VisAtt);
  logicTarget2->SetVisAttributes(VisAtt);
  logicTGVolume->SetVisAttributes(VisAtt);

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::setTarget1Material(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial)
     {
       fTarget1Material = pttoMaterial;
       G4cout << "New target1 material " << materialChoice << G4endl;
       if(logicTarget1) logicTarget1->SetMaterial(fTarget1Material);
     }
  else
       G4cout << "Material " << materialChoice << " is not found out!" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::setTarget2Material(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial)
     {
       fTarget2Material = pttoMaterial;
       G4cout << "New target2 material " << materialChoice << G4endl;
       if(logicTarget2) logicTarget2->SetMaterial(fTarget2Material);
     }
  else
       G4cout << "Material " << materialChoice << " is not found out!" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
