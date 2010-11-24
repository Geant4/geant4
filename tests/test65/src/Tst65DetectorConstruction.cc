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
// 070522 Add Nist Materals T. Koi
// 080901 Add TK Materals T. Koi
//

#include "Tst65DetectorConstruction.hh"

#include "Tst65DetectorMessenger.hh"


#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

Tst65DetectorConstruction::Tst65DetectorConstruction()
  : simpleBoxLog(0),selectedMaterial(0),Air(0),Al(0),Pb(0),U(0),
    elN(0),elO(0),elU(0),isoU(0)
{
  detectorMessenger = new Tst65DetectorMessenger(this);
  materialChoice = "Pb";
  defineNISTMaterials();
  defineTKMaterials();
}

Tst65DetectorConstruction::~Tst65DetectorConstruction()
{
  delete detectorMessenger;
  delete Pb;  delete Al;  delete Air;  delete U;
  delete isoU;  delete elU;  delete elO;  delete elN;
}

void Tst65DetectorConstruction::SelectMaterial(G4String val)
{
  materialChoice = val;
  SelectMaterialPointer();
  G4cout << "SimpleBox is now made of " << materialChoice << G4endl;
}

void Tst65DetectorConstruction::SelectMaterialPointer()
{
//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  if(!Air)
  {
    a = 14.01*g/mole;
    elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
    a = 16.00*g/mole;
    elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);
    density = 1.29e-03*g/cm3;
    Air = new G4Material(name="Air", density, nel=2);
    Air->AddElement(elN, .7);
    Air->AddElement(elO, .3);
  }

  if(!Al)
  {
    a = 26.98*g/mole;
    density = 2.7*g/cm3;
    Al = new G4Material(name="Aluminium", z=13., a, density);
  }

  if(!Pb)
  {
    a = 207.19*g/mole;
    density = 11.35*g/cm3;
    Pb = new G4Material(name="Lead", z=82., a, density);
  }

  if(!U)
  {
    a = 238.*g/mole;
    density = 21.*g/cm3;
    U  = new G4Material(name="Uranium", density, 1);
    elU = new G4Element(name="Uranium", "U138", 1);
    isoU = new G4Isotope("U238", 92, 238, a);
    elU->AddIsotope(isoU, 100.0);
    U->AddElement(elU, 1);
  }


  if(materialChoice=="Air")
  { selectedMaterial = Air; }
  else if(materialChoice=="Al")
  { selectedMaterial = Al; }
  else if(materialChoice=="Pb")
  { selectedMaterial = Pb; }
  else
  { 
    SetMaterial( materialChoice ); 
    if ( selectedMaterial == 0 ) { selectedMaterial = U; }
  }

  if(simpleBoxLog)
  { simpleBoxLog->SetMaterial(selectedMaterial); }
}

G4VPhysicalVolume* Tst65DetectorConstruction::Construct()
{
  SelectMaterialPointer();

  G4Box * mySimpleBox = new G4Box("SBox",5*cm, 5*cm, 20*m);
  simpleBoxLog = new G4LogicalVolume( mySimpleBox,
                                      selectedMaterial,"SLog",0,0,0);
  G4VPhysicalVolume* simpleBoxDetector = new G4PVPlacement(0,G4ThreeVector(),
                                        "SPhys",simpleBoxLog,0,false,0);

  return simpleBoxDetector;
}

#include "G4NistManager.hh"
void Tst65DetectorConstruction::defineNISTMaterials()
{
//--------- Material definition with NIST Manager ---------

   G4NistManager* nistMan = G4NistManager::Instance();

   nistMan->FindOrBuildMaterial( "G4_H" );
   nistMan->FindOrBuildMaterial( "G4_He" );
   nistMan->FindOrBuildMaterial( "G4_Li" );
   nistMan->FindOrBuildMaterial( "G4_Be" );
   nistMan->FindOrBuildMaterial( "G4_B" );
   nistMan->FindOrBuildMaterial( "G4_C" );
   nistMan->FindOrBuildMaterial( "G4_N" );
   nistMan->FindOrBuildMaterial( "G4_O" );
/*
   nistMan->FindOrBuildMaterial( "G4_F" );
   nistMan->FindOrBuildMaterial( "G4_Ne" );
   nistMan->FindOrBuildMaterial( "G4_Na" );
   nistMan->FindOrBuildMaterial( "G4_Mg" );
*/
   nistMan->FindOrBuildMaterial( "G4_Al" );
/*
   nistMan->FindOrBuildMaterial( "G4_Si" );
   nistMan->FindOrBuildMaterial( "G4_P" );
   nistMan->FindOrBuildMaterial( "G4_S" );
   nistMan->FindOrBuildMaterial( "G4_Cl" );
   nistMan->FindOrBuildMaterial( "G4_Ar" );
   nistMan->FindOrBuildMaterial( "G4_K" );
   nistMan->FindOrBuildMaterial( "G4_Ca" );
   nistMan->FindOrBuildMaterial( "G4_Sc" );
   nistMan->FindOrBuildMaterial( "G4_Ti" );
   nistMan->FindOrBuildMaterial( "G4_V" );
   nistMan->FindOrBuildMaterial( "G4_Cr" );
   nistMan->FindOrBuildMaterial( "G4_Mn" );
*/
   nistMan->FindOrBuildMaterial( "G4_Fe" );
/*
   nistMan->FindOrBuildMaterial( "G4_Co" );
   nistMan->FindOrBuildMaterial( "G4_Ni" );
*/
   nistMan->FindOrBuildMaterial( "G4_Cu" );
   nistMan->FindOrBuildMaterial( "G4_Zn" );
/*
   nistMan->FindOrBuildMaterial( "G4_Ga" );
*/
   nistMan->FindOrBuildMaterial( "G4_Ge" );
/*
   nistMan->FindOrBuildMaterial( "G4_As" );
   nistMan->FindOrBuildMaterial( "G4_Se" );
   nistMan->FindOrBuildMaterial( "G4_Br" );
   nistMan->FindOrBuildMaterial( "G4_Kr" );
   nistMan->FindOrBuildMaterial( "G4_Rb" );
   nistMan->FindOrBuildMaterial( "G4_Sr" );
   nistMan->FindOrBuildMaterial( "G4_Y" );
   nistMan->FindOrBuildMaterial( "G4_Zr" );
   nistMan->FindOrBuildMaterial( "G4_Nb" );
   nistMan->FindOrBuildMaterial( "G4_Mo" );
   nistMan->FindOrBuildMaterial( "G4_Tc" );
   nistMan->FindOrBuildMaterial( "G4_Ru" );
   nistMan->FindOrBuildMaterial( "G4_Rh" );
   nistMan->FindOrBuildMaterial( "G4_Pd" );
*/
   nistMan->FindOrBuildMaterial( "G4_Ag" );
/*
   nistMan->FindOrBuildMaterial( "G4_Cd" );
   nistMan->FindOrBuildMaterial( "G4_In" );
   nistMan->FindOrBuildMaterial( "G4_Sn" );
   nistMan->FindOrBuildMaterial( "G4_Sb" );
   nistMan->FindOrBuildMaterial( "G4_Te" );
   nistMan->FindOrBuildMaterial( "G4_I" );
   nistMan->FindOrBuildMaterial( "G4_Xe" );
   nistMan->FindOrBuildMaterial( "G4_Cs" );
   nistMan->FindOrBuildMaterial( "G4_Ba" );
   nistMan->FindOrBuildMaterial( "G4_La" );
   nistMan->FindOrBuildMaterial( "G4_Ce" );
   nistMan->FindOrBuildMaterial( "G4_Pr" );
   nistMan->FindOrBuildMaterial( "G4_Nd" );
   nistMan->FindOrBuildMaterial( "G4_Pm" );
   nistMan->FindOrBuildMaterial( "G4_Sm" );
   nistMan->FindOrBuildMaterial( "G4_Eu" );
   nistMan->FindOrBuildMaterial( "G4_Gd" );
   nistMan->FindOrBuildMaterial( "G4_Tb" );
   nistMan->FindOrBuildMaterial( "G4_Dy" );
   nistMan->FindOrBuildMaterial( "G4_Ho" );
   nistMan->FindOrBuildMaterial( "G4_Er" );
   nistMan->FindOrBuildMaterial( "G4_Tm" );
   nistMan->FindOrBuildMaterial( "G4_Yb" );
   nistMan->FindOrBuildMaterial( "G4_Lu" );
   nistMan->FindOrBuildMaterial( "G4_Hf" );
   nistMan->FindOrBuildMaterial( "G4_Ta" );
   nistMan->FindOrBuildMaterial( "G4_W" );
   nistMan->FindOrBuildMaterial( "G4_Re" );
   nistMan->FindOrBuildMaterial( "G4_Os" );
   nistMan->FindOrBuildMaterial( "G4_Ir" );
   nistMan->FindOrBuildMaterial( "G4_Pt" );
*/
   nistMan->FindOrBuildMaterial( "G4_Au" );
/*
   //nistMan->FindOrBuildMaterial( "G4_Hg" ); // No data HP
   nistMan->FindOrBuildMaterial( "G4_Tl" );
   nistMan->FindOrBuildMaterial( "G4_Pb" );
   nistMan->FindOrBuildMaterial( "G4_Bi" );
   nistMan->FindOrBuildMaterial( "G4_Po" );
   nistMan->FindOrBuildMaterial( "G4_At" );
   //nistMan->FindOrBuildMaterial( "G4_Rn" ); // No data HP
   //nistMan->FindOrBuildMaterial( "G4_Fr" ); // No data HP
   nistMan->FindOrBuildMaterial( "G4_Ra" );
   nistMan->FindOrBuildMaterial( "G4_Ac" );
   nistMan->FindOrBuildMaterial( "G4_Th" );
   nistMan->FindOrBuildMaterial( "G4_Pa" );
*/
   nistMan->FindOrBuildMaterial( "G4_U" );
   //nistMan->FindOrBuildMaterial( "G4_Np" );
   //nistMan->FindOrBuildMaterial( "G4_Pu" );
   //nistMan->FindOrBuildMaterial( "G4_Am" );
   //nistMan->FindOrBuildMaterial( "G4_Cm" );
   //nistMan->FindOrBuildMaterial( "G4_Bk" );
   //nistMan->FindOrBuildMaterial( "G4_Cf" ); // No data HP even with Am

}

void Tst65DetectorConstruction::defineTKMaterials()
{

   G4NistManager* nistMan = G4NistManager::Instance();

   G4String symbol;            
   G4double a, density;         //a=mass of a mole;
   G4int iz, n;                 //iz=number of protons  in an isotope;
   G4int ncomponents, natoms;
   G4double abundance, fractionmass;

   G4Isotope* isoH1 = new G4Isotope("H1", iz=1, n=1, a=1.0078250321*g/mole);
   G4Element* elH1  = new G4Element("Hydrogen1",symbol="H1" , ncomponents=1);
   elH1->AddIsotope(isoH1, abundance=100.*perCent);
   G4Material* Hydrogen1 = new G4Material("Hydrogen1" , density= 1.0*g/cm3 , ncomponents = 1 , kStateSolid , 0.*kelvin , 1*atmosphere );
   Hydrogen1->AddElement( elH1 , fractionmass = 1 );

   // Create Element for Thermal Scattering
   G4Element* elTSHW= new G4Element( "TS_H_of_Water" , "H_WATER" , 1.0 , 1.0079*g/mole );
   G4Element* elTSH= new G4Element( "TS_H_of_Polyethylene" , "H_POLYETHYLENE" , 1.0 , 1.0079*g/mole );
   // Create Materials from the elements
   G4Material* matH2O_TS = new G4Material( "Water_TS" , density = 1.0*g/cm3 , ncomponents= 2 );
   matH2O_TS -> AddElement(elTSHW,natoms=2);
   matH2O_TS -> AddElement(nistMan->FindOrBuildElement( "O" ), natoms = 1);
   G4Material* matCH2_TS = new G4Material( "Polyethylene_TS" , density = 0.94*g/cm3 , ncomponents= 2 );
   matCH2_TS -> AddElement(elTSH,natoms=2);
   matCH2_TS -> AddElement(nistMan->FindOrBuildElement( "C" ), natoms= 1);

}

void Tst65DetectorConstruction::SetMaterial( G4String materialChoice )
{
  // search the material by its name
  selectedMaterial = G4Material::GetMaterial( materialChoice );
}

