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
// $Id$
// ------------------------------------------------------------
// Geant4 class implementation file
//
// 03/09/2008, by T.Nikitina
// ------------------------------------------------------------

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

#include "AXPETMaterial.hh"

AXPETMaterial::AXPETMaterial ():
aluminum(0),matAir(0), LYSO(0), EJ280(0),Vacuum(0)
{;}

AXPETMaterial::~AXPETMaterial ()
{
  delete LYSO;
  delete EJ280;  
  delete matAir;  
  delete Vacuum;  
}

void AXPETMaterial::DefineMaterials()
{
  G4double a; //atomic mass
  G4double z; //atomic number
  G4double d; //density
  G4int nComponents; //number of elements in the material

  a=1.01*g/mole;
  G4Element* elH = new G4Element("Hydrogen", "H", z=1 , a);

  a = 14.01*g/mole;
  G4Element* elN = new G4Element("Nitrogen","N",z = 7.,a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O",z = 8.,a);

  a = 12.011*g/mole;
  G4Element* elC = new G4Element("Carbon","C",z = 6.,a);

  a = 28.09*g/mole;
  G4Element* elSi = new G4Element("Silicon", "Si", z = 14., a);

  a = 88.91*g/mole;
  G4Element* elY = new G4Element("Yttrium", "Y", z = 39., a);

  a = 174.97*g/mole;
  G4Element* elLu = new G4Element("Lutetium", "Lu", z = 71., a);


  //Aluminum
  aluminum =  new G4Material("Aluminium", z=13., a=26.98*g/mole, d=2.700*g/cm3);


  // AIR
  d = 1.290*mg/cm3;
  nComponents = 2;
  matAir = new G4Material("AIR", d, nComponents);
  matAir -> AddElement(elN,0.7);
  matAir -> AddElement(elO,0.3);

  // LYSO
  d = 7.1*g/cm3;
  nComponents=4;
  LYSO = new G4Material("LYSO", d, nComponents);
  LYSO->AddElement(elLu,0.72);
  LYSO->AddElement(elY, 0.04);
  LYSO->AddElement(elSi,0.06);
  LYSO->AddElement(elO, 0.18);

  // wls similar to VINYLTOLUENE
  d=1.02*g/cm3;
  nComponents = 2;
  EJ280 = new G4Material("WLSEJ280", d, nComponents);
  EJ280->AddElement(elH, 0.085);
  EJ280->AddElement(elC, 0.915);

  Vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole, d = universe_mean_density,
                           kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Adding scintillating material properties
  
  const G4int nEntries = 11;
  
  
 // G4double PhotonEnergy[nEntries] = {500, 490, 480, 470, 460, 450, 440, 430, 420, 410, 400};
  G4double PhotonEnergy[nEntries] = {2.478*eV, 2.53*eV, 2.58*eV, 2.636*eV, 2.69*eV, 2.75*eV, 2.816*eV, 2.88*eV, 2.95*eV, 3.022*eV, 3.097*eV};

  // Air properties
  G4double RIAir[nEntries] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}; // Is Isotropic
  G4MaterialPropertiesTable* MPTAir = new G4MaterialPropertiesTable();
  MPTAir->AddProperty("RINDEX", PhotonEnergy, RIAir, nEntries);
  matAir->SetMaterialPropertiesTable(MPTAir);
  Vacuum->SetMaterialPropertiesTable(MPTAir);

  // LYSO properties

  G4double RILyso[nEntries] = {1.82,1.82,1.82,1.82,1.82,1.82,1.82,1.82,1.82,1.82,1.82};// Is isotropic
  G4double ABLyso[nEntries] = {100.*cm, 100.*cm, 100.*cm, 100.*cm, 100.*cm, 100.*cm, 100.*cm, 80.* cm, 50. *cm, 30.*cm, 2.*cm}; // To be updated
  G4double SFLyso[nEntries] = {0., 0., 0., 0., 0., 0., 0.,0.1, 0.8, 0.1, 0.}; // To be updated
  
  G4MaterialPropertiesTable* MPTLyso = new G4MaterialPropertiesTable();
  MPTLyso->AddProperty("RINDEX",        PhotonEnergy, RILyso, nEntries);
  MPTLyso->AddProperty("ABSLENGTH",     PhotonEnergy, ABLyso, nEntries);
  MPTLyso->AddProperty("FASTCOMPONENT", PhotonEnergy, SFLyso, nEntries);
  MPTLyso->AddProperty("SLOWCOMPONENT", PhotonEnergy, SFLyso, nEntries);

  MPTLyso->AddConstProperty("SCINTILLATIONYIELD", 20000./MeV);
  MPTLyso->AddConstProperty("RESOLUTIONSCALE",    1.04);//R = resolutionscale*sqrt(yield)
  MPTLyso->AddConstProperty("FASTTIMECONSTANT", 40.*ns);
  MPTLyso->AddConstProperty("SLOWTIMECONSTANT",10000.*ns);
  MPTLyso->AddConstProperty("YIELDRATIO",1.0);

  LYSO->SetMaterialPropertiesTable(MPTLyso);


  // WLS properties


  G4double RIWls[nEntries] = {1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58}; 
  G4double ABWls[nEntries] = {200.*mm, 150.*mm, 60.*mm, 30.*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5 *mm}; // To be updated
  G4double EMWls[nEntries] = {0.15,0.6,0.15,0.1,0.,0.,0.,0.,0.,0.,0.}; 
  
  G4MaterialPropertiesTable* MPTWls = new G4MaterialPropertiesTable();
  MPTWls->AddProperty("RINDEX",               PhotonEnergy, RIWls, nEntries);
  MPTWls->AddProperty("WLSABSLENGTH",         PhotonEnergy, ABWls, nEntries);
  MPTWls->AddProperty("WLSCOMPONENT",         PhotonEnergy, EMWls, nEntries);
  MPTWls->AddConstProperty("WLSTIMECONSTANT", 8.5 *ns);

  EJ280->SetMaterialPropertiesTable(MPTWls);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

G4Material* AXPETMaterial::GetMat(const G4String& material)
{
  G4Material* pttomaterial = G4Material::GetMaterial(material);
  return pttomaterial;
}
