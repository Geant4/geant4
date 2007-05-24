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
// $Id: Materials.cc,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation - material definitions.
//
#include "Materials.hh"

#include "ConfigData.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"

namespace Materials {

  void Initialise()
  {
    G4NistManager* man = G4NistManager::Instance();
    
    new G4Material("User_Titanium", 22, man->FindOrBuildElement("Ti")->GetA(), 4.54*g/cm3);
    new G4Material("User_Silicon", 14,  man->FindOrBuildElement("Si")->GetA(), 2.33*g/cm3);
    new G4Material("User_Lead", 82,  man->FindOrBuildElement("Pb")->GetA(), 11.34*g/cm3);
    
    G4Material* air = new G4Material("User_Air", 1.205e-03*g/cm3, 4);
    air->AddElement(man->FindOrBuildElement("C"), 0.000124);
    air->AddElement(man->FindOrBuildElement("N"), 0.756);
    air->AddElement(man->FindOrBuildElement("O"), 0.232);
    air->AddElement(man->FindOrBuildElement("Ar"), 0.0128);
    
    G4Material* steel = new G4Material("User_Steel", 8.06*g/cm3, 6);
    steel->AddElement(man->FindOrBuildElement("C"), 0.001);
    steel->AddMaterial(G4Material::GetMaterial("User_Silicon"), 0.007);
    steel->AddElement(man->FindOrBuildElement("Cr"), 0.180);
    steel->AddElement(man->FindOrBuildElement("Mn"), 0.010);
    steel->AddElement(man->FindOrBuildElement("Fe"), 0.712);
    steel->AddElement(man->FindOrBuildElement("Ni"), 0.090);
    
    ConfigData::SetChamberWindowMaterial(steel);
    ConfigData::SetAirGap1Material(G4Material::GetMaterial("User_Air"));
    ConfigData::SetMonitorMaterial(G4Material::GetMaterial("User_Silicon"));
    ConfigData::SetAirGap2Material(G4Material::GetMaterial("User_Air"));
    ConfigData::SetBeamWindowMaterial(G4Material::GetMaterial("User_Titanium"));
    ConfigData::SetBeamPipeMaterial(man->FindOrBuildMaterial("G4_Galactic"));
  }
}
