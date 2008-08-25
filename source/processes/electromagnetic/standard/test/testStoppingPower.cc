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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// History: Test for calc. Stopping Power of ions
// 
//
// 11.08.08, A. Ivantchenko
// 

#include "globals.hh"
#include <iostream>
#include <vector>
//#include "G4MaterialStopping.hh"
//#include "G4SimpleMaterialStopping.hh"
#include "G4IronStopping.hh"
//#include "G4WaterStoppingRange.hh"

int main(){  
  
  /*  G4MaterialStopping mS;
  //  or G4MaterialStopping* mS = new G4MaterialStopping();
  G4double dedx = mS.GetDEDX(3,0,.025*MeV*6.941);
  G4String matN = mS.GetMaterialName(0);
  G4double dens = mS.GetDensity(0);
  G4cout << "DEDX for " << matN << " is "<< dedx*mm/MeV << " MeV/mm   Density= " 
	 << dens*cm3/gram << " g/cm3 " << G4endl;
  */ 
  /*  
  G4SimpleMaterialStopping mS;
  G4double dedx = mS.GetDEDX(3,0,.03*MeV*6.941);
  G4String matN = mS.GetMaterialName(0);
  G4double dens = mS.GetDensity(0);
  G4cout << "DEDX for " << matN << " is "<< dedx*mm/MeV << " MeV/mm   Density= " 
  << dens*cm3/gram << " g/cm3 " << G4endl;*/
  
  /* G4WaterStoppingRange mS1;
  G4double range = mS1.GetRange(18,.03*MeV*39.948);
  G4cout << "Range in Water is "<< range/mm << " mm  " << G4endl;*/
     
  G4IronStopping mS2;
  G4double dedx = mS2.GetDEDX(3,.025*MeV*55.847);
  G4String matN = mS2.GetMaterialName(0);
  G4double dens = mS2.GetDensity(0);
  G4cout << "DEDX for " << matN << " is "<< dedx*mm/MeV << " MeV/mm   Density= " 
  << dens*cm3/gram << " g/cm3 " << G4endl;
  return 0;
} 
 
