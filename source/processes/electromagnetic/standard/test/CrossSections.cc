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
// $Id: CrossSections.cc,v 1.2 2005-05-13 16:55:29 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// ------------------------------------------------------------
//
//  To print the cross sections
//
#include "G4Material.hh"
#include "G4PhotoElectricEffect52.hh"
#include "G4ComptonScattering52.hh"
#include "G4GammaConversion52.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4Gamma.hh"
#include "G4RegionStore.hh"

int main() {

  G4UnitDefinition::BuildUnitsTable();

  // define materials
  //
  G4double Z;

  new G4Material("Iodine", Z=53., 126.90*g/mole, 4.93*g/cm3);
  G4NistManager::Instance()->SetVerbose(1);

  // initialise processes
  // G4v52
  G4PhotoElectricEffect52* phot = new G4PhotoElectricEffect52();
  G4ComptonScattering52*   comp = new G4ComptonScattering52();
  G4GammaConversion52*     conv = new G4GammaConversion52();

  // Standard
  G4PhotoElectricEffect* phot = new G4PhotoElectricEffect();
  G4ComptonScattering*   comp = new G4ComptonScattering();
  G4GammaConversion*     conv = new G4GammaConversion();

  G4ParticleDefinition* gamma = G4Gamma::Gamma(); 
  G4RegionStore::GetInstance()->FindOrCreateRegion("dummy");

  comp->PreparePhysicsTable(*gamma);
  conv->PreparePhysicsTable(*gamma);
  phot->PreparePhysicsTable(*gamma); 

  // print cross section per atom
  //
  G4double Emin = 1.0*MeV, Emax = 2.01*MeV, dE = 100*keV;

  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
	   << "\tcomp= " << comp->ComputeCrossSectionPerAtom(Energy,Z)/barn
	   << "\tconv= " << conv->ComputeCrossSectionPerAtom(Energy,Z)/barn
	   << "\tphot= " << phot->ComputeCrossSectionPerAtom(Energy,Z)/barn 
	   << G4endl;
  }

  G4cout << G4endl;
                           
return EXIT_SUCCESS;
}
