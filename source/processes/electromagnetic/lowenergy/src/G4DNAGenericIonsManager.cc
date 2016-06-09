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
// $Id: G4DNAGenericIonsManager.cc,v 1.6 2009/06/10 13:32:36 mantero Exp $
// GEANT4 tag $Name: geant4-09-03 $

#include "G4DNAGenericIonsManager.hh"
#include "G4Alpha.hh"
#include "G4Ions.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAGenericIonsManager * G4DNAGenericIonsManager :: Instance(void)
{
 if (!theInstance)
  theInstance=new G4DNAGenericIonsManager;
 
 return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ParticleDefinition * G4DNAGenericIonsManager :: GetIon(const G4String & name)
{
 IonsMap::const_iterator i(map.find(name));
 
 if (i==map.end())
  return 0;
  
 return i->second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAGenericIonsManager :: G4DNAGenericIonsManager()
{
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //         shortlived          subType  anti_encoding
  //         excitation   
  
 G4Ions *helium;
 G4Ions *hydrogen;
 G4Ions *alphaPlus;
 G4Ions *positronium1s;
 G4Ions *positronium2s;
 
 
 helium=     new G4Ions(
			"helium",    3.727417*GeV,       0.0*MeV,  +0.0*eplus,
			0,              +1,             0,
			0,               0,             0,
			"nucleus",              +2,            +4,           0,
			true,                -1.0,             0,       
			false,		      "",               0,             
			0.0);

 alphaPlus=  new G4Ions("alpha+",    3.727417*GeV,       0.0*MeV,  +1.0*eplus,
                               1,              +1,             0,
                               0,               0,             0,
                       "nucleus",              +1,            +4,           0,
			true,            -1.0,             0, false,
			      "",               0,             0.0);

 hydrogen= new G4Ions("hydrogen",   0.9382723*GeV,       0.0*MeV,  +0.0*eplus,
                               0,              +1,             0,
                               0,               0,             0,
                       "nucleus",              +1,            +1,           0,
		        true,            -1.0,             0, false,
			      "",               0,             0.0);

 positronium1s= new G4Ions("Ps-1s",   2*electron_mass_c2,      0.0*MeV,  +0.0*eplus,
                               0,               0,             0,
                               0,               0,             0,
                       "nucleus",               0,             0,           0,
		            true,            -1.0,             0, false,
			      "",               0,             0.0);

 positronium2s= new G4Ions("Ps-2s",   2*electron_mass_c2,      0.0*MeV,  +0.0*eplus,
                               0,               0,             0,
                               0,               0,             0,
                       "nucleus",               0,             0,           0,
		            true,            -1.0,             0, false,
			      "",               0,             0.0);


 /*
 // molechules construction

 G4Ions* oxonium; // H3O -- it will become H3O+
 G4Ions* hydroxyl; // OH -- it will produce OH- too
 G4Ions* molHydrogen; // H2
 //G4Ions* hydroxide; // OH-
 G4Ions* hydroPeroxide; // H2O2
 G4Ions* water; // H2O -- it will become also H2O+


 G4double mass = 19.02*g/Avogadro - 11*electron_mass_c2;

 oxonium = new G4Ions("H3O",        mass,             0,  +11.0*eplus,
		      0,               0,             0,
		      0,               0,             0,
		      "molecule",      0,             0,           0,
		      true,         -1.0,             0, 
		      false,          "",             0,             
		      0.0); 
 
 mass = 17.00734*g/Avogadro - 9*electron_mass_c2;
 
 hydroxyl = new G4Ions("OH",        mass,             0,  +9.0*eplus,
		      0,               0,             0,
		      0,               0,             0,
		      "molecule",      0,             0,           0,
		      true,         -1.0,             0, 
		      false,          "",             0,             
		      0.0); 

 mass = 2.01588*g/Avogadro - 2*electron_mass_c2;

 molHydrogen = new G4Ions("H2",         mass,             0,  +2.0*eplus,
			  0,               0,             0,
			  0,               0,             0,
			  "molecule",      0,             0,           0,
			  true,         -1.0,             0, 
			  false,          "",             0,             
			  0.0); 

 mass = 34.01468*g/Avogadro - 18*electron_mass_c2;

 hydroPeroxide = new G4Ions("H2O2",      mass,             0, +18.0*eplus,
			   0,               0,             0,
			   0,               0,             0,
			   "molecule",      0,             0,           0,
			   true,         -1.0,             0, 
			   false,          "",             0,             
			   0.0); 
 
 mass = 18.015*g/Avogadro - 10*electron_mass_c2;

 water = new G4Ions("H2O",        mass,             0,  +10.0*eplus,
		    0,               0,             0,
		    0,               0,             0,
		    "molecule",      0,             0,           0,
		    true,         -1.0,             0, 
		    false,          "",             0,             
		    0.0); 

 map["H3O" ]  =oxonium;
 map["OH"  ]  =hydroxyl;
 map["H2"  ]  =molHydrogen;
 map["H2O2"]  =hydroPeroxide;
 map["H2O" ]  =water; 
 */


 map["helium"  ]=helium;
 map["hydrogen"]=hydrogen;
 map["alpha+"  ]=alphaPlus;
 map["alpha++" ]=G4Alpha::Alpha();
 map["Ps-1s"   ]=positronium1s;
 map["Ps-2s"   ]=positronium2s;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAGenericIonsManager * G4DNAGenericIonsManager::theInstance(0);
   
