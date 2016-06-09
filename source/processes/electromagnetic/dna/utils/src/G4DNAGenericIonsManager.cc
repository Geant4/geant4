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

#include "G4DNAGenericIonsManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
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
 
 G4Ions *carbon;
 G4Ions *nitrogen;
 G4Ions *oxygen;
 G4Ions *iron;

 iron=     new G4Ions(
			"iron",    52.5672*GeV,       0.0*MeV,  +26.0*eplus,
			0,              +1,             0,
			0,               0,             0,
			"nucleus",              +26,            +56,           0,
			true,                -1.0,             0,       
			false,		      "",               0,             
			0.0);

 oxygen=   new G4Ions(
			"oxygen",    15.0074*GeV,       0.0*MeV,  +8.0*eplus,
			0,              +1,             0,
			0,               0,             0,
			"nucleus",              +8,            +16,           0,
			true,                -1.0,             0,       
			false,		      "",               0,             
			0.0);


 nitrogen= new G4Ions(
			"nitrogen",    13.132*GeV,       0.0*MeV,  +7.0*eplus,
			0,              +1,             0,
			0,               0,             0,
			"nucleus",              +7,            +14,           0,
			true,                -1.0,             0,       
			false,		      "",               0,             
			0.0);

 carbon=   new G4Ions(
			"carbon",    11.267025440*GeV,       0.0*MeV,  +6.0*eplus,
			0,              +1,             0,
			0,               0,             0,
			"nucleus",              +6,            +12,           0,
			true,                -1.0,             0,       
			false,		      "",               0,             
			0.0);
 
 helium=   new G4Ions(
			"helium",    3.727417*GeV,       0.0*MeV,  +0.0*eplus,
			0,              +1,             0,
			0,               0,             0,
			"nucleus",              +2,            +4,           0,
			true,                -1.0,             0,       
			false,		      "",               0,             
			0.0);

 alphaPlus= new G4Ions("alpha+",    3.727417*GeV,       0.0*MeV,  +1.0*eplus,
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


 map["helium"  ]=helium;
 map["hydrogen"]=hydrogen;
 map["alpha+"  ]=alphaPlus;
 map["alpha++" ]=G4Alpha::Alpha();
 map["Ps-1s"   ]=positronium1s;
 map["Ps-2s"   ]=positronium2s;
 map["carbon"  ]=carbon;
 map["nitrogen"]=nitrogen;
 map["oxygen"  ]=oxygen;
 map["iron"    ]=iron;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAGenericIonsManager * G4DNAGenericIonsManager::theInstance(0);
   
