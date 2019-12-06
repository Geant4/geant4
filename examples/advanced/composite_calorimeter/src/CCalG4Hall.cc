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
///////////////////////////////////////////////////////////////////////////////
// File: CCalG4Hall.cc
// Description: CCalG4Hall Geometry factory class to construct G4 geometry 
//              of the experimental hall
///////////////////////////////////////////////////////////////////////////////
#include "CCalG4Hall.hh"

#include "CCalMaterialFactory.hh"

#include "CCalG4Hcal.hh"
#include "CCalG4Ecal.hh"

#include "CCalutils.hh"

#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"

//#define debug
//#define pdebug

////////////////////////////////////////////////////////////////////
// CCalG4Hall constructor & destructor...
////////////////////////////////////////////////////////////////////

CCalG4Hall::CCalG4Hall(const G4String &name):
  CCalHall(name),CCalG4Able(name) {}

CCalG4Hall::~CCalG4Hall() {}

////////////////////////////////////////////////////////////////////
// CCalG4Hall methods...
////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* CCalG4Hall::constructIn(G4VPhysicalVolume* mother) {

  G4cout << "==>> Constructing CCalG4Hall..." << G4endl;

  ///////////////////////////////////////////////////////////////
  //Pointers to the Rotation Matrices and to the Materials
  CCalMaterialFactory* matfact       = CCalMaterialFactory::getInstance();

  //Experimental Hall. The mother volume....
#ifdef debug
  G4cout << tab << "Building experimental Hall geometry...." << G4endl;
#endif

  G4Box* solid = new G4Box(Name(), getDx_2Hall()*mm, getDy_2Hall()*mm,
                           getDy_2Hall()*mm);
#ifdef debug
  G4cout << solid->GetName() << " made of " << getMaterial() << " of dimension " 
       << getDx_2Hall()*mm << " " << getDy_2Hall()*mm << " " 
       << getDy_2Hall()*mm << " (all in mm)" << G4endl;
#endif

  G4Material* matter = matfact->findMaterial(getMaterial());
  G4LogicalVolume* glog = new G4LogicalVolume (solid, matter, Name(), 0, 0, 0);
  setVisType(CCalVisualisable::PseudoVolumes,glog);

  G4VPhysicalVolume* volume = new G4PVPlacement(0,G4ThreeVector(),Name(),
                                                glog,mother,false,0);
#ifdef pdebug
  G4String name("Null");
  if (mother != 0) name = mother->GetName();
  G4cout << glog->GetName() << " number 1 positioned in " << name
       << " at (0,0,0) with no rotation" << G4endl;
#endif  

  G4cout << "<<== End of CCalG4Hall construction ..." << G4endl;

  return volume;
}

void CCalG4Hall::constructDaughters(){
  //Hadron Calorimeter
  CCalG4Hcal* hcal = new CCalG4Hcal("HadronCalorimeter");
  addDetector(hcal);
  AddCCalG4Able(hcal);

  //Crystal matrix
  CCalG4Ecal* ecal = new CCalG4Ecal("CrystalMatrixModule");
  ecal->setType(CCalG4Ecal::module1);
  addDetector(ecal);
  AddCCalG4Able(ecal);
}
