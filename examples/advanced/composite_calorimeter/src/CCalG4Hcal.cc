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
// File: CCalG4Hcal.cc
// Description: CCalG4Hcal Factory class to construct the G4 geometry of the
//              hadron calorimeter
///////////////////////////////////////////////////////////////////////////////
#include <cmath>

#include "CCalG4Hcal.hh"

#include "CCalMaterialFactory.hh"
#include "CCalRotationMatrixFactory.hh"
#include "CCalSensitiveDetectors.hh"

#include "CCalutils.hh"

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

//#define debug
//#define ddebug
//#define pdebug
//#define sdebug

////////////////////////////////////////////////////////////////////
// CCalG4Hcal constructor & destructor...
////////////////////////////////////////////////////////////////////

CCalG4Hcal::CCalG4Hcal(const G4String &name):
  CCalHcal(name), CCalG4Able(name), sclLog(0), absLog(0) {}

CCalG4Hcal::~CCalG4Hcal(){
  if (sclLog)
    delete[] sclLog;
  if (absLog)
    delete[] absLog;
}

////////////////////////////////////////////////////////////////////
// CCalG4Hcal methods...
////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* CCalG4Hcal::constructIn(G4VPhysicalVolume* mother) {
  G4cout << "==>> Constructing CCalG4Hcal..." << G4endl;

  //Common logical volumes between methods.
#ifdef debug
  G4cout << tab << "Common logical volumes initialization: " 
       << getNScintillator() << " scintillaor and " << getNAbsorber()
       << " absorber layers." << G4endl;
#endif
  G4int i = 0;
  sclLog  = new ptrG4Log[getNScintillator()];
  absLog  = new ptrG4Log[getNAbsorber()];
  for (i=0; i < getNScintillator(); i++)
    sclLog[i]  = 0;
  for (i=0; i < getNAbsorber(); i++)
    absLog[i] = 0;

  //Pointers to the Materials
  CCalMaterialFactory* matfact       = CCalMaterialFactory::getInstance();

  //Mother volume
  G4Material* matter = matfact->findMaterial(getGenMat());
  G4VSolid*   solid  = new G4Box (Name(), getDx_2Cal()*mm, getDy_2Cal()*mm,
                                  getDy_2Cal()*mm);
  G4LogicalVolume* logh = new G4LogicalVolume(solid, matter, Name());
  setVisType(CCalVisualisable::PseudoVolumes,logh);
#ifdef debug
    G4cout << tab << Name() << " Box made of " << getGenMat()
         << " of dimension " << getDx_2Cal()*mm << " " << getDy_2Cal()*mm
         << " " << getDy_2Cal()*mm << G4endl;
#endif

  G4PVPlacement* hcal = new G4PVPlacement(0,G4ThreeVector(getXposCal()*mm,0,0),
                                          Name(), logh, mother, false, 1);
  G4String name("Null");
#ifdef pdebug
  if (mother != 0) name = mother->GetName();
  G4cout << Name() << " Number 1 positioned in " << name << " at ("
       << getXposCal()*mm << ",0,0) with no rotation" << G4endl;
#endif

  //Wall of the Boxes
  solid  = new G4Box (name, 0.5*getWallThickBox()*mm, getDy_2Box()*mm, 
                      getDy_2Box()*mm);
  matter = matfact->findMaterial(getBoxMat());
  name   = Name() + "Wall";
  G4LogicalVolume* logw = new G4LogicalVolume(solid, matter, name);
  setVisType(CCalVisualisable::Support,logw);
#ifdef debug
  G4cout << tab << name << " Box made of " << getBoxMat()
       << " of dimension " << 0.5*getWallThickBox()*mm << " " 
       << getDy_2Box()*mm << " " << getDy_2Box()*mm << G4endl;
#endif

  //Now the boxes
  ptrG4Log* logb = new ptrG4Log[getNBox()];
  matter = matfact->findMaterial(getGenMat());
  for (i=0; i<getNBox(); i++) {
    name   = Name() + "Box" + i;
    solid  = new G4Box (name, getDx_2Box()*mm, getDy_2Box()*mm, 
                        getDy_2Box()*mm);
    logb[i]= new G4LogicalVolume(solid, matter, name);
    setVisType(CCalVisualisable::PseudoVolumes,logb[i]);
#ifdef debug
    G4cout << tab << name << " Box made of " << getGenMat()
         << " of dimension " << getDx_2Box()*mm << " " << getDy_2Box()*mm
         << " " << getDy_2Box()*mm << G4endl;
#endif

    G4double xpos = -(getDx_2Box() - 0.5*getWallThickBox());
    new G4PVPlacement (0, G4ThreeVector(xpos*mm,0,0), logw, logw->GetName(), 
                       logb[i], false, 1);
#ifdef pdebug
    G4cout << logw->GetName() << " Number 1 positioned in " << name
         << " at (" << xpos*mm << ",0,0) with no rotation" << G4endl;
#endif
    xpos = (getDx_2Box() - 0.5*getWallThickBox());
    new G4PVPlacement (0, G4ThreeVector(xpos*mm,0,0), logw, logw->GetName(), 
                       logb[i], false, 2);
#ifdef pdebug
    G4cout << logw->GetName() << " Number 2 positioned in " << name
         << " at (" << xpos*mm << ",0,0) with no rotation" << G4endl;
#endif

    new G4PVPlacement (0, G4ThreeVector(getXposBox(i)*mm,0,0), logb[i], name, 
                       logh, false, i+1);
#ifdef pdebug
    G4cout << name << " Number " << i+1 << " positioned in " << logh->GetName()
         << " at (" << getXposBox(i)*mm << ",0,0) with no rotation" << G4endl;
#endif
  }

  //Loop over scintillator layers
  for (i=0; i<getNLayerScnt(); i++) {
    G4int lay = getTypeScnt(i);
    if (!sclLog[lay])
      sclLog[lay] = constructScintillatorLayer(lay);
    if (getMotherScnt(i) < 0 || getMotherScnt(i) >= getNScintillator()) {
      logw = logh;
    } else {
      logw = logb[getMotherScnt(i)];
    }
    G4double xpos = getXposScnt(i);
    new G4PVPlacement (0, G4ThreeVector(xpos*mm,0,0), sclLog[lay], 
                       sclLog[lay]->GetName(), logw, false, i+1);
#ifdef pdebug
    G4cout << sclLog[lay]->GetName() << " Number " << i+1 << " positioned in " 
         << logw->GetName() << " at (" << xpos*mm << ",0,0) with no rotation" 
         << G4endl;
#endif
  }

  //Loop over absorber layers
  for (i=0; i<getNLayerAbs(); i++) {
    G4int lay = getTypeAbs(i);
    if (!absLog[lay])
      absLog[lay] = constructAbsorberLayer(lay);
    if (getMotherAbs(i) < 0 || getMotherAbs(i) >= getNAbsorber()) {
      logw = logh;
    } else {
      logw = logb[getMotherAbs(i)];
    }
    G4double xpos = getXposAbs(i);
    new G4PVPlacement (0, G4ThreeVector(xpos*mm,0,0), absLog[lay], 
                       absLog[lay]->GetName(), logw, false, i+1);
#ifdef pdebug
    G4cout << absLog[lay]->GetName() << " Number " << i+1 << " positioned in " 
         << logw->GetName() << " at (" << xpos*mm << ",0,0) with no rotation" 
         << G4endl;
#endif
  }

  delete [] logb;

  G4cout << "<<== End of CCalG4Hcal construction ..." << G4endl;

  return hcal;
}


G4LogicalVolume* CCalG4Hcal::constructScintillatorLayer(G4int lay) {

  //Pointers to the Materials
  CCalMaterialFactory* matfact       = CCalMaterialFactory::getInstance();

  //The scintillator layer
  G4Material* matter = matfact->findMaterial(getGenMat());
  G4String    name   = Name() + "ScntLayer" + lay;
  G4VSolid*   solid  = new G4Box (name, getDx_2ScntLay(lay)*mm, 
                                  getDy_2ScntLay(lay)*mm,
                                  getDy_2ScntLay(lay)*mm);
  G4LogicalVolume* log = new G4LogicalVolume(solid, matter, name);
  setVisType(CCalVisualisable::PseudoVolumes,log);
#ifdef debug
  G4cout << tab << name << " Box made of " << getGenMat() << " of dimension " 
       << getDx_2ScntLay(lay)*mm << " " << getDy_2ScntLay(lay)*mm << " " 
       << getDy_2ScntLay(lay)*mm << G4endl;
#endif

  G4LogicalVolume* logd;
  G4double         xpos;
  //Wrappers if any
  if (getDx_2Wrap(lay) > 0) {
    name   = Name() + "ScntWrapper" + lay;
    matter = matfact->findMaterial(getWrapMat());
    solid  = new G4Box (name, getDx_2Wrap(lay)*mm, 
                        getDy_2ScntLay(lay)*mm, getDy_2ScntLay(lay)*mm);
    logd   = new G4LogicalVolume(solid, matter, name);
    setVisType(CCalVisualisable::Support,logd);
#ifdef debug
    G4cout << tab << name << " Box made of " << getWrapMat() << " of dimension " 
         << getDx_2Wrap(lay)*mm << " " << getDy_2ScntLay(lay)*mm << " " 
         << getDy_2ScntLay(lay)*mm << G4endl;
#endif
    xpos   =-(getDx_2ScntLay(lay)-getDx_2Wrap(lay));
    new G4PVPlacement(0, G4ThreeVector(xpos*mm,0,0), logd, name, log, false,1);
#ifdef pdebug
    G4cout << logd->GetName() << " Number 1 positioned in " << log->GetName() 
         << " at (" << xpos*mm << ",0,0) with no rotation" << G4endl;
#endif
    xpos   = (getDx_2ScntLay(lay)-getDx_2Wrap(lay));
    new G4PVPlacement(0, G4ThreeVector(xpos*mm,0,0), logd, name, log, false,2);
#ifdef pdebug
    G4cout << logd->GetName() << " Number 2 positioned in " << log->GetName() 
         << " at (" << xpos*mm << ",0,0) with no rotation" << G4endl;
#endif
  }

  //Plastic covers
  matter = matfact->findMaterial(getPlasMat());
  name   = Name() + "FrontPlastic" + lay;
  solid  = new G4Box (name, getDx_2FrontP(lay)*mm, getDy_2ScntLay(lay)*mm, 
                      getDy_2ScntLay(lay)*mm);
  logd   = new G4LogicalVolume(solid, matter, name);
  setVisType(CCalVisualisable::Cable,logd);
#ifdef debug
  G4cout << tab << name << " Box made of " << getPlasMat() << " of dimension " 
       << getDx_2FrontP(lay)*mm << " " << getDy_2ScntLay(lay)*mm << " " 
       << getDy_2ScntLay(lay)*mm << G4endl;
#endif
  xpos   =-getDx_2ScntLay(lay)+2.*getDx_2Wrap(lay)+getDx_2FrontP(lay);
  new G4PVPlacement(0, G4ThreeVector(xpos*mm,0,0), logd, name, log, false,1);
#ifdef pdebug
  G4cout << logd->GetName() << " Number 1 positioned in " << log->GetName() 
       << " at (" << xpos*mm << ",0,0) with no rotation" << G4endl;
#endif
  name   = Name() + "BackPlastic" + lay;
  solid  = new G4Box (name, getDx_2BackP(lay)*mm, getDy_2ScntLay(lay)*mm, 
                      getDy_2ScntLay(lay)*mm);
  logd   = new G4LogicalVolume(solid, matter, name);
  setVisType(CCalVisualisable::Cable,logd);
#ifdef debug
  G4cout << tab << name << " Box made of " << getPlasMat() << " of dimension " 
       << getDx_2BackP(lay)*mm << " " << getDy_2ScntLay(lay)*mm << " " 
       << getDy_2ScntLay(lay)*mm << G4endl;
#endif
  xpos   =(-getDx_2ScntLay(lay)+2.*getDx_2Wrap(lay)+2.*getDx_2FrontP(lay)+
           2.*getDx_2Scnt(lay)+getDx_2BackP(lay));
  new G4PVPlacement(0, G4ThreeVector(xpos*mm,0,0), logd, name, log, false,1);
#ifdef pdebug
  G4cout << logd->GetName() << " Number 1 positioned in " << log->GetName() 
       << " at (" << xpos*mm << ",0,0) with no rotation" << G4endl;
#endif

  //Now the scintillators
  matter = matfact->findMaterial(getScntMat());
  name   = Name() + "Scintillator" + lay;
  solid  = new G4Box (name, getDx_2Scnt(lay)*mm, getDy_2ScntLay(lay)*mm, 
                      getDy_2ScntLay(lay)*mm);
  logd   = new G4LogicalVolume(solid, matter, name);
  setVisType(CCalVisualisable::Sensitive,logd);
  allSensitiveLogs.push_back(logd);
#ifdef debug
  G4cout << tab << name << " Box made of " << getScntMat() << " of dimension " 
       << getDx_2Scnt(lay)*mm << " " << getDy_2ScntLay(lay)*mm << " " 
       << getDy_2ScntLay(lay)*mm << G4endl;
#endif
  xpos   =(-getDx_2ScntLay(lay)+2.*getDx_2Wrap(lay)+2.*getDx_2FrontP(lay)+
           getDx_2Scnt(lay));
  new G4PVPlacement(0, G4ThreeVector(xpos*mm,0,0), logd, name, log, false,1);
#ifdef pdebug
  G4cout << logd->GetName() << " Number 1 positioned in " << log->GetName() 
       << " at (" << xpos*mm << ",0,0) with no rotation" << G4endl;
#endif

  return log;
}


G4LogicalVolume* CCalG4Hcal::constructAbsorberLayer(G4int lay) {

  //Pointers to the Materials
  CCalMaterialFactory* matfact       = CCalMaterialFactory::getInstance();

  //Now the absorber layer
  G4Material* matter = matfact->findMaterial(getAbsMat());
  G4String    name   = Name() + "Absorber" + lay;
  G4VSolid*   solid  = new G4Box (name, getDx_2Abs(lay)*mm, getDy_2Abs()*mm,
                                  getDy_2Abs()*mm);
  G4LogicalVolume* log = new G4LogicalVolume(solid, matter, name);
  setVisType(CCalVisualisable::Absorber,log);
#ifdef debug
  G4cout << tab << name << " Box made of " << getAbsMat() << " of dimension " 
       << getDx_2Abs(lay)*mm << " " << getDy_2Abs()*mm << " " 
       << getDy_2Abs()*mm << G4endl;
#endif

  return log;
}


void CCalG4Hcal::constructSensitive(){

  if (allSensitiveLogs.size()>0) {
    CCalSensitiveDetectors* sensDets = CCalSensitiveDetectors::getInstance();
    G4String SDname = Name();
    for (std::vector<ptrG4Log>::iterator iter=allSensitiveLogs.begin(); 
         iter<allSensitiveLogs.end(); iter++) {
      sensDets->registerVolume(SDname, (*iter));
#ifdef sdebug
      G4cout << "Register volume " << (*iter)->GetName() << " for" << SDname 
           << G4endl;
#endif
    }
  } else {
    G4cerr << "CCalG4Hcal ERROR: Could not construct Sensitive Detector" 
           << G4endl;
  }
}

void CCalG4Hcal::constructDaughters() {}
