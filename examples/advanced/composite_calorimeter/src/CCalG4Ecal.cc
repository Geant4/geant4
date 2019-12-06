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
// File: CCalG4Ecal.cc
// Description: CCalG4Ecal Factory class to construct the G4 geometry of the
//              electromagnetic calorimeter
///////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "CCalG4Ecal.hh"

#include "CCalMaterialFactory.hh"
#include "CCalRotationMatrixFactory.hh"
#include "CCalSensitiveDetectors.hh"

#include "CCalutils.hh"

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4Trd.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

//#define debug
//#define ddebug
//#define pdebug
//#define sdebug

//Initialize static logical volumes
G4LogicalVolume* CCalG4Ecal::crystalmatrixLog = 0;

//Initialize static prefix name
G4String CCalG4Ecal::idName = "CrystalMatrix";

////////////////////////////////////////////////////////////////////
// CCalG4Ecal constructor & destructor...
////////////////////////////////////////////////////////////////////

CCalG4Ecal::CCalG4Ecal(const G4String &name):
  CCalEcal(name), CCalG4Able(name), type(module1) {}

CCalG4Ecal::~CCalG4Ecal() {}

////////////////////////////////////////////////////////////////////
// CCalG4Ecal methods...
////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* CCalG4Ecal::constructIn(G4VPhysicalVolume* mother) {
  G4cout << "==>> Constructing CCalG4Ecal..." << G4endl;

  ///////////////////////////////////////////////////////////////
  // Construction of global volume as a Box

  if (!crystalmatrixLog) {
    crystalmatrixLog = constructGlobal();
  }
  CCalRotationMatrixFactory* rotfact = CCalRotationMatrixFactory::getInstance();

  G4double x, y, z;
  if (mother != 0) {
    x = getXpos()*mm;
    y = getYpos()*mm;
    z = getZpos()*mm;
  } else {
    x = y = z = 0;
  }
    
  int num;
  if (type == module2) {
    num = 2;
  } else {
    num = 1;
  }
#ifdef pdebug
  G4String name("Null");
  if (mother != 0) name = mother->GetName();
  G4cout << crystalmatrixLog->GetName() << " Number " << num << " positioned in "
       << name << " at (" << x << ", " << y << ", " << z << ")";
#endif

  G4RotationMatrix* cmrot = 0;
  if (mother != 0) {
    G4String rotstr = idName + num;
    cmrot  = rotfact->findMatrix(rotstr);
    if (!cmrot) {
#ifdef ddebug
      G4cout << "Creating a new rotation: " << rotstr << tab 
           << getThetaX()*deg << "," << getPhiX()*deg << "," 
           << getThetaY()*deg << "," << getPhiY()*deg << "," 
           << getThetaZ()*deg << "," << getPhiZ()*deg << G4endl;
#endif
      cmrot = rotfact->AddMatrix(rotstr, getThetaX()*deg, getPhiX()*deg, 
                                 getThetaY()*deg, getPhiY()*deg,
                                 getThetaZ()*deg, getPhiZ()*deg);
    } // if !cmrot
#ifdef pdebug
    G4cout << " rotation by (" <<  getThetaX() << ", " << getPhiX() << ", " 
         << getThetaY() << "," << getPhiY() << ", "  << getThetaZ() << ", " 
         << getPhiZ() << ")" << G4endl;
#endif
  } else {
#ifdef pdebug
    G4cout << " without rotation..." << G4endl;
#endif
  }

  G4PVPlacement* crystalmatrix;
  if (mother != 0) {
    crystalmatrix = new G4PVPlacement(cmrot, G4ThreeVector(x,y,z),
                                      crystalmatrixLog, idName,
                                      mother->GetLogicalVolume(), false, num);
  } else {
    crystalmatrix = new G4PVPlacement(cmrot, G4ThreeVector(x,y,z),
                                      idName, crystalmatrixLog,
                                      mother, false, num);
  }
  G4cout << "<<== End of CCalG4Ecal construction ..." << G4endl;

  return crystalmatrix;
}


G4LogicalVolume* CCalG4Ecal::constructGlobal() {

  //Pointers to the Materials and Rotation Matrix factory
  CCalMaterialFactory* matfact       = CCalMaterialFactory::getInstance();
  CCalRotationMatrixFactory* rotfact = CCalRotationMatrixFactory::getInstance();
  
  G4Material* matter = matfact->findMaterial(getGenMat());
  G4VSolid* solid = new G4Box (idName, 0.5*getWidBox()*mm, 0.5*getWidBox()*mm,
                               0.5*getLengBox()*mm);
#ifdef debug
  G4cout << tab << idName << " Box made of " << getGenMat() << " of dimension " 
       << 0.5*getWidBox()*mm << ", " << 0.5*getWidBox()*mm << ", "
       << 0.5*getLengBox()*mm << G4endl;
#endif
  G4LogicalVolume* glog = new G4LogicalVolume (solid, matter, idName);
  setVisType(CCalVisualisable::PseudoVolumes,glog);

  //Now the layers
  G4String name = idName + "Layer";
  matter = matfact->findMaterial(getLayMat());
  solid  = new G4Trd(name, getLayPar(0)*mm, getLayPar(1)*mm, getLayPar(2)*mm, 
                     getLayPar(3)*mm, getLayPar(4)*mm);
#ifdef debug
  G4cout << tab << name << " Trd made of " << getLayMat() << " of dimension " 
       << getLayPar(0)*mm << ", " << getLayPar(1)*mm << ", " << getLayPar(2)*mm
       << ", " << getLayPar(3)*mm << ", " << getLayPar(4)*mm << G4endl;
#endif
  G4LogicalVolume* laylog = new G4LogicalVolume (solid, matter, name);
  setVisType(CCalVisualisable::OtherServices,laylog);

  G4int i = 0;
  G4String rotstr;
  G4double xp, yp, zp, angle;
  G4double zshift = -0.5 * (getLengBox() - getCrystLength()) + getLengFront();
  G4RotationMatrix* rot = 0;
  for (i = 0; i < getLayNum(); i++) {
    angle  = 0.5 * getLayAngle() * (2*i + 1 - getLayNum());
    xp     = angle * (getLayPar(4) + getLayRadius()) * mm;
    zp     = (zshift + getLayPar(0)*std::abs(std::sin(angle))) * mm;
    rotstr = idName + "Layer" + i;
    rot    = rotfact->findMatrix(rotstr);
    if (!rot) {
#ifdef ddebug
      G4cout << "Creating a new rotation: " << rotstr << tab 
           << (90.0*deg+angle) << "," << 0.0*deg << "," << 90.0*deg << "," 
           << 90.0*deg << "," << angle << "," << 0.0*deg << G4endl;
#endif
      rot = rotfact->AddMatrix(rotstr, (90.0*deg+angle), 0.0*deg, 90.0*deg,
                               90.0*deg, angle, 0.0*deg);
    }
    new G4PVPlacement(rot, G4ThreeVector(xp,0.,zp), laylog, name, glog,
                      false, i+1);
#ifdef pdebug
    G4cout << laylog->GetName() << " number " << i+1 << " positioned in " 
         << glog->GetName()  << " at (" << xp << ", 0," << zp
         << ") with rotation angle " << angle/deg << G4endl;
#endif
  }

  //Now the crystals
  name   = idName + "Crystal";
  matter = matfact->findMaterial(getCrystMat());
  solid  = new G4Trd(name, getCrystPar(0)*mm, getCrystPar(1)*mm, 
                     getCrystPar(2)*mm, getCrystPar(3)*mm, getCrystPar(4)*mm);
#ifdef debug
  G4cout << tab << name << " Trd made of " << getCrystMat() << " of dimension " 
       << getCrystPar(0)*mm << ", " << getCrystPar(1)*mm << ", " 
       << getCrystPar(2)*mm << ", " << getCrystPar(3)*mm << ", " 
       << getCrystPar(4)*mm << G4endl;
#endif

  G4LogicalVolume* detLog = new G4LogicalVolume (solid, matter, name);
  setVisType(CCalVisualisable::Sensitive,detLog);
  sensitiveLogs.push_back(detLog);
  for (i = 0; i < getCrystNum(); i++) {
    angle  = 0.5 * getLayAngle() * (2*i + 1 - getCrystNum());
    yp     = angle * (getCrystPar(4) + getLayRadius()) * mm;
    zp     = (getCrystPar(0)*std::abs(std::sin(angle)) - getCrystTol()) * mm;
    rotstr = idName + "Crystal" + i;
    rot    = rotfact->findMatrix(rotstr);
    if (!rot) {
#ifdef ddebug
      G4cout << "Creating a new rotation: " << rotstr << tab << 90.0*deg << ","
           << 0.0*deg << "," << (90.0*deg+angle) << "," << 0.0*deg << "," 
           << angle << "," << 90.0*deg << G4endl;
#endif
      rot = rotfact->AddMatrix(rotstr, 90.0*deg, 0.0*deg, (90.0*deg+angle),
                               90.0*deg, angle, 90.0*deg);
    }
    new G4PVPlacement(rot, G4ThreeVector(0,yp,zp), detLog, name, laylog,
                      false, i+1);
#ifdef pdebug
    G4cout << detLog->GetName() << " number " << i+1 << " positioned in " 
         << laylog->GetName()  << " at (0," << yp << "," << zp
         << ") with rotation angle " << angle/deg << G4endl;
#endif
  }

  //Support boxes
  name   = idName + "Support";
  matter = matfact->findMaterial(getSuppMat());
  solid  = new G4Box (name, 0.5*getDxSupp()*mm, 0.5*getDySupp()*mm, 
                      0.5*getDzSupp()*mm);
#ifdef debug
  G4cout << tab << name << " Box made of " << getSuppMat() << " of dimension " 
       << 0.5*getDxSupp()*mm << ", " << 0.5*getDySupp()*mm << ", "
       << 0.5*getDzSupp()*mm << G4endl;
#endif
  G4LogicalVolume* slog = new G4LogicalVolume (solid, matter, name);
  setVisType(CCalVisualisable::Support,slog);

  zp   = (-0.5 * getLengBox() + getCrystLength() + getLengFront() +
          0.5 * getDzSupp() + getDistSupp()) * mm;
  for (i = 0; i < getCrystNum(); i++) {
    yp   = getLayPar(1) * (2*i + 1 - getCrystNum()) * mm;
    new G4PVPlacement(0, G4ThreeVector(0,yp,zp), slog, name, glog,
                      false, i+1);
#ifdef pdebug
    G4cout << slog->GetName() << " number " << i+1 << " positioned in " 
         << glog->GetName()  << " at (0," << yp << "," << zp
         << ") with no rotation" << G4endl;
#endif
  }

  return glog;
}

void CCalG4Ecal::constructSensitive() {

#ifdef debug
  G4cout << "Now registering CrystalMatrix LogicalVolume's to SD's:" << G4endl;
#endif
  if (sensitiveLogs.size()>0) {
    CCalSensitiveDetectors* sensDets = CCalSensitiveDetectors::getInstance();
    G4String SDname = idName;
    for(std::vector<ptrG4Log>::iterator iter=sensitiveLogs.begin(); 
                                   iter<sensitiveLogs.end(); iter++) {
      sensDets->registerVolume(SDname, (*iter));
#ifdef sdebug
      G4cout << "Register volume " << (*iter)->GetName() << " for" << SDname 
           << G4endl;
#endif
    }
  }

}
  
