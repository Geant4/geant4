///////////////////////////////////////////////////////////////////////////////
// File: G4CrystalMatrix.cc
// Date: 08/99 S.B.
// Modifications: 10/99 Naming convention
//                03/00 In OSCAR
//                06/09/01 S.B. Top level in logical volume
///////////////////////////////////////////////////////////////////////////////

#include "G4CrystalMatrix.hh"

#include "CCalMaterialFactory.hh"
#include "CCalRotationMatrixFactory.hh"
#include "CCalSensitiveDetectors.hh"

#include "utils.hh"
#include <math.h>

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
G4LogicalVolume* G4CrystalMatrix::crystalmatrixLog = 0;

//Initialize static prefix name
G4String G4CrystalMatrix::idName = "CrystalMatrix";

////////////////////////////////////////////////////////////////////
// G4CrystalMatrix constructor & destructor...
////////////////////////////////////////////////////////////////////

G4CrystalMatrix::G4CrystalMatrix(const G4String &name):
  CrystalMatrix(name), G4Able(name), type(module1) {}

G4CrystalMatrix::~G4CrystalMatrix() {}

////////////////////////////////////////////////////////////////////
// G4CrystalMatrix methods...
////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* G4CrystalMatrix::constructIn(G4VPhysicalVolume* mother) {
  cout << "==>> Constructing G4Crystalmatrix..." << endl;

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
  cout << crystalmatrixLog->GetName() << " Number " << num << " positioned in "
       << name << " at (" << x << ", " << y << ", " << z << ")";
#endif

  G4RotationMatrix* cmrot = 0;
  if (mother != 0) {
    G4String rotstr = idName + num;
    cmrot  = rotfact->findMatrix(rotstr);
    if (!cmrot) {
#ifdef ddebug
      cout << "Creating a new rotation: " << rotstr << tab 
	   << getThetaX()*deg << "," << getPhiX()*deg << "," 
	   << getThetaY()*deg << "," << getPhiY()*deg << "," 
	   << getThetaZ()*deg << "," << getPhiZ()*deg << endl;
#endif
      cmrot = rotfact->AddMatrix(rotstr, getThetaX()*deg, getPhiX()*deg, 
				 getThetaY()*deg, getPhiY()*deg,
				 getThetaZ()*deg, getPhiZ()*deg);
    } // if !cmrot
#ifdef pdebug
    cout << " rotation by (" <<  getThetaX() << ", " << getPhiX() << ", " 
	 << getThetaY() << "," << getPhiY() << ", "  << getThetaZ() << ", " 
	 << getPhiZ() << ")" << endl;
#endif
  } else {
#ifdef pdebug
    cout << " without rotation..." << endl;
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
  cout << "<<== End of G4CrystalMatrix construction ..." << endl;

  return crystalmatrix;
}


G4LogicalVolume* G4CrystalMatrix::constructGlobal() {

  //Pointers to the Materials and Rotation Matrix factory
  CCalMaterialFactory* matfact       = CCalMaterialFactory::getInstance();
  CCalRotationMatrixFactory* rotfact = CCalRotationMatrixFactory::getInstance();
  
  G4Material* matter = matfact->findMaterial(getGenMat());
  G4VSolid* solid = new G4Box (idName, 0.5*getWidBox()*mm, 0.5*getWidBox()*mm,
			       0.5*getLengBox()*mm);
#ifdef debug
  cout << tab << idName << " Box made of " << getGenMat() << " of dimension " 
       << 0.5*getWidBox()*mm << ", " << 0.5*getWidBox()*mm << ", "
       << 0.5*getLengBox()*mm << endl;
#endif
  G4LogicalVolume* glog = new G4LogicalVolume (solid, matter, idName);
  setVisType(Visualisable::PseudoVolumes,glog);

  //Now the layers
  G4String name = idName + "Layer";
  matter = matfact->findMaterial(getLayMat());
  solid  = new G4Trd(name, getLayPar(0)*mm, getLayPar(1)*mm, getLayPar(2)*mm, 
		     getLayPar(3)*mm, getLayPar(4)*mm);
#ifdef debug
  cout << tab << name << " Trd made of " << getLayMat() << " of dimension " 
       << getLayPar(0)*mm << ", " << getLayPar(1)*mm << ", " << getLayPar(2)*mm
       << ", " << getLayPar(3)*mm << ", " << getLayPar(4)*mm << endl;
#endif
  G4LogicalVolume* laylog = new G4LogicalVolume (solid, matter, name);
  setVisType(Visualisable::OtherServices,laylog);

  G4int i = 0;
  G4String rotstr;
  G4double xp, yp, zp, angle;
  G4double zshift = -0.5 * (getLengBox() - getCrystLength()) + getLengFront();
  G4RotationMatrix* rot = 0;
  for (i = 0; i < getLayNum(); i++) {
    angle  = 0.5 * getLayAngle() * (2*i + 1 - getLayNum());
    xp     = angle * (getLayPar(4) + getLayRadius()) * mm;
    zp     = (zshift + getLayPar(0)*abs(sin(angle))) * mm;
    rotstr = idName + "Layer" + i;
    rot    = rotfact->findMatrix(rotstr);
    if (!rot) {
#ifdef ddebug
      cout << "Creating a new rotation: " << rotstr << tab 
	   << (90.0*deg+angle) << "," << 0.0*deg << "," << 90.0*deg << "," 
	   << 90.0*deg << "," << angle << "," << 0.0*deg << endl;
#endif
      rot = rotfact->AddMatrix(rotstr, (90.0*deg+angle), 0.0*deg, 90.0*deg,
			       90.0*deg, angle, 0.0*deg);
    }
    new G4PVPlacement(rot, G4ThreeVector(xp,0.,zp), laylog, name, glog,
		      false, i+1);
#ifdef pdebug
    cout << laylog->GetName() << " number " << i+1 << " positioned in " 
         << glog->GetName()  << " at (" << xp << ", 0," << zp
	 << ") with rotation angle " << angle/deg << endl;
#endif
  }

  //Now the crystals
  name   = idName + "Crystal";
  matter = matfact->findMaterial(getCrystMat());
  solid  = new G4Trd(name, getCrystPar(0)*mm, getCrystPar(1)*mm, 
		     getCrystPar(2)*mm, getCrystPar(3)*mm, getCrystPar(4)*mm);
#ifdef debug
  cout << tab << name << " Trd made of " << getCrystMat() << " of dimension " 
       << getCrystPar(0)*mm << ", " << getCrystPar(1)*mm << ", " 
       << getCrystPar(2)*mm << ", " << getCrystPar(3)*mm << ", " 
       << getCrystPar(4)*mm << endl;
#endif

  G4LogicalVolume* detLog = new G4LogicalVolume (solid, matter, name);
  setVisType(Visualisable::Sensitive,detLog);
  sensitiveLogs.push_back(detLog);
  for (i = 0; i < getCrystNum(); i++) {
    angle  = 0.5 * getLayAngle() * (2*i + 1 - getCrystNum());
    yp     = angle * (getCrystPar(4) + getLayRadius()) * mm;
    zp     = (getCrystPar(0)*abs(sin(angle)) - getCrystTol()) * mm;
    rotstr = idName + "Crystal" + i;
    rot    = rotfact->findMatrix(rotstr);
    if (!rot) {
#ifdef ddebug
      cout << "Creating a new rotation: " << rotstr << tab << 90.0*deg << ","
	   << 0.0*deg << "," << (90.0*deg+angle) << "," << 0.0*deg << "," 
	   << angle << "," << 90.0*deg << endl;
#endif
      rot = rotfact->AddMatrix(rotstr, 90.0*deg, 0.0*deg, (90.0*deg+angle),
			       90.0*deg, angle, 90.0*deg);
    }
    new G4PVPlacement(rot, G4ThreeVector(0,yp,zp), detLog, name, laylog,
		      false, i+1);
#ifdef pdebug
    cout << detLog->GetName() << " number " << i+1 << " positioned in " 
         << laylog->GetName()  << " at (0," << yp << "," << zp
	 << ") with rotation angle " << angle/deg << endl;
#endif
  }

  //Support boxes
  name   = idName + "Support";
  matter = matfact->findMaterial(getSuppMat());
  solid  = new G4Box (name, 0.5*getDxSupp()*mm, 0.5*getDySupp()*mm, 
		      0.5*getDzSupp()*mm);
#ifdef debug
  cout << tab << name << " Box made of " << getSuppMat() << " of dimension " 
       << 0.5*getDxSupp()*mm << ", " << 0.5*getDySupp()*mm << ", "
       << 0.5*getDzSupp()*mm << endl;
#endif
  G4LogicalVolume* slog = new G4LogicalVolume (solid, matter, name);
  setVisType(Visualisable::Support,slog);

  zp   = (-0.5 * getLengBox() + getCrystLength() + getLengFront() +
	  0.5 * getDzSupp() + getDistSupp()) * mm;
  for (i = 0; i < getCrystNum(); i++) {
    yp   = getLayPar(1) * (2*i + 1 - getCrystNum()) * mm;
    new G4PVPlacement(0, G4ThreeVector(0,yp,zp), slog, name, glog,
		      false, i+1);
#ifdef pdebug
    cout << slog->GetName() << " number " << i+1 << " positioned in " 
         << glog->GetName()  << " at (0," << yp << "," << zp
	 << ") with no rotation" << endl;
#endif
  }

  return glog;
}

void G4CrystalMatrix::constructSensitive() {

#ifdef debug
  cout << "Now registering CrystalMatrix LogicalVolume's to SD's:" << endl;
#endif
  if (sensitiveLogs.size()>0) {
    CCalSensitiveDetectors* sensDets = CCalSensitiveDetectors::getInstance();
    G4String SDname = idName;
    for(vector<ptrG4Log>::iterator iter=sensitiveLogs.begin(); 
	                           iter<sensitiveLogs.end(); iter++) {
      sensDets->registerVolume(SDname, (*iter));
#ifdef sdebug
      cout << "Register volume " << (*iter)->GetName() << " for" << SDname 
	   << endl;
#endif
    }
  }

}

			  
