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
// File: CCalRotationMatrixFactory.cc
// Description: CCalRotationFactory is a factory class to define all rotation
//              matrices used in geometry building
///////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <stdlib.h>

#include "CCalRotationMatrixFactory.hh"

#include "CCalutils.hh"

#include "G4SystemOfUnits.hh"

//#define debug
//#define ddebug

CCalRotationMatrixFactory * CCalRotationMatrixFactory::instance = 0;
G4String CCalRotationMatrixFactory::file="";

CCalRotationMatrixFactory* CCalRotationMatrixFactory::getInstance(const G4String & rotfile){
  if (rotfile=="" || rotfile==file)
    return getInstance();
  else if (file=="") {
    file=rotfile;
    return getInstance();
  } else {
    G4cerr << "ERROR: Trying to get Rotation Matrices from " << rotfile 
         << " when previously were retrieved from " << file <<"." << G4endl;
    return 0;
  }
}


CCalRotationMatrixFactory* CCalRotationMatrixFactory::getInstance(){
  if (file=="") {
    G4cerr << "ERROR: You haven't defined which file to use for materials in "
         << "CCalRotationMatrixFactory::getInstance(G4String)" << G4endl;
    return 0;
  }

  if (instance==0) {
    instance = new CCalRotationMatrixFactory;
    return instance;
  }
  else
    return instance;
}

void CCalRotationMatrixFactory::setFileName(const G4String& rotfile) {
  if (rotfile!=file && file!="") {
    G4cerr << "ERROR: Trying to change Rotation Matrices file name to " 
         << rotfile << "." << G4endl;
    G4cerr << "       Using previous file: " << file << G4endl;
  }
  file=rotfile;
}

CCalRotationMatrixFactory::~CCalRotationMatrixFactory(){
  G4RotationMatrixTableIterator i;
  for(i=theMatrices.begin(); i != theMatrices.end(); ++i) {
    delete (*i).second;
  };
  theMatrices.clear();
}

G4RotationMatrix* CCalRotationMatrixFactory::findMatrix(const G4String & rot) {
  G4RotationMatrix* retrot=0;
  //Rotation :NULL is no rotation so a null pointer is returned
  if (rot != ":NULL") {
    //retrot untouched if rot not found!!!
    G4RotationMatrixTableIterator it = theMatrices.find(rot);
    if (it != theMatrices.end())
      retrot = (*it).second;
  }
  
  return retrot; //!!!Maybe a treatment on not-found case needed.
}

G4RotationMatrix* CCalRotationMatrixFactory::AddMatrix(const G4String& name, 
                                                       G4double th1, 
                                                       G4double phi1, 
                                                       G4double th2, 
                                                       G4double phi2, 
                                                       G4double th3, 
                                                       G4double phi3){
  G4double sinth1, sinth2,  sinth3, costh1, costh2, costh3;
  G4double sinph1, sinph2, sinph3, cosph1, cosph2, cosph3;
  G4double TH1 = th1/deg, TH2 = th2/deg, TH3 = th3/deg;
  G4double PH1 = phi1/deg, PH2 = phi2/deg, PH3 = phi3/deg;
                
  if (TH1 == 0.0 || TH1 == 360) {
    sinth1 = 0.0; costh1 = 1.0;
  } else if (TH1 == 90.0 || TH1 == -270) {
    sinth1 = 1.0; costh1 = 0.0;
  } else if (TH1 == 180.0 || TH1 == -180.0) {
    sinth1 = 0.0; costh1 = -1.0;
  } else if (TH1 == 270.0 || TH1 == -90.0) {
    sinth1 = -1.0; costh1 = 0.0;
  } else {
    sinth1 = std::sin(th1); costh1 = std::cos(th1);
  }
  
  if (TH2 == 0.0 || TH2 == 360) {
    sinth2 = 0.0; costh2 = 1.0;
  } else if (TH2 == 90.0 || TH2 == -270) {
    sinth2 = 1.0; costh2 = 0.0;
  } else if (TH2 == 180.0 || TH2 == -180.0) {
    sinth2 = 0.0; costh2 = -1.0;
  } else if (TH2 == 270.0 || TH2 == -90.0) {
    sinth2 = -1.0; costh2 = 0.0;
  } else {
    sinth2 = std::sin(th2); costh2 = std::cos(th2);
  }
                
  if (TH3 == 0.0 || TH3 == 360) {
    sinth3 = 0.0; costh3 = 1.0;
  } else if (TH3 == 90.0 || TH2 == -270) {
    sinth3 = 1.0; costh3 = 0.0;
  } else if (TH3 == 180.0 || TH3 == -180.0) {
    sinth3 = 0.0; costh3 = -1.0;
  } else if (TH3 == 270.0 || TH3 == -90.0) {
    sinth3 = -1.0; costh3 = 0.0;
  } else {
    sinth3 = std::sin(th3); costh3 = std::cos(th3);
  }
      
  if (PH1 == 0.0 || PH1 == 360) {
    sinph1 = 0.0; cosph1 = 1.0;
  } else if (PH1 == 90.0 || PH1 == -270) {
    sinph1 = 1.0; cosph1 = 0.0;
  } else if (PH1 == 180.0 || PH1 == -180.0) {
    sinph1 = 0.0; cosph1 = -1.0;
  } else if (PH1 == 270.0 || PH1 == -90.0) {
    sinph1 = -1.0; cosph1 = 0.0;
  } else {
    sinph1 = std::sin(phi1); cosph1 = std::cos(phi1);
  }

  if (PH2 == 0.0 || PH2 == 360) {
    sinph2 = 0.0; cosph2 = 1.0;
  } else if (PH2 == 90.0 || PH2 == -270) {
    sinph2 = 1.0; cosph2 = 0.0;
  } else if (PH2 == 180.0 || PH2 == -180.0) {
    sinph2 = 0.0; cosph2 = -1.0;
  } else if (PH2 == 270.0 || PH2 == -90.0) {
    sinph2 = -1.0; cosph2 = 0.0;
  } else {
    sinph2 = std::sin(phi2); cosph2 = std::cos(phi2);
  }
                
  if (PH3 == 0.0 || PH3 == 360) {
    sinph3 = 0.0; cosph3 = 1.0;
  } else if (PH3 == 90.0 || PH3 == -270) {
    sinph3 = 1.0; cosph3 = 0.0;
  } else if (PH3 == 180.0 || PH3 == -180.0) {
    sinph3 = 0.0; cosph3 = -1.0;
  } else if (PH3 == 270.0 || PH3 == -90.0) {
    sinph3 = -1.0; cosph3 = 0.0;
  } else {
    sinph3 = std::sin(phi3); cosph3 = std::cos(phi3);
  }
                                    
  //xprime axis coordinates
  CLHEP::Hep3Vector xprime(sinth1*cosph1,sinth1*sinph1,costh1);
  //yprime axis coordinates
  CLHEP::Hep3Vector yprime(sinth2*cosph2,sinth2*sinph2,costh2);
  //zprime axis coordinates
  CLHEP::Hep3Vector zprime(sinth3*cosph3,sinth3*sinph3,costh3);

#ifdef ddebug
  G4cout << xprime << '\t';    G4cout << yprime << '\t';    G4cout << zprime << G4endl;
#endif
  G4RotationMatrix *rotMat = new G4RotationMatrix();
  rotMat->rotateAxes(xprime, yprime, zprime);
  if (*rotMat == G4RotationMatrix()) {
    // G4cerr << "WARNING: Matrix " << name << " will not be created as a rotation matrix." 
    // G4cerr << "WARNING: Matrix " << name << " is = identity matrix. It will not be created as a rotation matrix." << G4endl;
    delete rotMat;
    rotMat=0;
  } else {
    rotMat->invert();
    theMatrices[name]=rotMat;
#ifdef ddebug
    G4cout << *rotMat << G4endl;
#endif
  }

  return rotMat;
}

CCalRotationMatrixFactory::CCalRotationMatrixFactory():theMatrices(){

  G4String path = "NULL";
  if (std::getenv("CCAL_GLOBALPATH"))
    path = std::getenv("CCAL_GLOBALPATH");

  G4cout << " ==> Opening file " << file << "..." << G4endl;
  std::ifstream is;
  G4bool ok = openGeomFile(is, path, file);
  if (!ok) {
    G4ExceptionDescription ed;
    ed << "Could not open file " << file << " ... Exiting!" << G4endl;
    G4Exception("CCalRotationMatrixFactory::CCalRotationMatrixFactory()",
                "ccal002",
                FatalException,ed);
  }

  //////////////////////////////////////////////////
  // Find *DO ROTM
  findDO(is, G4String("ROTM"));

  char rubish[256];
  G4String name;

#ifdef debug
  G4cout << "     ==> Reading Rotation Matrices... " << G4endl;
  G4cout << "       Name\tTheta1\tPhi1\tTheta2\tPhi2\tTheta3\tPhi3"<< G4endl;
#endif
  
  is >> name;
  while(name!="*ENDDO") { 
    if (name.index("#.")==0) { //It is a comment.Skip line.
      is.getline(rubish,256,'\n');
    } else {
#ifdef debug
      G4cout << "       " << name <<'\t';
#endif
      G4double th1, phi1, th2, phi2, th3, phi3;
      //Get xprime axis angles
      is >> th1 >> phi1;
#ifdef debug
      G4cout << th1 << '\t' << phi1 << '\t';
#endif
      //Get yprime axis angles
      is >> th2 >> phi2;
#ifdef debug
      G4cout << th2 << '\t' << phi2 << '\t';
#endif
      //Get zprime axis angles
      is >> th3 >> phi3;
#ifdef debug
      G4cout << th3 << '\t' << phi3 << '\t';
#endif

      is.getline(rubish,256,'\n');
#ifdef debug
      G4cout << rubish << G4endl;
#endif

      AddMatrix(name, th1*deg, phi1*deg, th2*deg, phi2*deg, th3*deg, phi3*deg);
    }

    is >> name;
  };

  is.close();
  G4cout << "       "  << theMatrices.size() << " rotation matrices read in." << G4endl;
}


// 29-Jan-2004 A.R. : commented to avoid clashes with CLHEP.
//                    Streaming operators for rotation matrices are
//                    already defined in CLHEP::HepRotation.
// std::ostream& operator<<(std::ostream& os , const G4RotationMatrix & rot){
//   //  os << "( " << rot.xx() << tab << rot.xy() << tab << rot.xz() << " )" << G4endl;
//   //  os << "( " << rot.yx() << tab << rot.yy() << tab << rot.yz() << " )" << G4endl;
//   //  os << "( " << rot.zx() << tab << rot.zy() << tab << rot.zz() << " )" << G4endl;
// 
//   os << "[" 
//      << rot.thetaX()/deg << tab << rot.phiX()/deg << tab
//      << rot.thetaY()/deg << tab << rot.phiY()/deg << tab
//      << rot.thetaZ()/deg << tab << rot.phiZ()/deg << "]"
//      << G4endl;
// 
//   return os;
// }
