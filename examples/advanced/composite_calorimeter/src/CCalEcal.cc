///////////////////////////////////////////////////////////////////////////////
// File: CCalEcal.cc
// Description: CCalEcal Geometry factory class for crystal matrix
///////////////////////////////////////////////////////////////////////////////
#include "CCalEcal.hh"

#include "g4std/fstream"
#include "CCalutils.hh"

//#define debug

CCalEcal::~CCalEcal() {}

int CCalEcal::readFile() {
  ///////////////////////////////////////////////////////////////
  //Let's open the file
  G4cout << " ==> Opening file " << File() << " to read elements..."
       << G4endl;

  G4std::ifstream is;
  bool ok = openGeomFile(is, pathName, File());
  if (!ok)
    return 0;

  // Find *DO CrystalMatrix 
  findDO(is, G4String("CrystalMatrix"));

  //Let's read overall box dimensions and positions
  readName(is,genMat);
  is >> widBox >> lengBox >> xpos >> ypos >> zpos >> thetaX >> phiX
     >> thetaY >> phiY >> thetaZ >> phiZ >> jump;
#ifdef debug
  G4cout << tab << "General material: " << genMat << "\tBox dimensions "
       << widBox << ", " << lengBox << G4endl;
  G4cout << tab << "Positioned at (" << xpos << ", " << ypos << ", " << zpos
       << ") with rotation (" << thetaX << ", " << phiX << ", " << thetaY
       << ", " << phiY << ", " << thetaZ << ", " << phiZ << ")" << G4endl;
#endif

  //Then the layer positions
  int i=0;
  readName(is,layMat);
  is >> layNum >> layRadius >> layAngle >> lengFront;
  for (i=0; i<5; i++) 
    is >> layPar[i];
#ifdef debug
  G4cout << tab << "Layer material: " << layMat << " Number " << layNum
       << " Radius " << layRadius << " Angle " << layAngle/deg 
       << " front dist " << lengFront << " Parameters ";
  for (i=0; i<5; i++)
    G4cout << layPar[i] << " ";
  G4cout << G4endl;
#endif

  //Then the crystal positions
  readName(is,crystMat);
  is >> crystNum >> crystLength >> crystTol;
  for (i=0; i<5; i++) 
    is >> crystPar[i];
#ifdef debug
  G4cout << tab << "Crystal material: " << crystMat << " Number " << crystNum
       << " Length " << crystLength << " Tolerance " << crystTol
       << " Parameters ";
  for (i=0; i<5; i++)
    G4cout << crystPar[i] << " ";
  G4cout << G4endl;
#endif

  //Then the support material
  readName(is,suppMat);
  is >> dxSupp >> dySupp >> dzSupp >> distSupp >> jump;
#ifdef debug
  G4cout << tab << "Support material: " << suppMat << " Dimensions " << dxSupp
       << ", " << dySupp << ", " << dzSupp << " Distance " << distSupp << G4endl;
#endif
   
  ///////////////////////////////////////////////////////////////
  // Close the file
  G4cout << " ==> Closing file " << File() << G4endl;
  is.close();

  return 1;

}

void CCalEcal::constructDaughters() {}
