///////////////////////////////////////////////////////////////////////////////
// File: CrystalMatrix.cc
// Date:             09/99 S.B.
// Modifications: 06/09/99 S.B.
//                27/03/00 S.B. In OSCAR
///////////////////////////////////////////////////////////////////////////////
#include "CrystalMatrix.hh"

#include <fstream>
#include "utils.hh"

//#define debug

CrystalMatrix::~CrystalMatrix() {}

int CrystalMatrix::readFile() {
  ///////////////////////////////////////////////////////////////
  //Let's open the file
  cout << " ==> Opening file " << File() << " to read elements..."
       << endl;

  ifstream is;
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
  cout << tab << "General material: " << genMat << "\tBox dimensions "
       << widBox << ", " << lengBox << endl;
  cout << tab << "Positioned at (" << xpos << ", " << ypos << ", " << zpos
       << ") with rotation (" << thetaX << ", " << phiX << ", " << thetaY
       << ", " << phiY << ", " << thetaZ << ", " << phiZ << ")" << endl;
#endif

  //Then the layer positions
  int i=0;
  readName(is,layMat);
  is >> layNum >> layRadius >> layAngle >> lengFront;
  for (i=0; i<5; i++) 
    is >> layPar[i];
#ifdef debug
  cout << tab << "Layer material: " << layMat << " Number " << layNum
       << " Radius " << layRadius << " Angle " << layAngle/deg 
       << " front dist " << lengFront << " Parameters ";
  for (i=0; i<5; i++)
    cout << layPar[i] << " ";
  cout << endl;
#endif

  //Then the crystal positions
  readName(is,crystMat);
  is >> crystNum >> crystLength >> crystTol;
  for (i=0; i<5; i++) 
    is >> crystPar[i];
#ifdef debug
  cout << tab << "Crystal material: " << crystMat << " Number " << crystNum
       << " Length " << crystLength << " Tolerance " << crystTol
       << " Parameters ";
  for (i=0; i<5; i++)
    cout << crystPar[i] << " ";
  cout << endl;
#endif

  //Then the support material
  readName(is,suppMat);
  is >> dxSupp >> dySupp >> dzSupp >> distSupp >> jump;
#ifdef debug
  cout << tab << "Support material: " << suppMat << " Dimensions " << dxSupp
       << ", " << dySupp << ", " << dzSupp << " Distance " << distSupp << endl;
#endif
   
  ///////////////////////////////////////////////////////////////
  // Close the file
  cout << " ==> Closing file " << File() << endl;
  is.close();

  return 1;

}

void CrystalMatrix::constructDaughters() {}
