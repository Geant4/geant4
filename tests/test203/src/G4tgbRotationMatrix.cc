#include "G4tgbRotationMatrix.hh"
#include "G4RotationMatrix.hh"
#include "G4tgrMessenger.hh"

using namespace CLHEP;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbRotationMatrix::G4tgbRotationMatrix( G4tgrRotationMatrix* tgr )
{
  theTgrRM = tgr;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4RotationMatrix* G4tgbRotationMatrix::BuildG4RotMatrix()
{
  std::vector<double> values = theTgrRM->GetValues();

  if( values.size() == 3 ) {
    return BuildG4RotMatrixFrom3( values );
  } else if( values.size() == 6 ) {
    return BuildG4RotMatrixFrom6( values );
  }else if( values.size() == 9 ) {
    return BuildG4RotMatrixFrom9( values );
  } else {
    G4cerr << " G4tgbRotationMatrix::BuildG4RotMatrix. Number of values is " << values.size() << G4endl;
    G4Exception(" It should be 3, 6, or 9 ");
  }

  return 0;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4RotationMatrix* G4tgbRotationMatrix::BuildG4RotMatrixFrom3( std::vector<double>& values )
{
  G4RotationMatrix* rotMat = new G4RotationMatrix();

  rotMat->rotateX( values[0] );
  rotMat->rotateY( values[1] );
  rotMat->rotateZ( values[2] );
  //  rotMat->invert();

  //-  cout << "  G4tgbRotationMatrix::buildG4RotMatrix(). unit rm " << *rotMat << endl;

  return rotMat;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4RotationMatrix* G4tgbRotationMatrix::BuildG4RotMatrixFrom6( std::vector<double>& values )
{
  double thetaX = values[0];
  double phiX = values[1];
  double thetaY = values[2];
  double phiY = values[3];
  double thetaZ = values[4];
  double phiZ = values[5];

  // build the 3 axis from the values
  Hep3Vector colx(sin(thetaX)*cos(phiX),sin(thetaX)*sin(phiX),cos(thetaX));
  Hep3Vector coly(sin(thetaY)*cos(phiY),sin(thetaY)*sin(phiY),cos(thetaY));
  Hep3Vector colz(sin(thetaZ)*cos(phiZ),sin(thetaZ)*sin(phiZ),cos(thetaZ));

  // Now create a G4RotationMatrix (HepRotation), which can be left handed. 
  // This is not forseen in CLHEP, but can be achieved using the
  // constructor which does not check its input arguments!   
  HepRep3x3 rottemp(colx.x(),coly.x(),colz.x(),
		    colx.y(),coly.y(),colz.y(),
		    colx.z(),coly.z(),colz.z()); //matrix representation (inverted)
  
  G4RotationMatrix* rotMat = new G4RotationMatrix(rottemp);

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
     cout << "  G4tgbRotationMatrix::buildG4RotMatrix(). " 
       << " colx " << colx << " = " << rotMat->colX() 
       << " coly " << coly << " = " << rotMat->colY() 
       << " colz " << colz << " = " << rotMat->colZ() 
       << endl;
#endif

  //-    cout << "  G4tgbRotationMatrix::buildG4RotMatrix(). unit rm " << *rotMat << endl;
  //- cout << " xrot " << xrot << " yrot " << yrot << " zrot " << zrot << endl;

  return rotMat;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4RotationMatrix* G4tgbRotationMatrix::BuildG4RotMatrixFrom9( std::vector<double>& values )
{
  // build the 3 axis from the values
  Hep3Vector colx(values[0],values[1],values[2]);
  Hep3Vector coly(values[3],values[4],values[5]);
  Hep3Vector colz(values[6],values[7],values[8]);

  G4cout << " BuildG4RotMatrixFrom9 colx " << colx << " y " << coly << " z " << colz << G4endl;
  // Now create a G4RotationMatrix (HepRotation), which can be left handed. 
  // This is not forseen in CLHEP, but can be achieved using the
  // constructor which does not check its input arguments!   
  HepRep3x3 rottemp(colx.x(),coly.x(),colz.x(),
		    colx.y(),coly.y(),colz.y(),
		    colx.z(),coly.z(),colz.z()); //matrix representation (inverted)
  
  G4RotationMatrix* rotMat = new G4RotationMatrix(rottemp);

  //-  cout << "  G4tgbRotationMatrix::buildG4RotMatrix(). unit rm " << *rotMat << endl;

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
     cout << "  G4tgbRotationMatrix::buildG4RotMatrixFrom9. " 
       << " colx " << colx << " = " << rotMat->colX() 
       << " coly " << coly << " = " << rotMat->colY() 
       << " colz " << colz << " = " << rotMat->colZ() 
       << endl;
#endif

  return rotMat;

}
