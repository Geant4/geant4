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
//
// $Id: G4tgbRotationMatrix.cc 68052 2013-03-13 14:38:53Z gcosmo $
//
//
// class G4tgbRotationMatrix

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbRotationMatrix.hh"
#include "G4RotationMatrix.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrUtils.hh"
#include "G4UIcommand.hh"

// -------------------------------------------------------------------------
G4tgbRotationMatrix::G4tgbRotationMatrix()
  : theTgrRM(0)
{
}


// -------------------------------------------------------------------------
G4tgbRotationMatrix::~G4tgbRotationMatrix()
{
}


// -------------------------------------------------------------------------
G4tgbRotationMatrix::G4tgbRotationMatrix( G4tgrRotationMatrix* tgr )
  : theTgrRM(tgr)
{
}


// -------------------------------------------------------------------------
G4RotationMatrix* G4tgbRotationMatrix::BuildG4RotMatrix()
{
  std::vector<G4double> values = theTgrRM->GetValues();

  if( values.size() == 3 ) {
    return BuildG4RotMatrixFrom3( values );
  } else if( values.size() == 6 ) {
    return BuildG4RotMatrixFrom6( values );
  } else if( values.size() == 9 ) {
    return BuildG4RotMatrixFrom9( values );
  }
  else
  {
    G4String ErrMessage = "Number of values is: "
                        + G4UIcommand::ConvertToString(G4int(values.size()))
                        + G4String(". It should be 3, 6, or 9 !");
    G4Exception("G4tgbRotationMatrix::BuildG4RotMatrix()",
                "InvalidData", FatalException, ErrMessage);
  }
  return 0;
}


// -------------------------------------------------------------------------
G4RotationMatrix*
G4tgbRotationMatrix::BuildG4RotMatrixFrom3( std::vector<G4double>& values )
{
  G4RotationMatrix* rotMat = new G4RotationMatrix();

  rotMat->rotateX( values[0] );
  rotMat->rotateY( values[1] );
  rotMat->rotateZ( values[2] );

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
    G4cout << " Constructing new G4RotationMatrix from 3 numbers "
           <<  GetName() << " : " << *rotMat << G4endl;
  }
#endif

  return rotMat;
}


// -------------------------------------------------------------------------
G4RotationMatrix*
G4tgbRotationMatrix::BuildG4RotMatrixFrom6( std::vector<G4double>& values )
{
  G4double thetaX = values[0];
  G4double phiX = values[1];
  G4double thetaY = values[2];
  G4double phiY = values[3];
  G4double thetaZ = values[4];
  G4double phiZ = values[5];

  // build the 3 axis from the values
  G4ThreeVector colx(std::sin(thetaX)*std::cos(phiX),
                     std::sin(thetaX)*std::sin(phiX),std::cos(thetaX));
  G4ThreeVector coly(std::sin(thetaY)*std::cos(phiY),
                     std::sin(thetaY)*std::sin(phiY),std::cos(thetaY));
  G4ThreeVector colz(std::sin(thetaZ)*std::cos(phiZ),
                     std::sin(thetaZ)*std::sin(phiZ),std::cos(thetaZ));

  // Now create a G4RotationMatrix (HepRotation), which can be left handed. 
  // This is not foreseen in CLHEP, but can be achieved using the
  // constructor which does not check its input arguments!   

  G4Rep3x3 rottemp(colx.x(),coly.x(),colz.x(), // matrix representation
                   colx.y(),coly.y(),colz.y(), // (inverted)
                   colx.z(),coly.z(),colz.z());
  
  G4RotationMatrix* rotMat = new G4RotationMatrix(rottemp);

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
    G4cout << " Constructing new G4RotationMatrix from 6 numbers "
           <<  GetName() << " : " << *rotMat << G4endl;
  }
#endif

  return rotMat;
}

// -------------------------------------------------------------------------
G4RotationMatrix*
G4tgbRotationMatrix::BuildG4RotMatrixFrom9( std::vector<G4double>& values )
{
  // build the 3 axis from the values
  G4ThreeVector colx(values[0],values[1],values[2]);
  G4ThreeVector coly(values[3],values[4],values[5]);
  G4ThreeVector colz(values[6],values[7],values[8]);

  // Now create a G4RotationMatrix (HepRotation), which can be left handed. 
  // This is not foreseen in CLHEP, but can be achieved using the
  // constructor which does not check its input arguments!

  G4Rep3x3 rottemp(colx.x(),coly.x(),colz.x(),  // matrix representation
                   colx.y(),coly.y(),colz.y(),  // (inverted)
                   colx.z(),coly.z(),colz.z());
  
  G4RotationMatrix* rotMat = new G4RotationMatrix(rottemp);

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
    G4cout << " Constructing new G4RotationMatrix from 9 numbers "
           <<  GetName() << " : " << *rotMat << G4endl;
  }
#endif

  return rotMat;
}
