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
// $Id: G4ThreeMat.cc,v 1.8 2006-06-29 18:42:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ThreeMat.cc
//
// ----------------------------------------------------------------------

#include "G4ThreeMat.hh"

G4ThreeMat::G4ThreeMat()
{  
  //  default (null) constructor 
  for ( G4int i = 0; i < 3 ; i++ )
  {
    row[i]    = G4Vector3D( 0., 0., 0. );
    column[i] = G4Vector3D( 0., 0., 0. );
    
    for ( G4int j = 0; j < 3 ; j++ ) 
      element[i][j] = 0.;
  }
}


G4ThreeMat::G4ThreeMat( G4double a[3][3] )
{
  //  constructor to make matrix from array
  for ( G4int i = 0; i < 3 ; i++ )
  {
    row[i]    = G4Vector3D( a[i][0], a[i][1], a[i][2] );
    column[i] = G4Vector3D( a[0][i], a[1][i], a[2][i] );
    
    for ( G4int j = 0; j < 3 ; j++ ) 
      element[i][j] = a[i][j];
  }
}


G4ThreeMat::~G4ThreeMat()
{
}


G4ThreeMat::G4ThreeMat( const G4ThreeMat& m )
{ 
  //  copy constructor
  for ( G4int i = 0; i < 3 ; i++ )
  {
    row[i]    = m.row[i];
    column[i] = m.column[i];
    
    for ( G4int j = 0; j < 3 ; j++ ) 
      element[i][j] = m.element[i][j];
  }
}


const char* G4ThreeMat::NameOf() const
{
  return "G4ThreeMat";
}


std::ostream& operator<<( std::ostream& os, const G4ThreeMat& m )
{
  // overwrite output operator << to Print out G4ThreeMat objects
  // using the PrintOn function defined below
  m.PrintOn( os );
  return os;
}


void G4ThreeMat::PrintOn( std::ostream& os ) const
{
  // printing function using C++ std::ostream class
  os << "[ " << element[0][0] << "\t" 
     << element[0][1] << "\t"
     << element[0][2] << "\n  "
     << element[1][0] << "\t"
     << element[1][1] << "\t"
     << element[1][2] << "\n  "
     << element[2][0] << "\t"
     << element[2][1] << "\t"
     << element[2][2] << " ]\n";
  /*
    for ( G4int i = 0; i < 3; i++ ) {
    os << "row   [" << i << "] " << row[i] << "\n"
    << "column[" << i << "] " << column[i] << "\n";
    }
    */
}


G4int G4ThreeMat::operator==( const G4ThreeMat& m ) const
{
  for ( G4int i = 0; i < 3 ; i++ ) 
  {
    for ( G4int j = 0; j < 3 ; j++ ) 
    {
      if ( element[i][j] != m.element[i][j] )
	return 0;
    }
  }

  return 1;
}


G4ThreeMat& G4ThreeMat::operator=( const G4ThreeMat& m )
{ 
  if (&m == this) return *this;
  for ( G4int i = 0; i < 3 ; i++ )
  {
    row[i]    = m.row[i];
    column[i] = m.column[i];

    for ( G4int j = 0; j < 3 ; j++ ) 
      element[i][j] = m.element[i][j];
  }
  return *this;
}


G4double G4ThreeMat::Determinant() const
{ 
  //  Determinant of a three by three matrix
  return element[0][0] * ( element[1][1] * element[2][2]
				    -  element[2][1] * element[1][2] )
    - element[0][1] * ( element[1][0] * element[2][2]
			-  element[2][0] * element[1][2] )
    + element[0][2] * ( element[1][0] * element[2][1]
			-  element[2][0] * element[1][1] );
}
