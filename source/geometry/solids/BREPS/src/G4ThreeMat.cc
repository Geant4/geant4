// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ThreeMat.cc,v 1.4 2000-08-28 15:00:38 gcosmo Exp $
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


G4std::ostream& operator<<( G4std::ostream& os, const G4ThreeMat& m )
{
  // overwrite output operator << to Print out G4ThreeMat objects
  // using the PrintOn function defined below
  m.PrintOn( os );
  return os;
}


void G4ThreeMat::PrintOn( G4std::ostream& os ) const
{
  // printing function using C++ G4std::ostream class
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


void G4ThreeMat::operator=( const G4ThreeMat& m )
{ 
 //  assignment operator
  for ( G4int i = 0; i < 3 ; i++ )
  {
    row[i]    = m.row[i];
    column[i] = m.column[i];

    for ( G4int j = 0; j < 3 ; j++ ) 
      element[i][j] = m.element[i][j];
  }
}
      
		
G4ThreeMat G4ThreeMat::operator-()
{  
  //  unary - operator, change Sign of all elements
  G4double a[3][3];
  
  for ( G4int i = 0; i < 3 ; i++ )
  {
    for ( G4int j = 0; j < 3 ; j++ ) 
      a[i][j] = -element[i][j];
  }
  
  return G4ThreeMat( a );
}


G4ThreeMat G4ThreeMat::operator+=( const G4ThreeMat& m2 )
{ 
  //  overload += operator
  for ( G4int i = 0; i < 3 ; i++ )
  {
    this->row[i] += m2.row[i];
    this->column[i] += m2.column[i];
    
    for ( G4int j = 0; j < 3 ; j++ ) 
      this->element[i][j] += m2.element[i][j];
  }
  
  return *this;
}


G4ThreeMat G4ThreeMat::operator-=( const G4ThreeMat& m2 )
{
  //  overload -= operator
  for ( G4int i = 0; i < 3 ; i++ )
  {
    this->row[i] -= m2.row[i];
    this->column[i] -= m2.column[i];
    
    for ( G4int j = 0; j < 3 ; j++ ) 
      this->element[i][j] -= m2.element[i][j];
  }

  return *this;
}


G4ThreeMat operator+( const G4ThreeMat& m1, const G4ThreeMat& m2 )
{  
  //  overload binary + operator
  G4double a[3][3];

  for ( G4int i = 0; i < 3 ; i++ )
  {
    for ( G4int j = 0; j < 3 ; j++ )
      a[i][j] = m1.element[i][j] + m2.element[i][j];
	       
  }
	
  return G4ThreeMat( a );
}


G4ThreeMat operator-( const G4ThreeMat& m1, const G4ThreeMat& m2 )
{ 
  //  overload binary - operator
  G4double a[3][3];
  
  for ( G4int i = 0; i < 3 ; i++ )
  {
    for ( G4int j = 0; j < 3 ; j++ ) 
      a[i][j] = m1.element[i][j] - m2.element[i][j];  
  }
  
  return G4ThreeMat( a );
}


G4ThreeMat operator*( G4double x, const G4ThreeMat& m )
{
  //  overload binary * operator for constant Times matrix 
  G4double a[3][3];
  for ( G4int i = 0; i < 3 ; i++ ){
    for ( G4int j = 0; j < 3 ; j++ ) {
      a[i][j] = x * m.element[i][j]; 
    }
  }
  return G4ThreeMat( a );
}


G4Vector3D operator*( const G4ThreeMat& m, const G4Vector3D& v )
{
  //  overload binary * operator for matrix Times vector
  G4double a[3]; 
  for ( G4int i = 0; i < 3 ; i++ ){
    a[i] = m.row[i] * v; 
  }
  return G4Vector3D( a[0], a[1], a[2] );
}


G4ThreeMat operator*( const G4ThreeMat& m1, const G4ThreeMat& m2 )
{ 
 //  overload binary * operator for matrix Times matrix 
  G4double a[3][3];
  for ( G4int i = 0; i < 3 ; i++ ){
    for ( G4int j = 0; j < 3 ; j++ ) {
      a[i][j] = m1.row[i] * m2.column[j]; 
		}
  }
  return G4ThreeMat( a );
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
