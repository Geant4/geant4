// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Sort.cc,v 1.3 2000-08-28 08:57:59 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 source file
//
// G4Sort.cc
//
// ----------------------------------------------------------------------

#include "G4Sort.hh"


void sort_double( G4double v[], G4int left, G4int right )
{
  //  G4Sort elements in array from v[left] to v[right]  
  //  used recursively  
  //  algorithm comes from Kernighan and Ritchie, "The C Programming
  //  Language", second edition, p.87  
  
  G4int i, last;
  if ( left >= right )	// do nothing if array contains 
    return;		// fewer than two elements
  
  swap_double( v, left, ( left + right ) / 2 );   // move part. elt. 
  last = left;				          // to v[0] 

  for ( i = left+1; i <= right; i++ )	// partition 
    if ( v[i] < v[left] )
      swap_double( v, ++last, i );
  
  swap_double( v, left, last );	// restore partition element
  
  sort_double( v, left, last-1 );
  sort_double( v, last+1, right );
  return;
}


void swap_double( G4double v[], G4int i, G4int j )
{
  /*  interchange elements i and j in an array  */
  G4double temp;
  temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}
