// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Sort.cc,v 1.2 1999-12-15 14:50:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  File: G4Sort.cc
//  Author:  Alan Breakstone
//
//  Description
//
//	Routines to G4Sort arrays of various kinds of numbers
//

#include "G4Sort.hh"


void G4Sort_double( G4double v[], int left, int right )
{
  //  G4Sort elements in array from v[left] to v[right]  
  //  used recursively  
  //  algorithm comes from Kernighan and Ritchie, "The C Programming
  //  Language", second edition, p.87  
  
  int i, last;
  if ( left >= right )	// do nothing if array contains 
    return;		// fewer than two elements
  
  swap_double( v, left, ( left + right ) / 2 );   // move part. elt. 
  last = left;				          // to v[0] 

  for ( i = left+1; i <= right; i++ )	// partition 
    if ( v[i] < v[left] )
      swap_double( v, ++last, i );
  
  swap_double( v, left, last );	// restore partition element
  
  G4Sort_double( v, left, last-1 );
  G4Sort_double( v, last+1, right );
  return;
}


void swap_double( G4double v[], int i, int j )
{
  /*  interchange elements i and j in an array  */
  G4double temp;
  temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}




