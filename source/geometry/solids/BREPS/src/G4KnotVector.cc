// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KnotVector.cc,v 1.4 2000-08-28 15:00:38 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4KnotVector.cc
//
// ----------------------------------------------------------------------

#include "G4KnotVector.hh"

G4KnotVector::G4KnotVector()
{
  knots=(G4double*)0;
}

G4KnotVector::G4KnotVector(G4int sz)
{
  k_size=sz; 
  knots = new G4double[k_size]; 
  for(G4int a=0;a<k_size;a++)
    knots[a]=0;
}


G4KnotVector::~G4KnotVector()
{
  delete [] knots;
}  

//Copy constructor
G4KnotVector::G4KnotVector(const G4KnotVector& old_kv)
{
  k_size = old_kv.GetSize();
  knots = new G4double[k_size];
  for(register G4int a=0; a < old_kv.k_size; a++)knots[a]=old_kv.knots[a];
}

G4int G4KnotVector::GetKnotIndex(G4double k_value, G4int order) const
{
  G4int	   i, knot_index;
  G4double knt;
 
  if ( k_value < ( knt = knots[order - 1])) 
    if (ApxEq( k_value, knt)) 
      k_value = knt;
    else
      return -1;
  
  if ( k_value > ( knt = knots[k_size - order + 1]))  
    if (ApxEq( k_value, knt)) 
      k_value = knt;
    else
      return -1; 

  if ( k_value == knots[k_size - order + 1] )
    knot_index = k_size - order - 1;
  else if ( k_value == knots[ order - 1] )
    knot_index = order - 1;
  else
  {
    knot_index = 0;
    
    for ( i = 0; i < k_size - 1; i++)
      if((knots[i]<k_value) && (k_value <= knots[i+1]))
	knot_index = i;  
  }

  return knot_index;
}


G4KnotVector* G4KnotVector::MultiplyKnotVector(G4int num,
					       G4double val)
{
  G4int n;
  G4double* knots_to_add;
  
  n            = CheckKnotVector( val );
  knots_to_add = new G4double[num-n];	
  
  for (G4int i = 0; i < num - n; i++)
    knots_to_add[i] = val;
  
  G4KnotVector* new_kv = new G4KnotVector();
  new_kv->k_size = num - n + GetSize();
  new_kv->knots  = MergeKnotVector(knots_to_add, num-n);
  
  delete [] knots_to_add;
  
  return new_kv;
}


G4double* G4KnotVector::MergeKnotVector(const G4double* knots_to_add, 
					G4int add_size           )
{
  G4double *newknots;
  G4int kv1_ptr = 0, 
        kv2_ptr = 0, 
        newptr;
  G4int old_size = k_size;

  newknots = new G4double[k_size + add_size];
    
  for ( newptr = 0; newptr < k_size+add_size; newptr++) 
    if ( kv1_ptr >= add_size )
      newknots[newptr] = knots[kv2_ptr++];
    else if ( kv2_ptr >= old_size )
      newknots[newptr] = knots_to_add[kv1_ptr++];
    else if ( knots_to_add[kv1_ptr] < knots[kv2_ptr])
      newknots[newptr] = knots_to_add[kv1_ptr++];
    else
      newknots[newptr] = knots[kv2_ptr++];
  
  return newknots;
}


G4int G4KnotVector::CheckKnotVector(G4double val) const
{
  G4int	num = 0;
  
  for ( G4int i = 0; i < k_size; i++) 
  {
    //	  if ( abs(val - knots[i]) < kCarTolerance)
    if ( val == knots[i] )
      num++;
  }

  return num;
}


void G4KnotVector::ExtractKnotVector( G4KnotVector* kv, 
				      G4int upper, G4int lower)
{
  delete[] kv->knots;
  kv->k_size = upper-lower;
  kv->knots  = new G4double[kv->k_size];
  
  for ( G4int i = lower; i < upper; i++)
    kv->knots[i-lower] = knots[i];
}
