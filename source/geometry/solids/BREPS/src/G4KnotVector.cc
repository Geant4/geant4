// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KnotVector.cc,v 1.1 1999-01-07 16:07:43 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4KnotVector.hh"

G4KnotVector::G4KnotVector()
{
  knots=(G4double*)0;
}

G4KnotVector::G4KnotVector(const int sz)
{
  k_size=sz; 
  knots = new G4double[k_size]; 
  for(int a=0;a<k_size;a++)
    knots[a]=0;
}


G4KnotVector::~G4KnotVector() { delete [] knots;}  

//Copy constructor
G4KnotVector::G4KnotVector(const G4KnotVector& old_kv)
{
  k_size = old_kv.GetSize();
  knots = new G4double[k_size];
  for(register int a=0; a < old_kv.k_size; a++)knots[a]=old_kv.knots[a];
}

int G4KnotVector::GetKnotIndex(G4double k_value,const int order)
{
  int	   i, knot_index;
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


G4KnotVector* G4KnotVector::MultiplyKnotVector(const int num,
					       const G4double val)
{
  int n;
  G4double* knots_to_add;
  
  n            = CheckKnotVector( val );
  knots_to_add = new G4double[num-n];	
  
  for (int i = 0; i < num - n; i++)
    knots_to_add[i] = val;
  
  G4KnotVector* new_kv = new G4KnotVector();
  new_kv->k_size = num - n + GetSize();
  new_kv->knots  = MergeKnotVector(knots_to_add, num-n);
  
  delete [] knots_to_add;
  
  return new_kv;
}


G4double* G4KnotVector::MergeKnotVector(const G4double* knots_to_add, 
					const int add_size           )
{
  G4double *newknots;
  int kv1_ptr = 0, 
      kv2_ptr = 0, 
      newptr;
  int old_size = k_size;

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


int G4KnotVector::CheckKnotVector(const G4double val)
{
  int	num = 0;
  
  for ( int i = 0; i < k_size; i++) 
  {
    //	  if ( abs(val - knots[i]) < kCarTolerance)
    if ( val == knots[i] )
      num++;
  }

  return num;
}


void G4KnotVector::ExtractKnotVector( G4KnotVector* kv, 
				      const int upper,const  int lower)
{
  delete[] kv->knots;
  kv->k_size = upper-lower;
  kv->knots  = new G4double[kv->k_size];
  
  for ( int i = lower; i < upper; i++)
    kv->knots[i-lower] = knots[i];
}


