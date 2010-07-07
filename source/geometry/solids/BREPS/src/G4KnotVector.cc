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
// $Id: G4KnotVector.cc,v 1.12 2010-07-07 14:45:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4KnotVector.cc
//
// ----------------------------------------------------------------------

#include "G4KnotVector.hh"
#include "G4GeometryTolerance.hh"

G4KnotVector::G4KnotVector()
  : k_size(0), knots(0)
{
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

G4KnotVector::G4KnotVector(G4int sz)
{
  k_size=sz; 
  knots = new G4double[k_size]; 
  for(G4int a=0;a<k_size;a++)
    { knots[a]=0; }
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}


G4KnotVector::~G4KnotVector()
{
  delete [] knots;
}  

G4KnotVector::G4KnotVector(const G4KnotVector& orig)
{
  delete [] knots;
  k_size = orig.k_size;
  knots = new G4double[k_size];
  for(register G4int a=0; a < orig.k_size; a++)
    { knots[a] = orig.knots[a]; }
  kCarTolerance = orig.kCarTolerance;
}

G4KnotVector& G4KnotVector::operator=(const G4KnotVector& right)
{
  if (&right == this) return *this;
  delete [] knots;
  k_size = right.k_size;
  knots = new G4double[k_size];
  for(register G4int a=0; a < right.k_size; a++)
    knots[a] = right.knots[a];
  kCarTolerance = right.kCarTolerance;

  return *this;
}

G4int G4KnotVector::GetKnotIndex(G4double k_value, G4int order) const
{
  G4int	   i, knot_index;
  G4double knt;

  knt = knots[order - 1];
  if ( k_value < knt )
  {
    if (ApxEq( k_value, knt)) 
      k_value = knt;
    else
      return -1;
  }
  
  knt = knots[k_size - order + 1];
  if ( k_value > knt )
  {
    if (ApxEq( k_value, knt)) 
      k_value = knt;
    else
      return -1;
  }

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
    //	  if ( std::abs(val - knots[i]) < kCarTolerance)
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
